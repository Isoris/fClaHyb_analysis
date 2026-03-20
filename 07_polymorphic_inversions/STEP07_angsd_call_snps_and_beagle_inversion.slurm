#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --ntasks=1
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J inv_beagle_step7
#SBATCH -o inv_beagle_step7.%A_%a.out
#SBATCH -e inv_beagle_step7.%A_%a.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

RFLIST="${1:?ERROR: provide path to chunk_rf.list as argument 1}"

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
OUTDIR="${BASE}/inversion_localpca"

REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
BAMLIST="${BASE}/bamlist.pp.samechr.tlenP99.filtered.txt"
PEST_SFS="${OUTDIR}/global_sfs/catfish.inversion.folded.mean.pest"
CALLABLE_SITES="${BASE}/fClaHyb_Gar_LG.mask_regions.inversion_acgt_allcase.1-based.angsd"

THREADS="${SLURM_CPUS_PER_TASK}"

MINQ=20
MINMAPQ=60
MININD=203
SETMAXDEPTH=2850
MAF=0.05
SNP_PVAL=1e-6

RF_FILE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$RFLIST")"
if [[ -z "${RF_FILE}" ]]; then
  echo "[ERROR] Empty RF_FILE for task ${SLURM_ARRAY_TASK_ID} from ${RFLIST}" >&2
  exit 1
fi

TAG="$(basename "$RF_FILE")"
TAG="${TAG//[^A-Za-z0-9_.-]/_}"

mkdir -p "${OUTDIR}"/{02_snps_beagle,03_sites,logs}

OUT_PREFIX="${OUTDIR}/02_snps_beagle/catfish.${TAG}"

angsd -b "$BAMLIST" -ref "$REF" \
  -GL 1 \
  -doMajorMinor 1 -doMaf 1 \
  -pest "$PEST_SFS" \
  -doGlf 2 \
  -doCounts 1 -doDepth 1 \
  -minQ "$MINQ" -minMapQ "$MINMAPQ" -baq 1 -C 50 \
  -setMaxDepth "$SETMAXDEPTH" \
  -minInd "$MININD" \
  -SNP_pval "$SNP_PVAL" -minMaf "$MAF" \
  -remove_bads 1 \
  -uniqueOnly 1 \
  -only_proper_pairs 1 \
  -rf "$RF_FILE" \
  -sites "$CALLABLE_SITES" \
  -P "$THREADS" \
  -out "$OUT_PREFIX" > "${OUTDIR}/logs/07_angsd_snps_beagle.${TAG}.log" 2>&1

zcat "${OUT_PREFIX}.mafs.gz" | awk 'NR==1{next} {print $1"\t"$2}' \
  > "${OUTDIR}/03_sites/sites.inversion.${TAG}.tsv"

echo "[DONE] Step 7 finished for TAG=$TAG"
echo "  OUT_PREFIX=${OUT_PREFIX}"
echo "  BEAGLE=${OUT_PREFIX}.beagle.gz"
echo "  MAFS=${OUT_PREFIX}.mafs.gz"
