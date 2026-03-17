###
# 03_angsd_call_snps_step3_v2.slurm does this:
# Runs ANGSD per RF chunk to call candidate SNP sites using the global folded SFS
# as a prior, applies site-level filters, and writes per-chunk .mafs outputs
# for downstream site merging and thinning.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --ntasks=1
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J angsd_step3
#SBATCH -o angsd_step3.%A_%a.out
#SBATCH -e angsd_step3.%A_%a.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

###############################################################################
# Step 3: Call candidate SNP sites per chromosome/region (-rf) using global SFS
# Usage:
#   sbatch --array=0-$((N-1))%8 03_angsd_call_snps_step3.slurm "${OUTDIR}/chunk_rf.list"
###############################################################################

# ---------------------------- INPUT ARG --------------------------------------
RFLIST="${1:?ERROR: provide path to chunk_rf.list as argument 1}"

# ---------------------------- PATHS / INPUTS ---------------------------------
BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
OUTDIR="${BASE}/popstruct_global"

REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
BAMLIST="${BASE}/bamlist.pp.samechr.tlenP99.filtered.txt"

# Global folded SFS (the one you already computed by concatenating SAFs)
#PEST_SFS="${BASE}/popstruct_global/global_sfs/catfish.global.folded.sfs"
PEST_SFS="${BASE}/popstruct_global/global_sfs/catfish.global.folded.sfs.mean.pest"

# Callable mask in ANGSD native format (1-based) + indexed
CALLABLE_SITES="${BASE}/fClaHyb_Gar_LG.mask_regions.normalACGT.renamed.1-based.angsd"

THREADS="${SLURM_CPUS_PER_TASK}"

# ---------------------------- FILTER SETTINGS --------------------------------
MINQ=25
MINMAPQ=25
MINDEPTHIND=3
MAXDEPTHIND=57
MININD=200

MAF=0.05
SNP_PVAL=1e-6

REMOVE_FLAGS=1
UNIQUE_ONLY=1
ONLY_PROPER_PAIRS=1

# SNP-stat filters
#SB_PVAL=1e-4
#EDGE_PVAL=1e-4
#QSCORE_PVAL=1e-4

# ---------------------------- DERIVE RF + TAG --------------------------------
RF_FILE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$RFLIST")"
if [[ -z "${RF_FILE}" ]]; then
  echo "[ERROR] Empty RF_FILE for task ${SLURM_ARRAY_TASK_ID} from ${RFLIST}" >&2
  exit 1
fi

# Make a safe tag for filenames (Chr1, LG01, etc.)
TAG="$(basename "$RF_FILE")"
TAG="${TAG//[^A-Za-z0-9_.-]/_}"

# ---------------------------- SETUP ------------------------------------------
mkdir -p "${OUTDIR}"/{02_snps,03_sites,logs}

MININD_ARG=()
if [[ "$MININD" -gt 0 ]]; then
  MININD_ARG=(-minInd "$MININD")
fi

echo "[INFO] TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "[INFO] RF_FILE=$RF_FILE"
echo "[INFO] TAG=$TAG"
echo "[INFO] PEST_SFS=$PEST_SFS"
echo "[INFO] CALLABLE_SITES=$CALLABLE_SITES"

# ---------------------- 3) Call candidate SNP sites --------------------------
SNP_PREFIX="${OUTDIR}/02_snps/catfish.${TAG}"

#  -doSnpStat 1 \
#  -sb_pval "$SB_PVAL" \
#  -edge_pval "$EDGE_PVAL" \
#  -qscore_pval "$QSCORE_PVAL" \
set -x

angsd -b "$BAMLIST" -ref "$REF" \
  -GL 1 \
  -doMajorMinor 1 -doMaf 1 \
  -pest "$PEST_SFS" \
  -doCounts 1 -doDepth 1 \
  -minQ "$MINQ" -minMapQ "$MINMAPQ" -baq 1 -C 50 \
  -setMinDepthInd "$MINDEPTHIND" -setMaxDepthInd "$MAXDEPTHIND" \
  "${MININD_ARG[@]}" \
  -SNP_pval "$SNP_PVAL" -minMaf "$MAF" \
  -remove_bads "$REMOVE_FLAGS" \
  -uniqueOnly "$UNIQUE_ONLY" \
  -only_proper_pairs "$ONLY_PROPER_PAIRS" \
  -rf "$RF_FILE" \
  -sites "$CALLABLE_SITES" \
  -P "$THREADS" \
  -out "$SNP_PREFIX" > "${OUTDIR}/logs/03_angsd_snps.${TAG}.log" 2>&1

set +x
#  -out "$SNP_PREFIX" 2> "${OUTDIR}/logs/03_angsd_snps.${TAG}.log"

# Extract sites (chr pos) from .mafs.gz
#zcat "${SNP_PREFIX}.mafs.gz" | awk 'NR==1{next} {print $1"\t"$2}' \
#  > "${OUTDIR}/03_sites/sites.raw.${TAG}.tsv"

echo "[DONE] Step 3 finished for TAG=$TAG"
echo "  SNP_PREFIX=${SNP_PREFIX}"
echo "  RAW_SITES=${OUTDIR}/03_sites/sites.raw.${TAG}.tsv"
