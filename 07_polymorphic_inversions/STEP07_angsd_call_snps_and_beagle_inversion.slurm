#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --ntasks=1
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J inv_beagle_step7
#SBATCH -o logs/inv_beagle_step7.%A_%a.out
#SBATCH -e logs/inv_beagle_step7.%A_%a.err

set -euo pipefail

# =============================================================================
# STEP07_angsd_call_snps_and_beagle_inversion.slurm
#
# Call biallelic SNPs using the inversion-specific folded SFS prior via -pest,
# and export BEAGLE genotype likelihoods.
#
# IMPORTANT: -doSaf 1 is required here so that -pest is used as the SFS prior
# for -doMaf. Without -doSaf, ANGSD ignores -pest entirely.
#
# Usage:
#   sbatch --array=0-$((N-1))%8 STEP07_angsd_call_snps_and_beagle_inversion.slurm chunk_rf.list
# =============================================================================

source ~/.bashrc
mamba activate assembly

RFLIST="${1:?ERROR: provide path to chunk_rf.list as argument 1}"

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
OUTDIR="${BASE}/inversion_localpca"

REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
BAMLIST="${BASE}/het_roh/01_inputs_check/bamlist_qcpass.txt"
PEST_SFS="${OUTDIR}/global_sfs/catfish.inversion.folded.mean.pest"
CALLABLE_SITES="${BASE}/fClaHyb_Gar_LG.mask_regions.inversion_acgt_allcase.1-based.angsd"

THREADS="${SLURM_CPUS_PER_TASK}"

# ── Filters ─────────────────────────────────────────────────────────────────
MINQ=20
MINMAPQ=60
MININD=203
SETMAXDEPTH=2850
MAF=0.05
SNP_PVAL=1e-6

# ── Parse chunk RF file ────────────────────────────────────────────────────
RF_FILE="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$RFLIST")"
if [[ -z "${RF_FILE}" || ! -s "${RF_FILE}" ]]; then
  echo "[ERROR] Empty or missing RF_FILE for task ${SLURM_ARRAY_TASK_ID} from ${RFLIST}" >&2
  exit 1
fi

# Extract TAG the same way as STEP04: strip .rf.txt suffix from basename
TAG="$(basename "$RF_FILE" .rf.txt)"

# ── Validate inputs ────────────────────────────────────────────────────────
for f in "$BAMLIST" "$REF" "$PEST_SFS" "$CALLABLE_SITES"; do
  [[ -s "$f" ]] || { echo "[ERROR] Missing required file: $f" >&2; exit 1; }
done

# ── Output dirs ─────────────────────────────────────────────────────────────
mkdir -p "${OUTDIR}/02_snps_beagle" "${OUTDIR}/03_sites" "${OUTDIR}/logs"

OUT_PREFIX="${OUTDIR}/02_snps_beagle/catfish.${TAG}"

echo "[INFO] Running ANGSD SNP calling + BEAGLE for chunk: ${TAG}"

# NOTE on -doSaf 1:
# ANGSD requires -doSaf to activate the SFS prior specified by -pest.
# We use -doSaf 1 here (not -doSaf 5 as in STEP04) because in this step
# we need the prior for MAF estimation, not for building a new SFS.
# The -pest file provides the folded SFS prior from STEP06.
angsd -b "$BAMLIST" -ref "$REF" \
  -GL 1 \
  -doSaf 1 \
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

# Extract site list from mafs
zcat "${OUT_PREFIX}.mafs.gz" | awk 'NR>1 {print $1"\t"$2}' \
  > "${OUTDIR}/03_sites/sites.inversion.${TAG}.tsv"

echo "[DONE] Step 7 finished for TAG=${TAG}"
echo "  BEAGLE : ${OUT_PREFIX}.beagle.gz"
echo "  MAFS   : ${OUT_PREFIX}.mafs.gz"
echo "  Sites  : ${OUTDIR}/03_sites/sites.inversion.${TAG}.tsv"
