#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=80
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J inv_saf_chunks
#SBATCH -o logs/inv_saf_chunks.%A_%a.out
#SBATCH -e logs/inv_saf_chunks.%A_%a.err
set -euo pipefail

# =============================================================================
# STEP04_run_angsd_saf_chunks_inversion.slurm
#
# Run ANGSD SAF per RF chunk under inversion-specific filters.
# Uses -doSaf 5 -doMajorMinor 2 to estimate SAF for folded SFS.
#
# Usage:
#   sbatch --array=0-$((N-1))%8 STEP04_run_angsd_saf_chunks_inversion.slurm chunk_rf.list
# =============================================================================

source ~/.bashrc
mamba activate assembly

CHUNKLIST="${1:?ERROR: provide chunk_rf.list as argv1}"

IDX="${SLURM_ARRAY_TASK_ID}"
LINE=$((IDX + 1))

RF_FILE=$(sed -n "${LINE}p" "$CHUNKLIST")
if [[ -z "${RF_FILE}" || ! -s "${RF_FILE}" ]]; then
  echo "[ERROR] Could not read RF_FILE from ${CHUNKLIST} at line ${LINE}" >&2
  exit 1
fi

TAG=$(basename "$RF_FILE" .rf.txt)

# ── Paths ───────────────────────────────────────────────────────────────────
BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
OUTDIR="${BASE}/inversion_localpca"
P="${SLURM_CPUS_PER_TASK}"

# Use the QC-passed bamlist from the het_roh workflow for consistency
BAMLIST="${BASE}/het_roh/01_inputs_check/bamlist_qcpass.txt"

# Inversion-callable sites (broader mask including soft-masked)
CALLABLE_SITES="${BASE}/fClaHyb_Gar_LG.mask_regions.inversion_acgt_allcase.1-based.angsd"

# ── Filters ─────────────────────────────────────────────────────────────────
MINQ=20
MINMAPQ=60
MININD=203
SETMAXDEPTH=2850

# ── Validate inputs ────────────────────────────────────────────────────────
for f in "$BAMLIST" "$REF" "$CALLABLE_SITES"; do
  [[ -s "$f" ]] || { echo "[ERROR] Missing required file: $f" >&2; exit 1; }
done

if [[ ! -s "${CALLABLE_SITES}.idx" || ! -s "${CALLABLE_SITES}.bin" ]]; then
  echo "[ERROR] Missing sites index: run 'angsd sites index ${CALLABLE_SITES}' first" >&2
  exit 1
fi

# ── Output dirs ─────────────────────────────────────────────────────────────
CHUNKDIR="${OUTDIR}/chunks/${TAG}"
mkdir -p "${CHUNKDIR}/01_saf" "${CHUNKDIR}/02_sfs" "${CHUNKDIR}/logs"

SAF_PREFIX="${CHUNKDIR}/01_saf/catfish.${TAG}"
SFS_OUT="${CHUNKDIR}/02_sfs/catfish.${TAG}.sfs"

# ── Run ANGSD SAF ──────────────────────────────────────────────────────────
echo "[INFO] Running ANGSD SAF for chunk: ${TAG}"

angsd -b "${BAMLIST}" -ref "${REF}" \
  -GL 1 \
  -doSaf 5 -doMajorMinor 2 \
  -doCounts 1 -doDepth 1 \
  -minQ "${MINQ}" -minMapQ "${MINMAPQ}" -baq 1 -C 50 \
  -setMaxDepth "${SETMAXDEPTH}" \
  -minInd "${MININD}" \
  -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
  -rf "${RF_FILE}" \
  -sites "${CALLABLE_SITES}" \
  -P "${P}" \
  -out "${SAF_PREFIX}" \
  2> "${CHUNKDIR}/logs/01_angsd_saf.${TAG}.log"

# ── Run realSFS bootstrap for this chunk ───────────────────────────────────
echo "[INFO] Running realSFS bootstrap for chunk: ${TAG}"

realSFS "${SAF_PREFIX}.saf.idx" -bootstrap 100 -fold 1 -P "${P}" \
  > "${SFS_OUT}" \
  2> "${CHUNKDIR}/logs/02_realSFS.${TAG}.log"

echo "[DONE] ${TAG}"
echo "  SAF prefix: ${SAF_PREFIX}"
echo "  SFS output: ${SFS_OUT}"
