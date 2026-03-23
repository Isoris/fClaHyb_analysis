#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J merge_inv_saf
#SBATCH -o logs/merge_inv_saf.%j.out
#SBATCH -e logs/merge_inv_saf.%j.err

set -euo pipefail

# =============================================================================
# STEP05_merge_chunk_saf_to_global_inversion.slurm
#
# Merges chunk SAF files into a single global SAF using realSFS cat.
#
# Usage:
#   sbatch STEP05_merge_chunk_saf_to_global_inversion.slurm
# =============================================================================

source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
OUTDIR="${BASE}/inversion_localpca"
CHUNKDIR="${OUTDIR}/chunks"
GLOBALDIR="${OUTDIR}/global_sfs"
THREADS="${SLURM_CPUS_PER_TASK}"

mkdir -p "$GLOBALDIR"

# Collect all chunk SAF index files
SAF_LIST="${GLOBALDIR}/saf.idx.list"
GLOBAL_PREFIX="${GLOBALDIR}/catfish.inversion"

find "$CHUNKDIR" -type f -name "*.saf.idx" | sort > "$SAF_LIST"

N=$(wc -l < "$SAF_LIST")
if [[ "$N" -eq 0 ]]; then
    echo "[ERROR] No .saf.idx files found under $CHUNKDIR" >&2
    exit 1
fi

echo "[INFO] Found $N chunk SAF files"

# Remove old global SAF if present
rm -f "${GLOBAL_PREFIX}.saf.idx" "${GLOBAL_PREFIX}.saf.gz" "${GLOBAL_PREFIX}.saf.pos.gz"

# Merge all chunk SAFs
realSFS cat -P "$THREADS" -outnames "$GLOBAL_PREFIX" $(cat "$SAF_LIST")

echo "[INFO] Done."
echo "[INFO] SAF list     : $SAF_LIST  ($N chunks)"
echo "[INFO] Global SAF   : ${GLOBAL_PREFIX}.saf.idx"
