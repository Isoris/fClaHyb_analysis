#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J merge_inv_saf
#SBATCH -o merge_inv_saf.%j.out
#SBATCH -e merge_inv_saf.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
OUTDIR="${BASE}/inversion_localpca"
CHUNKDIR="${OUTDIR}/chunks"
GLOBALDIR="${OUTDIR}/global_sfs"
THREADS="${SLURM_CPUS_PER_TASK}"

mkdir -p "$GLOBALDIR"

SAF_LIST="${GLOBALDIR}/saf.idx.list"
GLOBAL_PREFIX="${GLOBALDIR}/catfish.inversion"

find "$CHUNKDIR" -type f -name "*.saf.idx" | sort > "$SAF_LIST"

N=$(wc -l < "$SAF_LIST")
if [[ "$N" -eq 0 ]]; then
    echo "[ERROR] No .saf.idx files found under $CHUNKDIR" >&2
    exit 1
fi

rm -f "${GLOBAL_PREFIX}.saf.idx" "${GLOBAL_PREFIX}.saf.gz" "${GLOBAL_PREFIX}.saf.pos.gz"

realSFS cat -P "$THREADS" -outnames "$GLOBAL_PREFIX" $(cat "$SAF_LIST")

echo "[INFO] Done."
echo "[INFO] SAF list     : $SAF_LIST"
echo "[INFO] Global SAF   : ${GLOBAL_PREFIX}.saf.idx"
