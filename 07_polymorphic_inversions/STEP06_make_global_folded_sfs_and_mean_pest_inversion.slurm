#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 124
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J realSFS_inv_folded
#SBATCH -o logs/realSFS_inv_folded.%j.out
#SBATCH -e logs/realSFS_inv_folded.%j.err

set -euo pipefail

# =============================================================================
# STEP06_make_global_folded_sfs_and_mean_pest_inversion.slurm
#
# Runs realSFS on the merged SAF to compute a folded bootstrap SFS,
# then averages bootstrap rows into a mean .pest file.
#
# Usage:
#   sbatch STEP06_make_global_folded_sfs_and_mean_pest_inversion.slurm
# =============================================================================

source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
OUTDIR="${BASE}/inversion_localpca/global_sfs"

SAF_IDX="${OUTDIR}/catfish.inversion.saf.idx"
OUT_SFS="${OUTDIR}/catfish.inversion.folded.sfs"
OUT_LOG="${OUTDIR}/realSFS.inversion.folded.log"
OUT_PEST="${OUTDIR}/catfish.inversion.folded.mean.pest"

THREADS=124
BOOTSTRAPS=100

mkdir -p "$OUTDIR"

[[ -f "$SAF_IDX" ]] || { echo "[ERROR] Missing SAF index: $SAF_IDX" >&2; exit 1; }

echo "[INFO] Running realSFS with $BOOTSTRAPS bootstraps (folded)"

realSFS "$SAF_IDX" -P "$THREADS" -bootstrap "$BOOTSTRAPS" -fold 1 \
  > "$OUT_SFS" \
  2> "$OUT_LOG"

# Average bootstrap rows to get mean pest prior
awk '{
  for(i=1; i<=NF; i++) s[i]+=$i
  n++
}
END{
  if(n==0){
    print "ERROR: no rows found in SFS file" > "/dev/stderr"
    exit 1
  }
  for(i=1; i<=NF; i++){
    printf "%.10f%s", s[i]/n, (i<NF ? OFS : ORS)
  }
}' OFS=" " "$OUT_SFS" > "$OUT_PEST"

echo "[INFO] Outputs:"
echo "  SFS  : $OUT_SFS"
echo "  PEST : $OUT_PEST"
