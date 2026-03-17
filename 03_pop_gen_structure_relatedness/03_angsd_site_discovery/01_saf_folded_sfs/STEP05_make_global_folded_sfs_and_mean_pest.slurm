#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 124
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J realSFS_folded
#SBATCH -o realSFS_folded.%j.out
#SBATCH -e realSFS_folded.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

###
# STEP05_make_global_folded_sfs_and_mean_pest.slurm
#
# Does this:
# 1. Runs realSFS on the merged global SAF to compute a folded bootstrap SFS
# 2. Averages bootstrap rows into one mean prior vector (.pest)
# 3. Prints quick QC
###

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
OUTDIR="${BASE}/popstruct_global/global_sfs"

SAF_IDX="${OUTDIR}/catfish.global.saf.idx"
OUT_SFS="${OUTDIR}/catfish.global.folded.sfs"
OUT_LOG="${OUTDIR}/realSFS.global.folded.log"
OUT_PEST="${OUTDIR}/catfish.global.folded.mean.pest"

THREADS=124
BOOTSTRAPS=100

mkdir -p "$OUTDIR"

[[ -f "$SAF_IDX" ]] || { echo "[ERROR] Missing SAF index: $SAF_IDX" >&2; exit 1; }

realSFS "$SAF_IDX" -P "$THREADS" -bootstrap "$BOOTSTRAPS" -fold 1 \
  > "$OUT_SFS" \
  2> "$OUT_LOG"

awk '{
  for(i=1;i<=NF;i++) s[i]+=$i
  n++
}
END{
  if(n==0){
    print "ERROR: no rows found in SFS file" > "/dev/stderr"
    exit 1
  }
  for(i=1;i<=NF;i++){
    printf "%.10f%s", s[i]/n, (i<NF?OFS:ORS)
  }
}' OFS=" " "$OUT_SFS" > "$OUT_PEST"

echo "[INFO] QC on folded bootstrap rows:"
awk '{
  sum=0
  for(i=1;i<=NF;i++) sum+=$i
  printf "row=%d\tNF=%d\tsum=%.10f\n", NR, NF, sum
}' "$OUT_SFS" | head -n 20

echo "[INFO] QC on mean pest:"
awk '{
  sum=0
  for(i=1;i<=NF;i++) sum+=$i
  printf "file=%s\tNF=%d\tsum=%.10f\n", FILENAME, NF, sum
}' "$OUT_PEST"

echo "[INFO] Outputs:"
echo "  SFS  : $OUT_SFS"
echo "  LOG  : $OUT_LOG"
echo "  PEST : $OUT_PEST"
