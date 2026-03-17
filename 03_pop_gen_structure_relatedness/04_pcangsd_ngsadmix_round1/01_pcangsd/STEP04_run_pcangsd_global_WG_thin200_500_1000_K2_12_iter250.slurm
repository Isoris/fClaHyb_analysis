#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=237G
#SBATCH -t 2-00:00:00
#SBATCH -J pcangsdWG
#SBATCH -o 06_pcangsdWG.%j.out
#SBATCH -e 06_pcangsdWG.%j.err

set -euo pipefail

source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin"

SAMPLES="${BASE}/list_of_samples_one_per_line_same_bamfile_list.tsv"
[[ -s "$SAMPLES" ]] || { echo "[ERROR] missing samples file: $SAMPLES" >&2; exit 1; }

BEAGLEDIR="${BASE}/04_beagle_byRF_majmin"
declare -A BEAGLE=(
  [200]="${BEAGLEDIR}/catfish.wholegenome.byRF.thin_200.beagle.gz"
  [500]="${BEAGLEDIR}/catfish.wholegenome.byRF.thin_500.beagle.gz"
  [1000]="${BEAGLEDIR}/catfish.wholegenome.byRF.thin_1000.beagle.gz"
)

OUTBASE="${BASE}/05_pcangsd_global_WGthin"
mkdir -p "${OUTBASE}/logs"

THINS=(200 500 1000)

THREADS="${THREADS:-64}"
MAF="${MAF:-0.05}"
EIG="${EIG:-4}"


for thin in "${THINS[@]}"; do
  IN="${BEAGLE[$thin]}"
  [[ -s "$IN" ]] || { echo "[ERROR] missing beagle: $IN" >&2; exit 1; }
  gzip -t "$IN" || { echo "[ERROR] corrupt gzip: $IN" >&2; exit 2; }

  tag="catfish.wg.byRF.thin_${thin}"

  for K in $(seq 2 12); do
    OUTDIR="${OUTBASE}/iter_250_thin_${thin}/K${K}"
    mkdir -p "$OUTDIR"
    OUT="${OUTDIR}/${tag}.pcangsd"

    echo "============================================================"
    echo "[INFO] thin=$thin  K=$K"
    echo "[INFO] IN=$IN"
    echo "[INFO] OUT=$OUT"
    echo "[INFO] threads=$THREADS maf=$MAF eig=$EIG"

    pcangsd \
      -b "$IN" \
      -e "$EIG" \
      -t "$THREADS" \
      --maf "$MAF" \
      --iter 250 \
      --tree \
      --tree-samples "$SAMPLES" \
      -o "$OUT" \
      > "${OUTBASE}/logs/pcangsd_${tag}_K${K}.log" 2>&1
  done
done
