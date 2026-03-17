###
# run_pcangsd_byLG_thin200_500_1000_K2_12.sh does this:
# Runs PCAngsd serially across all per-LG beagle files from the
# 200 bp, 500 bp, and 1000 bp thinned datasets, testing admixture K=2..12
# and writing outputs grouped by thinning, linkage group, and K.
###
#!/usr/bin/env bash
set -euo pipefail
source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin"
OUTBASE="${BASE}/05_pcangsd_byLG"
mkdir -p "${OUTBASE}"

THREADS=${THREADS:-32}

SAMPLES="${BASE}/list_of_samples_one_per_line_same_bamfile_list.tsv"
[[ -s "$SAMPLES" ]] || { echo "[ERROR] missing samples file: $SAMPLES" >&2; exit 1; }

L200="${OUTBASE}/beagle_LG_thin_200.list"
L500="${OUTBASE}/beagle_LG_thin_500.list"
L1000="${OUTBASE}/beagle_LG_thin_1000.list"
[[ -s "$L200" && -s "$L500" && -s "$L1000" ]] || {
  echo "[ERROR] missing list files: $L200 / $L500 / $L1000" >&2
  exit 1
}

run_one () {
  local thin="$1"
  local in="$2"
  local K="$3"

  [[ -s "$in" ]] || { echo "[ERROR] missing input: $in" >&2; exit 1; }

  local bn tag LG outdir out
  bn=$(basename "$in")
  tag="${bn%.beagle.gz}"

  LG=$(echo "$bn" | grep -oE 'LG[0-9]+' | head -n1)
  [[ -n "${LG:-}" ]] || LG="LGNA"

  outdir="${OUTBASE}/thin_${thin}/${LG}/K${K}"
  mkdir -p "$outdir"
  out="${outdir}/${tag}.pcangsd"

  echo "[INFO] thin=$thin LG=$LG K=$K"
  echo "[INFO] IN=$in"
  echo "[INFO] OUT=$out"

  pcangsd \
    --beagle "$in" \
    --eig 6 \
    --threads "$THREADS" \
    --maf 0.05 \
    --admix \
    --admix-K "$K" \
    --admix-seed 1 \
    --tree \
    --tree-samples "$SAMPLES" \
    --out "$out"
}

for thin in 200 500 1000; do
  list_var="L${thin}"
  LIST="${!list_var}"

  while read -r IN; do
    [[ -z "$IN" ]] && continue
    for K in $(seq 2 12); do
      run_one "$thin" "$BASE/$IN" "$K"
      # NOTE: your list files contain paths relative to BASE (04_beagle...).
      # If your list already contains absolute paths, replace "$BASE/$IN" with "$IN".
    done
  done < "$LIST"
done
