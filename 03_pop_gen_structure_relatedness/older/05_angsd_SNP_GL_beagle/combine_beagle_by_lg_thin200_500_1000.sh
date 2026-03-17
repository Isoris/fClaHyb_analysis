###
# combine_beagle_by_lg_thin200_500_1000.sh does this:
# Combines per-linkage-group beagle files into one whole-genome beagle file
# for each thinning level (200, 500, 1000 bp), with basic validation.
###
#!/usr/bin/env bash
set -euo pipefail

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/04_beagle_byRF_majmin"

for thin in 200 500 1000; do
  indir="${BASE}/thin_${thin}"
  [[ -d "$indir" ]] || { echo "[ERROR] missing dir: $indir" >&2; exit 1; }

  echo "[INFO] thin=$thin  dir=$indir"
  cd "$indir"

  OUT="../catfish.wholegenome.byRF.thin_${thin}.beagle.gz"
  LIST="beagle_thin_${thin}.list"
  TMP="catfish.wholegenome.byRF.thin_${thin}.beagle"

  # list LG files in order
  ls -1 catfish.C_gar_LG*.rf.txt.thin_${thin}.beagle.gz 2>/dev/null | sort -V > "$LIST"
  n=$(wc -l < "$LIST")
  [[ "$n" -gt 0 ]] || { echo "[ERROR] no LG beagles for thin=$thin in $indir" >&2; exit 1; }

  first=$(head -n1 "$LIST")

  # header once
  zcat "$first" | head -n1 > "$TMP"

  # append data lines from all files
  while read -r f; do
    zcat "$f" | tail -n +2 >> "$TMP"
  done < "$LIST"

  gzip -c "$TMP" > "$OUT"
  rm -f "$TMP"

  gzip -t "$OUT"
  echo "[OK] wrote $OUT  lines=$(zcat "$OUT" | wc -l)"
done
