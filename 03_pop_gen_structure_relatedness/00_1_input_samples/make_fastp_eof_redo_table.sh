###
# make_fastp_eof_redo_table.sh does this:
# Finds fastp outputs with unexpected EOF errors, extracts affected sample IDs,
# searches uploadku for the matching raw read pairs, and writes redo/missing tables.
###
#!/usr/bin/env bash
set -euo pipefail

# ---- paths ----
REPORT_DIR="/scratch/lt200308-agbsci/Quentin_project/00-samples/fastp/gzip_check_report_4139502_20260116_151011"
GZIP_TSV="${REPORT_DIR}/gzip_check_all.tsv"

FASTP_DIR="/scratch/lt200308-agbsci/Quentin_project/00-samples/fastp"
UPLOADKU_GLOB="/project/*agb*/uploadku"

OUT_TSV="${REPORT_DIR}/fastp_unexpected_eof_redo.tsv"
MISS_TSV="${REPORT_DIR}/fastp_unexpected_eof_notfound_in_uploadku.tsv"

# ---- sanity ----
[[ -s "$GZIP_TSV" ]] || { echo "[ERROR] Missing $GZIP_TSV" >&2; exit 1; }
[[ -d "$FASTP_DIR" ]] || { echo "[ERROR] Missing $FASTP_DIR" >&2; exit 1; }

# ---- collect uploadku dirs (expand glob safely) ----
mapfile -t UPDIRS < <(ls -d ${UPLOADKU_GLOB} 2>/dev/null || true)
if [[ "${#UPDIRS[@]}" -eq 0 ]]; then
  echo "[ERROR] No uploadku dirs found via: ${UPLOADKU_GLOB}" >&2
  exit 1
fi

echo "[INFO] Found uploadku dirs:"
printf "  - %s\n" "${UPDIRS[@]}"

# ---- helper: find raw R1/R2 by ID (CGAxxx.Exxxxx.L1) ----
find_raw_pair() {
  local id="$1"
  local r1="" r2=""

  for d in "${UPDIRS[@]}"; do
    r1=$(find "$d" -type f \
        \( -name "*${id}*R1*.f*q.gz" -o -name "*${id}*_R1*.f*q.gz" -o -name "*${id}*1.f*q.gz" -o -name "*${id}*_1.f*q.gz" \) \
        2>/dev/null | head -n 1 || true)
    r2=$(find "$d" -type f \
        \( -name "*${id}*R2*.f*q.gz" -o -name "*${id}*_R2*.f*q.gz" -o -name "*${id}*2.f*q.gz" -o -name "*${id}*_2.f*q.gz" \) \
        2>/dev/null | head -n 1 || true)

    if [[ -n "$r1" && -n "$r2" ]]; then
      echo -e "${r1}\t${r2}"
      return 0
    fi
  done

  echo -e "\t"
  return 0
}

# ---- extract EOF fastp files -> IDs ----
mapfile -t IDS < <(
  grep -F "unexpected end of file" "$GZIP_TSV" \
  | awk '{print $1}' \
  | xargs -n1 basename \
  | sed -E 's/\.R[12]\.fastp\.fq\.gz$//' \
  | sort -u
)

echo "[INFO] Unexpected-EOF IDs: ${#IDS[@]}"

# ---- write outputs ----
echo -e "id\tfastp_R1\tfastp_R2\traw_R1_uploadku\traw_R2_uploadku\tstatus" > "$OUT_TSV"
echo -e "id\tfastp_R1\tfastp_R2\tnote" > "$MISS_TSV"

ok=0
miss=0

for id in "${IDS[@]}"; do
  f1="${FASTP_DIR}/${id}.R1.fastp.fq.gz"
  f2="${FASTP_DIR}/${id}.R2.fastp.fq.gz"

  raw_pair="$(find_raw_pair "$id")"
  raw1="$(cut -f1 <<<"$raw_pair")"
  raw2="$(cut -f2 <<<"$raw_pair")"

  if [[ -n "$raw1" && -n "$raw2" ]]; then
    echo -e "${id}\t${f1}\t${f2}\t${raw1}\t${raw2}\tFOUND_RAW_PAIR" >> "$OUT_TSV"
    ok=$((ok+1))
  else
    echo -e "${id}\t${f1}\t${f2}\tRAW_NOT_FOUND_IN_UPLOADKU" >> "$OUT_TSV"
    echo -e "${id}\t${f1}\t${f2}\tRAW_NOT_FOUND_IN_UPLOADKU" >> "$MISS_TSV"
    miss=$((miss+1))
  fi
done

echo "[DONE] Wrote: $OUT_TSV"
echo "[DONE] Missing raw pairs: $MISS_TSV"
echo "[SUMMARY] FOUND_RAW_PAIR=$ok  RAW_NOT_FOUND=$miss"
