###
# assign_species_kmers_mash_meryl.sh does this:
# Builds Mash sketches and/or meryl k-mer sets for two reference species,
# scores each WGS sample, and assigns it to Species1, Species2, or Ambiguous.
###
#!/usr/bin/env bash
set -euo pipefail

########################
# Defaults
########################
REF_A=""
REF_B=""
READS_DIR="."
OUTDIR="species_assign_out"
MODE="both"          # mash | meryl | both
K=21
THREADS=16
JOBS=1
MASH_S=10000
MASH_MIN_COPIES=2
MIN_WGS_COUNT=2
CALL_LOG2_THR=1.0
CALL_MASH_MARGIN=0.01
MAX_REF_COUNT=0

########################
# Args
########################
usage() {
  cat <<EOF
Usage:
  $0 -a Species1.fa -b Species2.fa -r /path/to/reads -o outdir [options]

Required:
  -a  reference FASTA for Species1
  -b  reference FASTA for Species2
  -r  reads directory containing *.fq.gz/*.fastq.gz (recursive)
Optional:
  -o  output dir (default: $OUTDIR)
  --mode mash|meryl|both (default: $MODE)
  -k  k-mer size (default: $K)
  -t  threads (default: $THREADS)
  -j  jobs (currently only used for single-method modes) (default: $JOBS)
  --mash-s N (default: $MASH_S)
  --mash-m N (default: $MASH_MIN_COPIES)
  --min-wgs-count N (default: $MIN_WGS_COUNT)
  --call-log2 THR (default: $CALL_LOG2_THR)
  --call-mash-margin M (default: $CALL_MASH_MARGIN)
  --max-ref-count N (default: $MAX_REF_COUNT)

Outputs:
  \$OUTDIR/results.tsv
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a) REF_A="$2"; shift 2;;
    -b) REF_B="$2"; shift 2;;
    -r) READS_DIR="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    --mode) MODE="$2"; shift 2;;
    -k) K="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    -j) JOBS="$2"; shift 2;;
    --mash-s) MASH_S="$2"; shift 2;;
    --mash-m) MASH_MIN_COPIES="$2"; shift 2;;
    --min-wgs-count) MIN_WGS_COUNT="$2"; shift 2;;
    --call-log2) CALL_LOG2_THR="$2"; shift 2;;
    --call-mash-margin) CALL_MASH_MARGIN="$2"; shift 2;;
    --max-ref-count) MAX_REF_COUNT="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

if [[ -z "$REF_A" || -z "$REF_B" ]]; then
  echo "ERROR: -a and -b are required."
  usage
  exit 1
fi

mkdir -p "$OUTDIR"/{00_lists,01_refs,02_samples,logs,tmp}

need_tool() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing tool '$1' in PATH"; exit 1; }
}

if [[ "$MODE" == "mash" || "$MODE" == "both" ]]; then
  need_tool mash
fi
if [[ "$MODE" == "meryl" || "$MODE" == "both" ]]; then
  need_tool meryl
fi

echo "[1/5] Finding FASTQ files under: $READS_DIR"
mapfile -t ALLFQ < <(find "$READS_DIR" -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fastq" \) | sort)
if [[ ${#ALLFQ[@]} -eq 0 ]]; then
  echo "ERROR: no fastq files found in $READS_DIR"
  exit 1
fi

declare -A SAMPLE_FILES
declare -a SAMPLE_NAMES

normalize_sample_name() {
  local bn="$1"
  bn="${bn%.fastq.gz}"; bn="${bn%.fq.gz}"; bn="${bn%.fastq}"; bn="${bn%.fq}"
  bn="$(echo "$bn" | sed -E 's/([._-]R?[12])$//; s/([._-][12])$//')"
  echo "$bn"
}

for f in "${ALLFQ[@]}"; do
  bn="$(basename "$f")"
  s="$(normalize_sample_name "$bn")"
  if [[ -z "${SAMPLE_FILES[$s]+x}" ]]; then
    SAMPLE_FILES["$s"]="$f"
    SAMPLE_NAMES+=("$s")
  else
    SAMPLE_FILES["$s"]+=$'\n'"$f"
  fi
done

printf "%s\n" "${SAMPLE_NAMES[@]}" | sort > "$OUTDIR/00_lists/samples.txt"
echo "  Found ${#SAMPLE_NAMES[@]} samples."

PY_SCORE="$OUTDIR/tmp/score_call.py"
cat > "$PY_SCORE" <<'PY'
import math, sys
A = float(sys.argv[1]); B = float(sys.argv[2]); thr = float(sys.argv[3])
log2 = math.log((A+1.0)/(B+1.0), 2)
call = "Ambiguous"
if log2 > thr: call = "Species1"
elif log2 < -thr: call = "Species2"
print(f"{log2:.6f}\t{call}")
PY

REF_A_MSH="$OUTDIR/01_refs/Species1.msh"
REF_B_MSH="$OUTDIR/01_refs/Species2.msh"

do_mash_refs() {
  echo "[2/5] (Mash) Sketching references"
  mash sketch -k "$K" -s "$MASH_S" -p "$THREADS" -o "$OUTDIR/01_refs/Species1" "$REF_A" >"$OUTDIR/logs/mash_refA.log" 2>&1
  mash sketch -k "$K" -s "$MASH_S" -p "$THREADS" -o "$OUTDIR/01_refs/Species2" "$REF_B" >"$OUTDIR/logs/mash_refB.log" 2>&1
}

mash_one_sample() {
  local sample="$1"
  local listfile="$OUTDIR/02_samples/${sample}.files.txt"
  printf "%s\n" "${SAMPLE_FILES[$sample]}" > "$listfile"

  local msh="$OUTDIR/02_samples/${sample}.msh"
  mash sketch -r -k "$K" -s "$MASH_S" -m "$MASH_MIN_COPIES" -p "$THREADS" -o "$OUTDIR/02_samples/${sample}" -l "$listfile" \
    >"$OUTDIR/logs/mash_${sample}.log" 2>&1

  local distA distB
  distA="$(mash dist "$REF_A_MSH" "$msh" | awk 'NR==1{print $3}')"
  distB="$(mash dist "$REF_B_MSH" "$msh" | awk 'NR==1{print $3}')"

  local call="Ambiguous"
  call="$(awk -v a="$distA" -v b="$distB" -v m="$CALL_MASH_MARGIN" 'BEGIN{
    if (a + m < b) print "Species1";
    else if (b + m < a) print "Species2";
    else print "Ambiguous";
  }')"

  printf "%s\t%s\t%s\n" "$distA" "$distB" "$call" > "$OUTDIR/tmp/${sample}.mash.tsv"
}

REF_A_DB="$OUTDIR/01_refs/Species1.k${K}.meryl"
REF_B_DB="$OUTDIR/01_refs/Species2.k${K}.meryl"
A_ONLY_DB="$OUTDIR/01_refs/Species1_only.k${K}.meryl"
B_ONLY_DB="$OUTDIR/01_refs/Species2_only.k${K}.meryl"

meryl_stat_total() {
  local db="$1"
  meryl statistics "$db" | python3 - <<'PY'
import re, sys
txt = sys.stdin.read()
patterns = [
  r'Number of k-mers:\s*([0-9]+)',
  r'Total k-mers:\s*([0-9]+)',
  r'Total number of k-mers:\s*([0-9]+)',
  r'k-mers:\s*([0-9]+)\s*$',
]
for pat in patterns:
  m = re.search(pat, txt, re.MULTILINE)
  if m:
    print(m.group(1))
    break
else:
  print("0")
PY
}

do_meryl_refs() {
  echo "[2/5] (meryl) Counting reference k-mers"
  meryl count k="$K" threads="$THREADS" output "$REF_A_DB" "$REF_A" >"$OUTDIR/logs/meryl_refA.log" 2>&1
  meryl count k="$K" threads="$THREADS" output "$REF_B_DB" "$REF_B" >"$OUTDIR/logs/meryl_refB.log" 2>&1

  if [[ "$MAX_REF_COUNT" -gt 0 ]]; then
    echo "  (meryl) Filtering repetitive reference k-mers (max=$MAX_REF_COUNT)"
    meryl filter max="$MAX_REF_COUNT" output "$OUTDIR/01_refs/Species1.filtered.meryl" "$REF_A_DB" >"$OUTDIR/logs/meryl_refA_filter.log" 2>&1
    meryl filter max="$MAX_REF_COUNT" output "$OUTDIR/01_refs/Species2.filtered.meryl" "$REF_B_DB" >"$OUTDIR/logs/meryl_refB_filter.log" 2>&1
    REF_A_DB="$OUTDIR/01_refs/Species1.filtered.meryl"
    REF_B_DB="$OUTDIR/01_refs/Species2.filtered.meryl"
  fi

  echo "  (meryl) Building species-specific sets via difference"
  meryl difference output "$A_ONLY_DB" "$REF_A_DB" "$REF_B_DB" >"$OUTDIR/logs/meryl_Aonly.log" 2>&1
  meryl difference output "$B_ONLY_DB" "$REF_B_DB" "$REF_A_DB" >"$OUTDIR/logs/meryl_Bonly.log" 2>&1
}

meryl_one_sample() {
  local sample="$1"
  local sample_db="$OUTDIR/02_samples/${sample}.k${K}.meryl"
  local filt_db="$OUTDIR/02_samples/${sample}.k${K}.min${MIN_WGS_COUNT}.meryl"
  local a_hit_db="$OUTDIR/02_samples/${sample}.Aonly.hit.meryl"
  local b_hit_db="$OUTDIR/02_samples/${sample}.Bonly.hit.meryl"

  local listfile="$OUTDIR/02_samples/${sample}.files.txt"
  printf "%s\n" "${SAMPLE_FILES[$sample]}" > "$listfile"

  meryl count k="$K" threads="$THREADS" output "$sample_db" $(cat "$listfile") \
    >"$OUTDIR/logs/meryl_${sample}_count.log" 2>&1

  if [[ "$MIN_WGS_COUNT" -gt 1 ]]; then
    meryl filter min="$MIN_WGS_COUNT" output "$filt_db" "$sample_db" \
      >"$OUTDIR/logs/meryl_${sample}_filter.log" 2>&1
  else
    filt_db="$sample_db"
  fi

  meryl intersect output "$a_hit_db" "$filt_db" "$A_ONLY_DB" >"$OUTDIR/logs/meryl_${sample}_Aint.log" 2>&1
  meryl intersect output "$b_hit_db" "$filt_db" "$B_ONLY_DB" >"$OUTDIR/logs/meryl_${sample}_Bint.log" 2>&1

  local A B
  A="$(meryl_stat_total "$a_hit_db")"
  B="$(meryl_stat_total "$b_hit_db")"

  local scored log2ratio call
  scored="$(python3 "$PY_SCORE" "$A" "$B" "$CALL_LOG2_THR")"
  log2ratio="$(echo "$scored" | awk '{print $1}')"
  call="$(echo "$scored" | awk '{print $2}')"

  printf "%s\t%s\t%s\t%s\n" "$A" "$B" "$log2ratio" "$call" > "$OUTDIR/tmp/${sample}.meryl.tsv"
}

run_samples_single() {
  local fn="$1"
  if [[ "$JOBS" -gt 1 && -n "$(command -v parallel || true)" ]]; then
    echo "[WARN] GNU parallel for function calls is not enabled in this safe version; running sequentially."
  fi
  for s in "${SAMPLE_NAMES[@]}"; do
    "$fn" "$s"
  done
}

if [[ "$MODE" == "mash" || "$MODE" == "both" ]]; then
  do_mash_refs
fi
if [[ "$MODE" == "meryl" || "$MODE" == "both" ]]; then
  do_meryl_refs
fi

echo "[3/5] Scoring samples (mode=$MODE)"
if [[ "$MODE" == "mash" ]]; then
  run_samples_single mash_one_sample
elif [[ "$MODE" == "meryl" ]]; then
  run_samples_single meryl_one_sample
else
  for s in "${SAMPLE_NAMES[@]}"; do
    mash_one_sample "$s"
    meryl_one_sample "$s"
  done
fi

echo "[4/5] Writing results table"
RESULTS="$OUTDIR/results.tsv"
{
  echo -e "sample\tfastq_files\tmash_dist_to_Species1\tmash_dist_to_Species2\tmash_call\tmeryl_Aonly_hits\tmeryl_Bonly_hits\tmeryl_log2ratio\tmeryl_call\tfinal_call"
  for s in "${SAMPLE_NAMES[@]}"; do
    files="$(printf "%s" "${SAMPLE_FILES[$s]}" | paste -sd',' -)"
    mashA="NA"; mashB="NA"; mashC="NA"
    mA="NA"; mB="NA"; mL="NA"; mC="NA"

    if [[ -f "$OUTDIR/tmp/${s}.mash.tsv" ]]; then
      mashA="$(awk '{print $1}' "$OUTDIR/tmp/${s}.mash.tsv")"
      mashB="$(awk '{print $2}' "$OUTDIR/tmp/${s}.mash.tsv")"
      mashC="$(awk '{print $3}' "$OUTDIR/tmp/${s}.mash.tsv")"
    fi
    if [[ -f "$OUTDIR/tmp/${s}.meryl.tsv" ]]; then
      mA="$(awk '{print $1}' "$OUTDIR/tmp/${s}.meryl.tsv")"
      mB="$(awk '{print $2}' "$OUTDIR/tmp/${s}.meryl.tsv")"
      mL="$(awk '{print $3}' "$OUTDIR/tmp/${s}.meryl.tsv")"
      mC="$(awk '{print $4}' "$OUTDIR/tmp/${s}.meryl.tsv")"
    fi

    final="Ambiguous"
    if [[ "$mC" != "NA" && "$mashC" != "NA" && "$mC" == "$mashC" && "$mC" != "Ambiguous" ]]; then
      final="$mC"
    elif [[ "$mC" != "NA" && "$mC" != "Ambiguous" ]]; then
      final="$mC"
    elif [[ "$mashC" != "NA" && "$mashC" != "Ambiguous" ]]; then
      final="$mashC"
    fi

    echo -e "${s}\t${files}\t${mashA}\t${mashB}\t${mashC}\t${mA}\t${mB}\t${mL}\t${mC}\t${final}"
  done
} > "$RESULTS"

echo "[5/5] Done."
echo "Results: $RESULTS"
