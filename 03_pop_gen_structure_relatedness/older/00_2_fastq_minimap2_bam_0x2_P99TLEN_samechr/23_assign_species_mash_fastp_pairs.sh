###
# assign_species_mash_fastp_pairs.sh does this:
# Builds Mash reference sketches for Gar and Mac genomes, sketches each fastp R1/R2 pair,
# and assigns each run-lane unit to Gariepinus, Macrocephalus, or Ambiguous.
###
#!/usr/bin/env bash
set -euo pipefail

K=31
SKETCH=50000
MINCOPIES=2
THREADS=16
MARGIN=0.002

usage() {
  cat <<EOF
Usage:
  $0 -a fClaHyb_Gar_LG.fa -b fClaHyb_Mac_LG.fa -r /path/to/fastp_dir -o outdir [options]

Options:
  -k INT      k-mer size (default: $K)
  -s INT      sketch size (default: $SKETCH)
  -m INT      mash -m min copies (default: $MINCOPIES)
  -t INT      threads (default: $THREADS)
  --margin X  decision margin (default: $MARGIN)

Output:
  outdir/results.tsv
EOF
}

REF_A=""
REF_B=""
READS_DIR=""
OUTDIR="mash_species_assign"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a) REF_A="$2"; shift 2;;
    -b) REF_B="$2"; shift 2;;
    -r) READS_DIR="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    -k) K="$2"; shift 2;;
    -s) SKETCH="$2"; shift 2;;
    -m) MINCOPIES="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    --margin) MARGIN="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

if [[ -z "$REF_A" || -z "$REF_B" || -z "$READS_DIR" ]]; then
  echo "ERROR: -a, -b, -r are required." >&2
  usage
  exit 1
fi

command -v mash >/dev/null 2>&1 || { echo "ERROR: mash not found in PATH"; exit 1; }

mkdir -p "$OUTDIR"/{01_refs,02_samples,logs,tmp}

timestamp() { date '+%F %T'; }

echo "[$(timestamp)] [INFO] REF_A=$REF_A"
echo "[$(timestamp)] [INFO] REF_B=$REF_B"
echo "[$(timestamp)] [INFO] READS_DIR=$READS_DIR"
echo "[$(timestamp)] [INFO] K=$K SKETCH=$SKETCH MINCOPIES=$MINCOPIES THREADS=$THREADS MARGIN=$MARGIN"
echo "[$(timestamp)] [INFO] OUTDIR=$OUTDIR"

A_MSH="$OUTDIR/01_refs/gar_ref.msh"
B_MSH="$OUTDIR/01_refs/mac_ref.msh"

if [[ ! -s "$A_MSH" ]]; then
  echo "[$(timestamp)] [INFO] Building reference sketch A (Gar)..."
  mash sketch -k "$K" -s "$SKETCH" -p "$THREADS" -o "$OUTDIR/01_refs/gar_ref" "$REF_A" \
    >"$OUTDIR/logs/mash_refA.log" 2>&1
fi
if [[ ! -s "$B_MSH" ]]; then
  echo "[$(timestamp)] [INFO] Building reference sketch B (Mac)..."
  mash sketch -k "$K" -s "$SKETCH" -p "$THREADS" -o "$OUTDIR/01_refs/mac_ref" "$REF_B" \
    >"$OUTDIR/logs/mash_refB.log" 2>&1
fi

echo "[$(timestamp)] [INFO] Scanning fastp dir for pairs..."
mapfile -t R1S < <(find "$READS_DIR" -type f -name "*.R1.fastp.fq.gz" | sort)
if [[ ${#R1S[@]} -eq 0 ]]; then
  echo "ERROR: no *.R1.fastp.fq.gz found in $READS_DIR" >&2
  exit 1
fi

IDS_FILE="$OUTDIR/tmp/ids.txt"
: > "$IDS_FILE"
for r1 in "${R1S[@]}"; do
  bn="$(basename "$r1")"
  id="${bn%.R1.fastp.fq.gz}"
  r2="$READS_DIR/${id}.R2.fastp.fq.gz"
  if [[ -s "$r2" ]]; then
    echo "$id" >> "$IDS_FILE"
  else
    echo "[$(timestamp)] [WARN] Missing R2 for ID=$id (expected: $r2)" >&2
  fi
done

sort -u "$IDS_FILE" > "$OUTDIR/tmp/ids.unique.txt"
N=$(wc -l < "$OUTDIR/tmp/ids.unique.txt")
echo "[$(timestamp)] [INFO] Found $N paired run-lane units"

RESULTS="$OUTDIR/results.tsv"
echo -e "id\tR1\tR2\tdist_gar\tdist_mac\tcall\twinner_margin" > "$RESULTS"

while read -r id; do
  r1="$READS_DIR/${id}.R1.fastp.fq.gz"
  r2="$READS_DIR/${id}.R2.fastp.fq.gz"
  msh="$OUTDIR/02_samples/${id}.msh"
  slog="$OUTDIR/logs/mash_${id}.log"

  echo "[$(timestamp)] [RUN] $id"
  echo "[$(timestamp)] [PATH] R1=$r1" > "$slog"
  echo "[$(timestamp)] [PATH] R2=$r2" >> "$slog"

  if [[ ! -s "$msh" ]]; then
    mash sketch -r -k "$K" -s "$SKETCH" -m "$MINCOPIES" -p "$THREADS" \
      -o "$OUTDIR/02_samples/${id}" "$r1" "$r2" >> "$slog" 2>&1
  fi

  distA="$(mash dist "$A_MSH" "$msh" | awk 'NR==1{print $3}')"
  distB="$(mash dist "$B_MSH" "$msh" | awk 'NR==1{print $3}')"

  call="$(awk -v a="$distA" -v b="$distB" -v m="$MARGIN" 'BEGIN{
    if (a + m < b) print "Gariepinus";
    else if (b + m < a) print "Macrocephalus";
    else print "Ambiguous";
  }')"

  winner_margin="$(awk -v a="$distA" -v b="$distB" 'BEGIN{
    d=a-b; if(d<0) d=-d; printf "%.6f", d
  }')"

  echo -e "${id}\t${r1}\t${r2}\t${distA}\t${distB}\t${call}\t${winner_margin}" >> "$RESULTS"
done < "$OUTDIR/tmp/ids.unique.txt"

echo "[$(timestamp)] [DONE] Wrote: $RESULTS"
echo "Quick check:"
echo "  cut -f6 $RESULTS | tail -n +2 | sort | uniq -c"
