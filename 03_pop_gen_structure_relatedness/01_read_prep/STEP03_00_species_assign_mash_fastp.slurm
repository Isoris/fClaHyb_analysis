###
# species_assign_mash_fastp_v2.slurm does this:
# Builds Mash sketches for Gar and Mac references, sketches all fastp read pairs,
# assigns each run-lane unit by Mash distance in parallel, and writes both
# per-run and collapsed per-CGA species assignment tables.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -c 127
#SBATCH --mem=237GB
#SBATCH -t 1-00:00:00
#SBATCH -J mash_assign
#SBATCH -o mash_assign.%j.out
#SBATCH -e mash_assign.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

# -----------------------
# Paths
# -----------------------
BASE="/scratch/lt200308-agbsci/Quentin_project"
FASTP_DIR="${BASE}/00-samples/fastp"
REF_GAR="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
REF_MAC="${BASE}/00-samples/fClaHyb_Mac_LG.fa"

OUTDIR="${BASE}/00-samples/species_assign_mash_fastp"
LOGDIR="${OUTDIR}/logs"
TMPDIR="${OUTDIR}/tmp"
SAMPDIR="${OUTDIR}/mash_samples"
REFDIR="${OUTDIR}/mash_refs"

mkdir -p "$OUTDIR" "$LOGDIR" "$TMPDIR" "$SAMPDIR" "$REFDIR"

timestamp() { date '+%F %T'; }

# -----------------------
# Mash parameters
# -----------------------
MASH="$(command -v mash)"
K=31
SKETCH=50000
MINCOPIES=2
MARGIN=0.002

# Parallelism tuning (IMPORTANT)
TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-127}"

# Use 8 threads per sample by default (safe). Change to 4 if filesystem is slow.
THREADS_PER_SAMPLE="${THREADS_PER_SAMPLE:-8}"

# Number of concurrent samples
NJOBS=$(( TOTAL_CPUS / THREADS_PER_SAMPLE ))
if [[ "$NJOBS" -lt 1 ]]; then NJOBS=1; fi

echo "[$(timestamp)] [INFO] Host=$(hostname)  SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"
echo "[$(timestamp)] [INFO] MASH=$MASH"
echo "[$(timestamp)] [INFO] FASTP_DIR=$FASTP_DIR"
echo "[$(timestamp)] [INFO] REF_GAR=$REF_GAR"
echo "[$(timestamp)] [INFO] REF_MAC=$REF_MAC"
echo "[$(timestamp)] [INFO] OUTDIR=$OUTDIR"
echo "[$(timestamp)] [INFO] K=$K SKETCH=$SKETCH MINCOPIES=$MINCOPIES MARGIN=$MARGIN"
echo "[$(timestamp)] [INFO] TOTAL_CPUS=$TOTAL_CPUS THREADS_PER_SAMPLE=$THREADS_PER_SAMPLE NJOBS=$NJOBS"

# -----------------------
# Build reference sketches (once)
# -----------------------
GAR_MSH="${REFDIR}/gar_ref.msh"
MAC_MSH="${REFDIR}/mac_ref.msh"

if [[ ! -s "$GAR_MSH" ]]; then
  echo "[$(timestamp)] [INFO] Building Gar reference sketch..."
  "$MASH" sketch -k "$K" -s "$SKETCH" -p "$THREADS_PER_SAMPLE" -o "${REFDIR}/gar_ref" "$REF_GAR" \
    > "${LOGDIR}/mash_ref_gar.log" 2>&1
fi

if [[ ! -s "$MAC_MSH" ]]; then
  echo "[$(timestamp)] [INFO] Building Mac reference sketch..."
  "$MASH" sketch -k "$K" -s "$SKETCH" -p "$THREADS_PER_SAMPLE" -o "${REFDIR}/mac_ref" "$REF_MAC" \
    > "${LOGDIR}/mash_ref_mac.log" 2>&1
fi

echo "[$(timestamp)] [INFO] Reference sketches ready."

# -----------------------
# Discover all fastp pairs
# -----------------------
IDS="${TMPDIR}/ids.txt"
: > "$IDS"

echo "[$(timestamp)] [INFO] Scanning FASTP_DIR for *.R1.fastp.fq.gz ..."
while IFS= read -r r1; do
  bn="$(basename "$r1")"
  id="${bn%.R1.fastp.fq.gz}"
  r2="${FASTP_DIR}/${id}.R2.fastp.fq.gz"
  if [[ -s "$r2" ]]; then
    echo "$id" >> "$IDS"
  else
    echo "[$(timestamp)] [WARN] Missing R2 for $id" >&2
  fi
done < <(find "$FASTP_DIR" -type f -name "*.R1.fastp.fq.gz" | sort)

sort -u "$IDS" -o "$IDS"
N=$(wc -l < "$IDS")
echo "[$(timestamp)] [INFO] Found $N run-lane units with pairs."

if [[ "$N" -eq 0 ]]; then
  echo "[$(timestamp)] [ERROR] No FASTP pairs found. Check FASTP_DIR=$FASTP_DIR" >&2
  exit 1
fi

# -----------------------
# Output tables
# -----------------------
RESULTS="${OUTDIR}/results.tsv"
echo -e "id\tR1\tR2\tdist_gar\tdist_mac\tcall\tabs_diff" > "$RESULTS"

# -----------------------
# Per-sample worker (runs in parallel)
# -----------------------
run_one() {
  id="$1"
  r1="${FASTP_DIR}/${id}.R1.fastp.fq.gz"
  r2="${FASTP_DIR}/${id}.R2.fastp.fq.gz"
  msh="${SAMPDIR}/${id}.msh"
  slog="${LOGDIR}/mash_${id}.log"

  {
    echo "[$(timestamp)] [INFO] ID=$id"
    echo "[$(timestamp)] [PATH] R1=$r1"
    echo "[$(timestamp)] [PATH] R2=$r2"
    echo "[$(timestamp)] [PATH] MSH=$msh"
    echo "[$(timestamp)] [INFO] THREADS_PER_SAMPLE=${THREADS_PER_SAMPLE}"
  } > "$slog"

  # sketch reads (skip if exists)
  if [[ ! -s "$msh" ]]; then
    "$MASH" sketch -r -k "$K" -s "$SKETCH" -m "$MINCOPIES" -p "$THREADS_PER_SAMPLE" \
      -o "${SAMPDIR}/${id}" "$r1" "$r2" >> "$slog" 2>&1
  else
    echo "[$(timestamp)] [SKIP] Sketch exists." >> "$slog"
  fi

  # distances
  dist_gar="$("$MASH" dist "$GAR_MSH" "$msh" | awk 'NR==1{print $3}')"
  dist_mac="$("$MASH" dist "$MAC_MSH" "$msh" | awk 'NR==1{print $3}')"

  # call with margin
  call="$(awk -v a="$dist_gar" -v b="$dist_mac" -v m="$MARGIN" 'BEGIN{
    if (a + m < b) print "Gariepinus";
    else if (b + m < a) print "Macrocephalus";
    else print "Ambiguous";
  }')"

  absdiff="$(awk -v a="$dist_gar" -v b="$dist_mac" 'BEGIN{
    d=a-b; if(d<0) d=-d; printf "%.6f", d
  }')"

  echo "[$(timestamp)] [DIST] gar=$dist_gar  mac=$dist_mac  call=$call  |diff|=$absdiff" >> "$slog"

  # IMPORTANT: print one line to stdout; parent appends to RESULTS
  echo -e "${id}\t${r1}\t${r2}\t${dist_gar}\t${dist_mac}\t${call}\t${absdiff}"
}

export -f run_one
export FASTP_DIR SAMPDIR LOGDIR MASH K SKETCH MINCOPIES GAR_MSH MAC_MSH MARGIN THREADS_PER_SAMPLE
export -f timestamp

# -----------------------
# Run all ids in parallel and append to RESULTS
# -----------------------
echo "[$(timestamp)] [INFO] Running mash per id in parallel..."
xargs -a "$IDS" -n 1 -P "$NJOBS" bash -lc 'run_one "$0"' >> "$RESULTS"

echo "[$(timestamp)] [INFO] Wrote per-run table: $RESULTS"

# -----------------------
# Collapse to per-CGA (median dist + majority call)
# -----------------------
PER_CGA="${OUTDIR}/results.per_CGA.tsv"
python3 - <<PY
import csv, statistics
from collections import defaultdict, Counter

infile = "${RESULTS}"
outfile = "${PER_CGA}"

rows = []
with open(infile, newline="") as f:
    r = csv.DictReader(f, delimiter="\\t")
    for row in r:
        rows.append(row)

by_cga = defaultdict(list)
for row in rows:
    cga = row["id"].split(".")[0]
    by_cga[cga].append(row)

def median(xs):
    xs = [float(x) for x in xs if x not in ("NA","")]
    return statistics.median(xs) if xs else float("nan")

with open(outfile, "w", newline="") as out:
    w = csv.writer(out, delimiter="\\t")
    w.writerow(["CGA_sample","n_runs","median_dist_gar","median_dist_mac","majority_call","calls_breakdown"])
    for cga, lst in sorted(by_cga.items()):
        n = len(lst)
        med_g = median([x["dist_gar"] for x in lst])
        med_m = median([x["dist_mac"] for x in lst])
        calls = [x["call"] for x in lst]
        c = Counter(calls)
        nonamb = [x for x in calls if x != "Ambiguous"]
        maj = Counter(nonamb).most_common(1)[0][0] if nonamb else "Ambiguous"
        breakdown = ",".join([f"{k}:{v}" for k,v in c.most_common()])
        w.writerow([cga,n,f"{med_g:.6f}",f"{med_m:.6f}",maj,breakdown])
PY

echo "[$(timestamp)] [INFO] Wrote per-CGA table: $PER_CGA"

echo "[$(timestamp)] [DONE] Quick summary:"
echo "  Per-run calls:"
cut -f6 "$RESULTS" | tail -n +2 | sort | uniq -c || true
echo "  Per-CGA calls:"
cut -f5 "$PER_CGA" | tail -n +2 | sort | uniq -c || true
