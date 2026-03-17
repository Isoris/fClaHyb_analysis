###
# fastp_redo_fullnode.slurm does this:
# Re-runs fastp on a TSV-defined set of read pairs (typically redo cases such as
# unexpected EOF), using parallel per-sample jobs on a full node and writing
# fresh fastp outputs plus per-sample redo logs.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -c 127
#SBATCH --mem=237GB
#SBATCH -t 1-00:00:00
#SBATCH -J fastp_redo
#SBATCH -o fastp_redo.%j.out
#SBATCH -e fastp_redo.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

FASTP="$(command -v fastp)"

BASE="/scratch/lt200308-agbsci/Quentin_project"
INDIR="${BASE}/00-samples"
OUTFQDIR="${INDIR}/fastp"
LOGDIR="${OUTFQDIR}/logs"
mkdir -p "$OUTFQDIR" "$LOGDIR"

# REQUIRED: set TSV at submit time
# sbatch --export=ALL,TSV=/path/to/cga_fastq.batch.00.redo_unexpected_eof.tsv fastp_redo_fullnode.slurm
: "${TSV:?ERROR: set TSV in sbatch --export (TSV=/path/to/file.tsv)}"

timestamp(){ date '+%F %T'; }

# Total CPUs Slurm gave you
CPUS="${SLURM_CPUS_PER_TASK:-1}"

# Tune parallelism:
# - Each fastp uses TFASTP threads (-w)
# - We run NJOBS fastp processes concurrently
# Make sure NJOBS*TFASTP <= CPUS
TFASTP="${TFASTP:-8}"                 # threads per fastp
NJOBS="${NJOBS:-$(( CPUS / TFASTP ))}" # number of concurrent fastp

if [[ "$NJOBS" -lt 1 ]]; then NJOBS=1; fi

echo "[$(timestamp)] [INFO] host=$(hostname)"
echo "[$(timestamp)] [INFO] FASTP=$FASTP"
echo "[$(timestamp)] [INFO] TSV=$TSV"
echo "[$(timestamp)] [INFO] SLURM_CPUS_PER_TASK=$CPUS"
echo "[$(timestamp)] [INFO] NJOBS=$NJOBS TFASTP=$TFASTP total_threads=$((NJOBS*TFASTP))"
echo "[$(timestamp)] [INFO] OUTFQDIR=$OUTFQDIR LOGDIR=$LOGDIR"

fastp_one () {
  local SAMPLE="$1" RUN="$2" LANE="$3" R1="$4" R2="$5"
  local ID="${SAMPLE}.${RUN}.${LANE}"

  local OUT1="${OUTFQDIR}/${ID}.R1.fastp.fq.gz"
  local OUT2="${OUTFQDIR}/${ID}.R2.fastp.fq.gz"
  local LOG="${LOGDIR}/${ID}.fastp.redo.log"

  {
    echo "[$(timestamp)] [START] $ID"
    echo "[$(timestamp)] [IN]  R1=$R1"
    echo "[$(timestamp)] [IN]  R2=$R2"
    echo "[$(timestamp)] [OUT] O1=$OUT1"
    echo "[$(timestamp)] [OUT] O2=$OUT2"

    [[ -s "$R1" ]] || { echo "[$(timestamp)] [ERROR] missing R1"; exit 2; }
    [[ -s "$R2" ]] || { echo "[$(timestamp)] [ERROR] missing R2"; exit 2; }

    # For unexpected EOF cases: remove old corrupted outputs and redo
    rm -f "$OUT1" "$OUT2"

    fastp \
      -i "$R1" -I "$R2" \
      -o "$OUT1" -O "$OUT2" \
      --detect_adapter_for_pe \
      --trim_poly_g --trim_poly_x \
      -5 -3 -q 20 \
      -f 5 -F 5 -t 5 -T 5 \
      -n 0 -l 30 \
      -w "$TFASTP"

    echo "[$(timestamp)] [DONE] OK"
  } > "$LOG" 2>&1
}

export -f fastp_one timestamp
export TFASTP OUTFQDIR LOGDIR

NROWS=$(( $(wc -l < "$TSV") - 1 ))
echo "[$(timestamp)] [INFO] rows(excl header)=$NROWS"

i=0
tail -n +2 "$TSV" | while IFS=$'\t' read -r SAMPLE RUN LANE R1 R2; do
  i=$((i+1))
  while [[ $(jobs -rp | wc -l) -ge "$NJOBS" ]]; do
    sleep 1
  done
  echo "[$(timestamp)] [QUEUE] ($i/$NROWS) ${SAMPLE}.${RUN}.${LANE}"
  fastp_one "$SAMPLE" "$RUN" "$LANE" "$R1" "$R2" &
done

wait
echo "[$(timestamp)] [INFO] finished fastp redo for TSV=$TSV"
