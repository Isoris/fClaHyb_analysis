###
# dedup_all_batches_array.slurm does this:
# Runs a SLURM array over BAMs from batches 00-02, performs name-sort,
# fixmate, coordinate sort, and samtools markdup, then writes deduplicated
# or duplicate-marked BAMs with flagstat and completion checkpoints.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=120GB
#SBATCH -t 2-00:00:00
#SBATCH -J dedup_all
#SBATCH -o dedup_all.%A_%a.out
#SBATCH -e dedup_all.%A_%a.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project/01-minimap2-bams"
QCBASE="${BASE}/qc_dedup"
mkdir -p "$QCBASE"

# If you want to REMOVE duplicates, set REMOVE_DUPLICATES=1 when submitting
REMOVE_DUPLICATES="${REMOVE_DUPLICATES:-0}"

THREADS="${SLURM_NTASKS:-32}"

# Memory per thread for samtools sort
SORT_MEM_PER_THREAD="2500M"   # 2.5G * 32 = 80G (safe under 120G)

timestamp(){ date '+%F %T'; }

# -----------------------
# Build file list once (shared)
# -----------------------
LIST="${QCBASE}/bam_list.batches00-02.txt"
if [[ ! -s "$LIST" ]]; then
  echo "[$(timestamp)] [INFO] Building BAM list at $LIST"
  find "$BASE"/batch0{0,1,2} -type f -name "*.sr.bam" | sort > "$LIST"
fi

N=$(wc -l < "$LIST")
if [[ "${SLURM_ARRAY_TASK_ID:-}" == "" ]]; then
  echo "[$(timestamp)] [ERROR] This script is meant to run as an array job."
  echo "Hint: sbatch --array=1-$N dedup_all_batches_array.slurm"
  exit 1
fi

BAM_IN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
if [[ -z "$BAM_IN" ]]; then
  echo "[$(timestamp)] [ERROR] No BAM for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

BN="$(basename "$BAM_IN" .bam)"
BATCH="$(basename "$(dirname "$BAM_IN")")"
OUTDIR="${QCBASE}/${BATCH}/${BN}"
mkdir -p "$OUTDIR"

# Final outputs
if [[ "$REMOVE_DUPLICATES" -eq 1 ]]; then
  FINAL_BAM="${OUTDIR}/${BN}.rmdup.bam"
else
  FINAL_BAM="${OUTDIR}/${BN}.markdup.bam"
fi
FINAL_BAI="${FINAL_BAM}.bai"
DONE="${OUTDIR}/${BN}.DONE.txt"
FLAGSTAT="${OUTDIR}/${BN}.flagstat.txt"

# Skip if already done
if [[ -s "$FINAL_BAM" && -s "$FINAL_BAI" && -s "$DONE" ]]; then
  echo "[$(timestamp)] [SKIP] Already finished: $FINAL_BAM"
  exit 0
fi

# tmpdir (avoid /tmp permission issues)
TMPDIR="${OUTDIR}/tmp_${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

echo "[$(timestamp)] [INFO] Host=$(hostname) JOB=${SLURM_JOB_ID:-NA} TASK=${SLURM_ARRAY_TASK_ID}"
echo "[$(timestamp)] [INFO] BAM_IN=$BAM_IN"
echo "[$(timestamp)] [INFO] OUTDIR=$OUTDIR"
echo "[$(timestamp)] [INFO] THREADS=$THREADS MEM/thread=$SORT_MEM_PER_THREAD"
echo "[$(timestamp)] [INFO] REMOVE_DUPLICATES=$REMOVE_DUPLICATES"
echo "[$(timestamp)] [INFO] TMPDIR=$TMPDIR"

# Intermediates (will be deleted)
NSORT_BAM="${TMPDIR}/${BN}.namesort.bam"
FIXM_BAM="${TMPDIR}/${BN}.fixmate.bam"
CSORT_BAM="${TMPDIR}/${BN}.possort.bam"

# 1) name-sort
echo "[$(timestamp)] [STEP] name-sort"
samtools sort -n -@ "$THREADS" -m "$SORT_MEM_PER_THREAD" -T "${TMPDIR}/nsort" \
  -o "$NSORT_BAM" "$BAM_IN"

# 2) fixmate
echo "[$(timestamp)] [STEP] fixmate"
samtools fixmate -m -@ "$THREADS" "$NSORT_BAM" "$FIXM_BAM"

# 3) coord-sort
echo "[$(timestamp)] [STEP] coord-sort"
samtools sort -@ "$THREADS" -m "$SORT_MEM_PER_THREAD" -T "${TMPDIR}/csort" \
  -o "$CSORT_BAM" "$FIXM_BAM"

# 4) markdup (or remove)
if [[ "$REMOVE_DUPLICATES" -eq 1 ]]; then
  echo "[$(timestamp)] [STEP] markdup + remove"
  samtools markdup -@ "$THREADS" -r -s "$CSORT_BAM" "$FINAL_BAM"
else
  echo "[$(timestamp)] [STEP] markdup (mark only)"
  samtools markdup -@ "$THREADS" -s "$CSORT_BAM" "$FINAL_BAM"
fi

# 5) index + flagstat
echo "[$(timestamp)] [STEP] index"
samtools index -@ "$THREADS" "$FINAL_BAM"

echo "[$(timestamp)] [STEP] flagstat"
samtools flagstat -@ "$THREADS" "$FINAL_BAM" > "$FLAGSTAT"

# checkpoint
{
  echo "DONE $(timestamp)"
  echo "BAM_IN=$BAM_IN"
  echo "FINAL_BAM=$FINAL_BAM"
  echo "THREADS=$THREADS"
  echo "REMOVE_DUPLICATES=$REMOVE_DUPLICATES"
} > "$DONE"

# cleanup intermediates
echo "[$(timestamp)] [STEP] cleanup tmp"
rm -rf "$TMPDIR"

echo "[$(timestamp)] [DONE] $FINAL_BAM"
