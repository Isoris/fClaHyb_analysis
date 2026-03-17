###
# merge_markdup_clip_qc_array.slurm does this:
# Merges per-run BAMs into one BAM per sample, validates RG/SM consistency,
# performs fixmate + markdup, optionally runs clipOverlap, generates QC outputs,
# and appends a summary row for each sample in a SLURM array workflow.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=120GB
#SBATCH -t 2-00:00:00
#SBATCH -J merge_md_clip
#SBATCH -o merge_md_clip.%A_%a.out
#SBATCH -e merge_md_clip.%A_%a.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

############################
# CONFIG
############################
BAMBASE="/scratch/lt200308-agbsci/Quentin_project/01-minimap2-bams"
OUTBASE="/scratch/lt200308-agbsci/Quentin_project/02-merged_per_sample"

# Optional (only used if you later add extra metrics)
REF="/scratch/lt200308-agbsci/Quentin_project/00-samples/fClaHyb_Gar_LG.fa"

# Threads
THREADS="${SLURM_NTASKS:-16}"

# Toggle clipOverlap (1=on, 0=off)
DO_CLIP="${DO_CLIP:-1}"

# Tools
SAMTOOLS="$(command -v samtools)"
BAMUTIL="$(command -v bam)"   # BamUtil clipOverlap is "bam"

timestamp(){ date '+%F %T'; }

mkdir -p "$OUTBASE"/{lists,logs,tmp,results}

############################
# SAFE SORT MEMORY
############################
# samtools sort memory = THREADS * (-m) + overhead.
# Compute per-thread memory from SLURM --mem and clamp it.
mem_mb="${SLURM_MEM_PER_NODE:-0}"
if [[ "$mem_mb" -eq 0 ]]; then
  # fallback: parse from /proc if SLURM env missing
  mem_mb=$(awk '/MemTotal/ {print int($2/1024)}' /proc/meminfo)
fi

# Reserve 20% for overhead + other steps, then divide by threads.
# Clamp between 256M and 1500M per thread.
usable_mb=$(( mem_mb * 80 / 100 ))
m_per_thread=$(( usable_mb / THREADS ))
if (( m_per_thread < 256 )); then m_per_thread=256; fi
if (( m_per_thread > 1500 )); then m_per_thread=1500; fi
SORT_MEM_PER_THREAD="${m_per_thread}M"

############################
# BUILD INPUT LISTS (once)
############################
ALLBAMS="${OUTBASE}/lists/all_run_bams.txt"
if [[ ! -s "$ALLBAMS" ]]; then
  find "$BAMBASE"/batch0{0,1,2,3,4,5,6,7} -type f -name "*.sr.bam" ! -name "*.tmp.*" \
    | sort > "$ALLBAMS"
fi

SAMPLES="${OUTBASE}/lists/samples.txt"
if [[ ! -s "$SAMPLES" ]]; then
  awk -F'/' '{print $NF}' "$ALLBAMS" \
    | sed 's/\.sr\.bam$//' \
    | awk -F'.' '{print $1}' \
    | sort -u > "$SAMPLES"
fi

NSAMP=$(wc -l < "$SAMPLES")
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "[$(timestamp)] [ERROR] Run as array: sbatch --array=1-${NSAMP} $0"
  exit 1
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES")
[[ -n "$SAMPLE" ]] || { echo "[$(timestamp)] [ERROR] Empty SAMPLE for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

############################
# PER-SAMPLE FILE SET
############################
SDIR="${OUTBASE}/${SAMPLE}"
mkdir -p "$SDIR"

RUNLIST="${SDIR}/${SAMPLE}.run_bams.txt"
grep -F "/${SAMPLE}." "$ALLBAMS" > "$RUNLIST" || true
if [[ ! -s "$RUNLIST" ]]; then
  echo "[$(timestamp)] [ERROR] No run BAMs found for SAMPLE=$SAMPLE"
  exit 2
fi

TMPDIR="${OUTBASE}/tmp/${SAMPLE}.tmp_${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
export TMPDIR

LOG="${SDIR}/${SAMPLE}.pipeline.log"

FINAL_BAM="${SDIR}/${SAMPLE}.merged.markdup${DO_CLIP:+.clip}.bam"
FINAL_BAI="${FINAL_BAM}.bai"

QC_ROW="${OUTBASE}/results/QC_ROWS.tsv"

# Skip if done
if [[ -s "$FINAL_BAM" && -s "$FINAL_BAI" ]]; then
  echo "[$(timestamp)] [SKIP] $SAMPLE already finished" | tee -a "$LOG"
  exit 0
fi

echo "[$(timestamp)] [INFO] SAMPLE=$SAMPLE task=${SLURM_ARRAY_TASK_ID}/${NSAMP}" | tee -a "$LOG"
echo "[$(timestamp)] [INFO] THREADS=$THREADS mem_mb=$mem_mb sort_mem_per_thread=$SORT_MEM_PER_THREAD" | tee -a "$LOG"
echo "[$(timestamp)] [INFO] Using $(wc -l < "$RUNLIST") run BAMs" | tee -a "$LOG"

############################
# 0) Validate RG: SM == SAMPLE
############################
echo "[$(timestamp)] [STEP] validate @RG SM tags" | tee -a "$LOG"
BADSM=0
while read -r bam; do
  sm_list=$("$SAMTOOLS" view -H "$bam" | awk -F'\t' '
    $1=="@RG"{
      for(i=1;i<=NF;i++){
        if($i ~ /^SM:/){sub(/^SM:/,"",$i); print $i}
      }
    }' | sort -u)

  if [[ -z "$sm_list" ]]; then
    echo "[$(timestamp)] [ERROR] No @RG SM found in $bam" | tee -a "$LOG"
    BADSM=1
    continue
  fi

  while read -r sm; do
    if [[ "$sm" != "$SAMPLE" ]]; then
      echo "[$(timestamp)] [ERROR] SM mismatch in $bam : SM=$sm expected=$SAMPLE" | tee -a "$LOG"
      BADSM=1
    fi
  done <<< "$sm_list"
done < "$RUNLIST"

if [[ "$BADSM" -ne 0 ]]; then
  echo "[$(timestamp)] [FAIL] RG/SM mismatch. Fix RG first then rerun." | tee -a "$LOG"
  exit 3
fi

############################
# 1) Merge run BAMs (keep RG)
############################
MERGED="${TMPDIR}/${SAMPLE}.merged.bam"
echo "[$(timestamp)] [STEP] samtools merge" | tee -a "$LOG"
"$SAMTOOLS" merge -@ "$THREADS" -f -b "$RUNLIST" -o "$MERGED"

############################
# 2) name-sort -> fixmate -> coord-sort
############################
NSORT="${TMPDIR}/${SAMPLE}.namesort.bam"
FIXM="${TMPDIR}/${SAMPLE}.fixmate.bam"
CSORT="${TMPDIR}/${SAMPLE}.possort.bam"

echo "[$(timestamp)] [STEP] name-sort" | tee -a "$LOG"
"$SAMTOOLS" sort -n -@ "$THREADS" -m "$SORT_MEM_PER_THREAD" -T "${TMPDIR}/nsort" -o "$NSORT" "$MERGED"

echo "[$(timestamp)] [STEP] fixmate" | tee -a "$LOG"
"$SAMTOOLS" fixmate -m -@ "$THREADS" "$NSORT" "$FIXM"

echo "[$(timestamp)] [STEP] coord-sort" | tee -a "$LOG"
"$SAMTOOLS" sort -@ "$THREADS" -m "$SORT_MEM_PER_THREAD" -T "${TMPDIR}/csort" -o "$CSORT" "$FIXM"

############################
# 3) markdup (mark only)
############################
DUP="${TMPDIR}/${SAMPLE}.markdup.bam"
echo "[$(timestamp)] [STEP] markdup (mark only)" | tee -a "$LOG"
"$SAMTOOLS" markdup -@ "$THREADS" -s "$CSORT" "$DUP" 2> "${SDIR}/${SAMPLE}.markdup.stderr.txt"

############################
# 4) clipOverlap (optional)
############################
if [[ "$DO_CLIP" -eq 1 ]]; then
  echo "[$(timestamp)] [STEP] clipOverlap" | tee -a "$LOG"
  "$BAMUTIL" clipOverlap --in "$DUP" --out "$FINAL_BAM" --storeOrig OC --stats --unmapped \
    > "${SDIR}/${SAMPLE}.clipOverlap.stdout.txt" 2> "${SDIR}/${SAMPLE}.clipOverlap.stats.txt"
else
  echo "[$(timestamp)] [STEP] skip clipOverlap (DO_CLIP=0)" | tee -a "$LOG"
  mv "$DUP" "$FINAL_BAM"
fi

############################
# 5) QC: index + quickcheck + flagstat + stats
############################
echo "[$(timestamp)] [STEP] index" | tee -a "$LOG"
"$SAMTOOLS" index -@ "$THREADS" "$FINAL_BAM"

echo "[$(timestamp)] [STEP] quickcheck" | tee -a "$LOG"
QC_STATUS="OK"
if ! "$SAMTOOLS" quickcheck -v "$FINAL_BAM" > "${SDIR}/${SAMPLE}.quickcheck.txt" 2>&1; then
  QC_STATUS="BAD"
fi

echo "[$(timestamp)] [STEP] flagstat" | tee -a "$LOG"
"$SAMTOOLS" flagstat -@ "$THREADS" "$FINAL_BAM" > "${SDIR}/${SAMPLE}.flagstat.txt"

echo "[$(timestamp)] [STEP] samtools stats" | tee -a "$LOG"
"$SAMTOOLS" stats -@ "$THREADS" "$FINAL_BAM" > "${SDIR}/${SAMPLE}.samtools.stats.txt"

# Parse quick metrics
total_reads=$(awk '/in total/ {print $1; exit}' "${SDIR}/${SAMPLE}.flagstat.txt" || echo "NA")
mapped_reads=$(awk '/ mapped \(/ && $0 !~ /primary/ {print $1; exit}' "${SDIR}/${SAMPLE}.flagstat.txt" || echo "NA")
pct_mapped=$(awk '/ mapped \(/ && $0 !~ /primary/ {gsub(/[()%]/,""); print $5; exit}' "${SDIR}/${SAMPLE}.flagstat.txt" || echo "NA")
ins_mean=$(awk -F'\t' '$1=="SN" && $2~/(insert size average)/ {print $3; exit}' "${SDIR}/${SAMPLE}.samtools.stats.txt" || echo "NA")
err_rate=$(awk -F'\t' '$1=="SN" && $2~/(error rate)/ {print $3; exit}' "${SDIR}/${SAMPLE}.samtools.stats.txt" || echo "NA")

# Append row (safe)
mkdir -p "$(dirname "$QC_ROW")"
{
  flock -x 9
  if [[ ! -s "$QC_ROW" ]]; then
    printf "bam\tsample\tn_run_bams\tstatus\tthreads\tsort_mem_per_thread\tn_reads\tn_mapped\tpct_mapped\tinsert_mean\terror_rate\n" > "$QC_ROW"
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$FINAL_BAM" "$SAMPLE" "$(wc -l < "$RUNLIST")" "$QC_STATUS" \
    "$THREADS" "$SORT_MEM_PER_THREAD" \
    "$total_reads" "$mapped_reads" "$pct_mapped" "$ins_mean" "$err_rate" \
    >> "$QC_ROW"
} 9>>"${OUTBASE}/results/.lock"

echo "[$(timestamp)] [DONE] SAMPLE=$SAMPLE status=$QC_STATUS" | tee -a "$LOG"

# Cleanup
rm -rf "$TMPDIR"
