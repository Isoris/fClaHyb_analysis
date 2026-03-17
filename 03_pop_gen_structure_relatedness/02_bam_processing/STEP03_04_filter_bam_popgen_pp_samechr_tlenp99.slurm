###
# bam_filter_array_redo_0x2_TLEN.slurm does this:
# Runs a SLURM array over merged markdup+clip BAMs, filters reads by MAPQ,
# proper pairing, same-chromosome mate, and absolute TLEN bounds, then
# writes filtered BAMs plus index and flagstat outputs.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64GB
#SBATCH -t 1-00:00:00
#SBATCH -J bam_filter_pp_tlen
#SBATCH -o /scratch/lt200308-agbsci/Quentin_project/02-merged_per_sample/results/bam_filter_pp_tlen.%A_%a.out
#SBATCH -e /scratch/lt200308-agbsci/Quentin_project/02-merged_per_sample/results/bam_filter_pp_tlen.%A_%a.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

SAMTOOLS="$(command -v samtools)"

BASE="/scratch/lt200308-agbsci/Quentin_project"
IN_GLOB="${BASE}/02-merged_per_sample/*/*.merged.markdup.clip.bam"
RESULTS_DIR="${BASE}/02-merged_per_sample/results"
THREADS="${SLURM_NTASKS:-16}"

# ----------------------------
# Filtering params
# ----------------------------
MINMAPQ=60

# Strict popgen: drop unmapped, mate-unmapped, secondary, QCfail, duplicates, supplementary
# 0x4 + 0x8 + 0x100 + 0x200 + 0x400 + 0x800 = 0xF0C
REMOVE_MASK=0xF0C

# Enforce proper pairs
REQUIRE_PROPER_PAIRS=1

# Enforce same-chr mate (RNEXT "=" or equals RNAME)
ENFORCE_SAME_CHR=1

# TLEN (abs) band: use empirical p99
ENFORCE_TLEN=1
MIN_TLEN=0
MAX_TLEN=514   # p99(abs(TLEN)) = 514 bp (p95=446)

timestamp(){ date '+%F %T'; }
mkdir -p "$RESULTS_DIR"

mapfile -t BAMS < <(ls -1 $IN_GLOB 2>/dev/null | sort)
NBAM="${#BAMS[@]}"
if [[ "$NBAM" -eq 0 ]]; then
  echo "[$(timestamp)] [ERROR] No BAMs found with: $IN_GLOB" >&2
  exit 2
fi

BAM="${BAMS[$SLURM_ARRAY_TASK_ID]}"
SDIR="$(dirname "$BAM")"
BN="$(basename "$BAM")"
SAMPLE="${BN%%.*}"

OUT="${SDIR}/${SAMPLE}.merged.markdup.clip.pp.samechr.tlenP99.filtered.bam"
LOG="${SDIR}/${SAMPLE}.filter.pp.samechr.tlenP99.pipeline.log"

echo "[$(timestamp)] [INFO] SAMPLE=$SAMPLE" | tee "$LOG"
echo "[$(timestamp)] [INFO] BAM=$BAM" | tee -a "$LOG"
echo "[$(timestamp)] [INFO] OUT=$OUT" | tee -a "$LOG"
echo "[$(timestamp)] [INFO] THREADS=$THREADS MINMAPQ=$MINMAPQ REMOVE_MASK=$REMOVE_MASK REQUIRE_PROPER_PAIRS=$REQUIRE_PROPER_PAIRS ENFORCE_SAME_CHR=$ENFORCE_SAME_CHR ENFORCE_TLEN=$ENFORCE_TLEN MIN_TLEN=$MIN_TLEN MAX_TLEN=$MAX_TLEN" | tee -a "$LOG"

[[ -s "$BAM" ]] || { echo "[$(timestamp)] [ERROR] missing BAM=$BAM" | tee -a "$LOG"; exit 3; }

echo "[$(timestamp)] [STEP] samtools view + awk filters" | tee -a "$LOG"

VIEW_ARGS=(-@ "$THREADS" -h -F "$REMOVE_MASK" -q "$MINMAPQ")
if [[ "$REQUIRE_PROPER_PAIRS" -eq 1 ]]; then
  VIEW_ARGS+=(-f 0x2)
fi

"$SAMTOOLS" view "${VIEW_ARGS[@]}" "$BAM" \
| awk -v samechr="$ENFORCE_SAME_CHR" -v dotlen="$ENFORCE_TLEN" -v min="$MIN_TLEN" -v max="$MAX_TLEN" '
  BEGIN{OFS="\t"}
  /^@/ {print; next}
  {
    if (samechr==1) {
      if (!($7=="=" || $7==$3)) next
    }
    if (dotlen==1) {
      t=$9; if (t<0) t=-t
      if (t<min || t>max) next
    }
    print
  }' \
| "$SAMTOOLS" view -@ "$THREADS" -bh -o "$OUT" - 2>>"$LOG"

echo "[$(timestamp)] [STEP] index" | tee -a "$LOG"
"$SAMTOOLS" index -@ "$THREADS" "$OUT" >> "$LOG" 2>&1

echo "[$(timestamp)] [STEP] flagstat" | tee -a "$LOG"
"$SAMTOOLS" flagstat -@ "$THREADS" "$OUT" > "${SDIR}/${SAMPLE}.pp.samechr.tlenP99.filtered.flagstat.txt"

echo "[$(timestamp)] [DONE] OK" | tee -a "$LOG"
