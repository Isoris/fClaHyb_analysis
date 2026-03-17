###
# map_minimap2_batch00.slurm does this:
# Maps fastp-processed paired-end reads from batch00 with minimap2, adds read-group
# information, sorts and indexes each per-run BAM, and writes one log per run.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237GB
#SBATCH -t 5-00:00:00
##SBATCH -A lt200308
#SBATCH -J map_mm2_b00
#SBATCH -o map_mm2_b00.%j.out
#SBATCH -e map_mm2_b00.%j.err

set -euo pipefail

source ~/.bashrc
mamba activate assembly

MINIMAP2="/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/minimap2"
SAMTOOLS="$(command -v samtools)"

BASE="/scratch/lt200308-agbsci/Quentin_project"
INDIR="${BASE}/00-samples"
REF="${INDIR}/fClaHyb_Gar_LG.fa"
TSV="${INDIR}/cga_fastq.batch.00.tsv"
OUTDIR="${BASE}/01-minimap2-bams/batch00"
LOGDIR="${OUTDIR}/logs"

mkdir -p "$OUTDIR" "$LOGDIR"

NJOBS=8
TMAP=16

timestamp() { date '+%F %T'; }

echo "[$(timestamp)] [INFO] Job started on host=$(hostname)"
echo "[$(timestamp)] [INFO] MINIMAP2=$MINIMAP2"
echo "[$(timestamp)] [INFO] SAMTOOLS=$SAMTOOLS"
echo "[$(timestamp)] [INFO] REF=$REF"
echo "[$(timestamp)] [INFO] TSV=$TSV"
echo "[$(timestamp)] [INFO] OUTDIR=$OUTDIR"
echo "[$(timestamp)] [INFO] NJOBS=$NJOBS  TMAP=$TMAP  (approx threads used = $((NJOBS*TMAP)))"
echo "[$(timestamp)] [INFO] SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"

# Build index once
if [[ ! -s "${REF}.mmi" ]]; then
  echo "[$(timestamp)] [INFO] Building minimap2 index: ${REF}.mmi"
  "$MINIMAP2" -d "${REF}.mmi" "$REF"
  echo "[$(timestamp)] [INFO] Index built"
else
  echo "[$(timestamp)] [INFO] Index exists: ${REF}.mmi"
fi

map_one () {
  local SAMPLE="$1" RUN="$2" LANE="$3" R1="$4" R2="$5"

  # override raw FASTQs with fastp outputs
  local FASTP_DIR="${INDIR}/fastp"
  R1="${FASTP_DIR}/${SAMPLE}.${RUN}.${LANE}.R1.fastp.fq.gz"
  R2="${FASTP_DIR}/${SAMPLE}.${RUN}.${LANE}.R2.fastp.fq.gz"

  local SM="$SAMPLE"
  local ID="${SAMPLE}.${RUN}.${LANE}"
  local LB="$SAMPLE"
  local PU="${RUN}.${LANE}"

  local OUTBAM="${OUTDIR}/${ID}.sr.bam"
  local ROWLOG="${LOGDIR}/${ID}.log"

  {
    echo "[$(timestamp)] [START] ID=$ID  SM=$SM  RUN=$RUN  LANE=$LANE"
    echo "[$(timestamp)] [PATHS] R1=$R1"
    echo "[$(timestamp)] [PATHS] R2=$R2"
    echo "[$(timestamp)] [RG] @RG\\tID:${ID}\\tSM:${SM}\\tLB:${LB}\\tPU:${PU}\\tPL:ILLUMINA"
    echo "[$(timestamp)] [OUT]  OUTBAM=$OUTBAM"

    if [[ -s "$OUTBAM" && -s "${OUTBAM}.bai" ]]; then
      echo "[$(timestamp)] [SKIP] BAM+BAI already exist"
      exit 0
    fi

    echo "[$(timestamp)] [CMD] minimap2 -> samtools sort"
    "$MINIMAP2" -ax sr \
      -Y -c \
      --secondary-seq --sam-hit-only \
      -t "$TMAP" \
      -R "@RG\\tID:${ID}\\tSM:${SM}\\tLB:${LB}\\tPU:${PU}\\tPL:ILLUMINA" \
      "${REF}.mmi" "$R1" "$R2" \
    | "$SAMTOOLS" sort -@ "$TMAP" -o "$OUTBAM" -

    echo "[$(timestamp)] [CMD] samtools index"
    "$SAMTOOLS" index -@ "$TMAP" "$OUTBAM"

    echo "[$(timestamp)] [DONE] OK"
  } > "$ROWLOG" 2>&1
}

export -f map_one
export MINIMAP2 SAMTOOLS REF OUTDIR LOGDIR TMAP
export -f timestamp

# Count rows (minus header) for quick progress context
NROWS=$(( $(wc -l < "$TSV") - 1 ))
echo "[$(timestamp)] [INFO] TSV rows to process (excluding header) = $NROWS"

i=0
tail -n +2 "$TSV" \
| while IFS=$'\t' read -r SAMPLE RUN LANE R1 R2; do
    i=$((i+1))
    while [[ $(jobs -rp | wc -l) -ge "$NJOBS" ]]; do
      echo "[$(timestamp)] [THROTTLE] running=$(jobs -rp | wc -l) >= NJOBS=$NJOBS ; sleeping..."
      sleep 3
    done

    ID="${SAMPLE}.${RUN}.${LANE}"
    echo "[$(timestamp)] [QUEUE] ($i/$NROWS) launching $ID  threads=$TMAP  running_now=$(jobs -rp | wc -l)"
    map_one "$SAMPLE" "$RUN" "$LANE" "$R1" "$R2" &
  done

wait
echo "[$(timestamp)] [INFO] All mappings finished for batch00"
echo "[$(timestamp)] [INFO] Per-run logs in: $LOGDIR"
