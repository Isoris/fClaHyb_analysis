#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1 -n 127
#SBATCH --mem=237GB
#SBATCH -t 08:00:00
#SBATCH -J bam_extract_discord
#SBATCH -o bam_extract_discord.%j.out
#SBATCH -e bam_extract_discord.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

# =============================================================================
# extract_discordant_longspan_from_BAMminimap2.slurm
#
# Extract same-chromosome, discordant, long-span read pairs from the ORIGINAL
# BAMminimap2 column in the manifest.
#
# Manifest columns:
#   Sample  BAMminimap2  BAMP99TLENMAPQ30
#
# Keep class:
#   - paired
#   - NOT proper-pair
#   - both read and mate mapped
#   - primary
#   - not QC-fail
#   - optional duplicate removal
#   - not supplementary
#   - same chromosome (RNEXT "=" or equals RNAME)
#   - opposite strands
#   - abs(TLEN) >= MIN_TLEN
#   - MAPQ >= MIN_MAPQ
#
# Usage:
#   sbatch extract_discordant_longspan_from_BAMminimap2.slurm \
#     <sample_manifest.tsv> <outdir> <min_tlen> [min_mapq] [dedup_yes_no]
#
# Example:
#   sbatch extract_discordant_longspan_from_BAMminimap2.slurm \
#     /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv \
#     /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/pa_roary_results/discordant_longspan_from_minimap2 \
#     800 0 yes
# =============================================================================

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <sample_manifest.tsv> <outdir> <min_tlen> [min_mapq] [dedup_yes_no]"
  echo "Manifest columns: Sample BAMminimap2 BAMP99TLENMAPQ30"
  exit 1
fi

SAMTOOLS="$(command -v samtools)"

MANIFEST="$1"
OUTDIR="$2"
MIN_TLEN="$3"
MIN_MAPQ="${4:-0}"
DEDUP="${5:-yes}"   # yes = remove duplicates, no = keep duplicates

THREADS="${SLURM_NTASKS:-16}"

timestamp(){ date '+%F %T'; }

mkdir -p "$OUTDIR"/{bam,logs,flagstat}

echo "[$(timestamp)] [INFO] MANIFEST=$MANIFEST"
echo "[$(timestamp)] [INFO] OUTDIR=$OUTDIR"
echo "[$(timestamp)] [INFO] MIN_TLEN=$MIN_TLEN"
echo "[$(timestamp)] [INFO] MIN_MAPQ=$MIN_MAPQ"
echo "[$(timestamp)] [INFO] DEDUP=$DEDUP"
echo "[$(timestamp)] [INFO] THREADS=$THREADS"

[[ -s "$MANIFEST" ]] || { echo "[$(timestamp)] [ERROR] missing manifest: $MANIFEST" >&2; exit 2; }

while IFS=$'\t' read -r SAMPLE BAMMINIMAP2 BAMP99TLENMAPQ30; do
  [[ "${SAMPLE:-}" == "Sample" ]] && continue
  [[ -z "${SAMPLE:-}" ]] && continue

  BAM="$BAMMINIMAP2"
  [[ -s "$BAM" ]] || { echo "[$(timestamp)] [ERROR] missing BAM for $SAMPLE: $BAM" >&2; exit 3; }

  OUT="$OUTDIR/bam/${SAMPLE}.discordant_longspan.from_minimap2.bam"
  LOG="$OUTDIR/logs/${SAMPLE}.discordant_longspan.from_minimap2.log"
  FLAGSTAT="$OUTDIR/flagstat/${SAMPLE}.discordant_longspan.from_minimap2.flagstat.txt"

  echo "[$(timestamp)] [INFO] SAMPLE=$SAMPLE" | tee "$LOG"
  echo "[$(timestamp)] [INFO] BAM=$BAM" | tee -a "$LOG"
  echo "[$(timestamp)] [INFO] OUT=$OUT" | tee -a "$LOG"
  echo "[$(timestamp)] [INFO] THREADS=$THREADS MIN_TLEN=$MIN_TLEN MIN_MAPQ=$MIN_MAPQ DEDUP=$DEDUP" | tee -a "$LOG"

  echo "[$(timestamp)] [STEP] samtools view + awk filters" | tee -a "$LOG"

  "$SAMTOOLS" view -@ "$THREADS" -h "$BAM" \
  | awk -v min_tlen="$MIN_TLEN" -v min_mapq="$MIN_MAPQ" -v dedup="$DEDUP" '
      BEGIN{OFS="\t"}
      /^@/ {print; next}
      {
        flag=$2
        mapq=$5 + 0
        rname=$3
        rnext=$7
        tlen=$9 + 0

        paired        = and(flag, 0x1)
        proper        = and(flag, 0x2)
        unmapped      = and(flag, 0x4)
        mate_unmapped = and(flag, 0x8)
        rev           = and(flag, 0x10)
        mate_rev      = and(flag, 0x20)
        secondary     = and(flag, 0x100)
        qcfail        = and(flag, 0x200)
        duplicate     = and(flag, 0x400)
        supplementary = and(flag, 0x800)

        abs_tlen = tlen
        if (abs_tlen < 0) abs_tlen = -abs_tlen

        samechr = (rnext == "=" || rnext == rname)

        keep = 1
        if (!paired) keep=0
        if (proper) keep=0
        if (unmapped || mate_unmapped) keep=0
        if (secondary || qcfail || supplementary) keep=0
        if (!samechr) keep=0
        if (rev == mate_rev) keep=0
        if (abs_tlen < min_tlen) keep=0
        if (mapq < min_mapq) keep=0
        if (dedup == "yes" && duplicate) keep=0

        if (keep) print
      }' \
  | "$SAMTOOLS" view -@ "$THREADS" -bh -o "$OUT" - 2>>"$LOG"

  echo "[$(timestamp)] [STEP] index" | tee -a "$LOG"
  "$SAMTOOLS" index -@ "$THREADS" "$OUT" >> "$LOG" 2>&1

  echo "[$(timestamp)] [STEP] flagstat" | tee -a "$LOG"
  "$SAMTOOLS" flagstat -@ "$THREADS" "$OUT" > "$FLAGSTAT"

  {
    echo -e "retained_alignments\t$("$SAMTOOLS" view -@ "$THREADS" -c "$OUT")"
    echo -e "retained_qnames\t$("$SAMTOOLS" view -@ "$THREADS" "$OUT" | cut -f1 | sort | uniq | wc -l)"
  } | tee -a "$LOG"

  echo "[$(timestamp)] [DONE] OK" | tee -a "$LOG"
  echo | tee -a "$LOG"

done < "$MANIFEST"

echo "[$(timestamp)] [DONE] all samples completed"
