#!/usr/bin/env bash
# =============================================================================
# 05_run_all_plots.sh — Run all plotting + statistics scripts
# =============================================================================
set -euo pipefail
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPTDIR}/00_config.sh"
hr_init_dirs

hr_log "=== Step 05: Run all plots and statistics ==="

# ── Check R packages ─────────────────────────────────────────────────────
Rscript -e 'for(p in c("data.table","ggplot2")) if(!require(p, character.only=TRUE)) stop(paste("Missing R package:", p))' || {
  hr_die "Missing required R packages: data.table, ggplot2"
}

MASTER="${DIR_TABLES}/master_summary.tsv"
PER_CHR="${DIR_TABLES}/per_chr_roh_summary.tsv"
BINS="${DIR_TABLES}/catfish_roh.per_sample_roh_bins_long.tsv"
HET_TSV="${DIR_HET}/04_summary/genomewide_heterozygosity.tsv"
CHROM_SIZES="${DIR_INPUTS}/chrom_sizes.tsv"

hr_check_file "${MASTER}" "master_summary.tsv"

# ── 1. Core heterozygosity plots ────────────────────────────────────────
hr_log "Running plot_heterozygosity_core.R..."
if [[ -f "${HET_TSV}" ]]; then
  Rscript "${SCRIPTDIR}/plot_heterozygosity_core.R" \
    "${HET_TSV}" \
    "${DIR_PLOTS_CORE}" \
    ${ANCESTRY_LABELS:+"${ANCESTRY_LABELS}"} \
    2>&1 | tee "${DIR_LOGS}/plot_het_core.log"
else
  hr_log "  Skipping (no het summary TSV)"
fi

# ── 2. Core ROH plots ──────────────────────────────────────────────────
hr_log "Running plot_roh_core.R..."
if [[ -f "${BINS}" ]]; then
  Rscript "${SCRIPTDIR}/plot_roh_core.R" \
    "${MASTER}" \
    "${BINS}" \
    "${DIR_PLOTS_CORE}" \
    2>&1 | tee "${DIR_LOGS}/plot_roh_core.log"
else
  # Try without bins
  Rscript "${SCRIPTDIR}/plot_roh_core.R" \
    "${MASTER}" \
    "/dev/null" \
    "${DIR_PLOTS_CORE}" \
    2>&1 | tee "${DIR_LOGS}/plot_roh_core.log" || true
fi

# ── 3. Scatter plots ───────────────────────────────────────────────────
hr_log "Running plot_scatter_stats.R..."
Rscript "${SCRIPTDIR}/plot_scatter_stats.R" \
  "${MASTER}" \
  "${DIR_PLOTS_CORE}" \
  ${SAMPLE_QC_TABLE:+"${SAMPLE_QC_TABLE}"} \
  ${ANCESTRY_LABELS:+"${ANCESTRY_LABELS}"} \
  2>&1 | tee "${DIR_LOGS}/plot_scatter.log"

# ── 4. Chromosome-level plots ─────────────────────────────────────────
hr_log "Running plot_roh_by_chromosome.R..."
if [[ -f "${PER_CHR}" ]]; then
  Rscript "${SCRIPTDIR}/plot_roh_by_chromosome.R" \
    "${PER_CHR}" \
    "${DIR_PLOTS_CORE}" \
    ${CHROM_SIZES:+"${CHROM_SIZES}"} \
    ${ANCESTRY_LABELS:+"${ANCESTRY_LABELS}"} \
    2>&1 | tee "${DIR_LOGS}/plot_chr.log"
else
  hr_log "  Skipping (no per_chr_roh_summary.tsv)"
fi

# ── 5. Theta ideogram plots ──────────────────────────────────────────
hr_log "Running plot_theta_ideogram.R..."

if [[ -f "${CHROM_SIZES}" ]]; then
  # Main theta directory
  if [[ -d "${DIR_HET}/03_theta" ]]; then
    hr_log "  Theta plots: main directory"
    Rscript "${SCRIPTDIR}/plot_theta_ideogram.R" \
      "${DIR_HET}/03_theta" \
      "${CHROM_SIZES}" \
      "${DIR_PLOTS_CORE}/theta_main" \
      "${SAMPLE_LIST}" \
      10 \
      2>&1 | tee "${DIR_LOGS}/plot_ideogram_main.log"
  else
    hr_log "  Skipping main theta plots (directory missing)"
  fi

  # Multiscale theta directory
  if [[ -d "${DIR_HET}/03_theta/multiscale" ]]; then
    hr_log "  Theta plots: multiscale directory"
    Rscript "${SCRIPTDIR}/plot_theta_ideogram.R" \
      "${DIR_HET}/03_theta/multiscale" \
      "${CHROM_SIZES}" \
      "${DIR_PLOTS_CORE}/theta_multiscale" \
      "${SAMPLE_LIST}" \
      10 \
      2>&1 | tee "${DIR_LOGS}/plot_ideogram_multiscale.log"
  else
    hr_log "  Skipping multiscale theta plots (directory missing)"
  fi
else
  hr_log "  Skipping theta plots (chrom_sizes missing)"
fi

# ── 6. Metadata overlay plots ────────────────────────────────────────
hr_log "Running plot_roh_metadata_overlays.R..."
META_ARGS=()
if [[ -n "${ANCESTRY_LABELS}" && -f "${ANCESTRY_LABELS}" ]]; then
  META_ARGS+=(--ancestry "${ANCESTRY_LABELS}")
fi
if [[ -n "${COVTREE_ORDER}" && -f "${COVTREE_ORDER}" ]]; then
  META_ARGS+=(--order "${COVTREE_ORDER}")
fi
if [[ -n "${PRUNED81_SAMPLES}" && -f "${PRUNED81_SAMPLES}" ]]; then
  META_ARGS+=(--pruned81 "${PRUNED81_SAMPLES}")
fi
if [[ -n "${NGSRELATE_PAIRS}" && -f "${NGSRELATE_PAIRS}" ]]; then
  META_ARGS+=(--kinship "${NGSRELATE_PAIRS}")
fi

if [[ -f "${PER_CHR}" ]]; then
  Rscript "${SCRIPTDIR}/plot_roh_metadata_overlays.R" \
    "${MASTER}" \
    "${PER_CHR}" \
    "${BINS:-/dev/null}" \
    "${DIR_PLOTS_META}" \
    "${META_ARGS[@]}" \
    2>&1 | tee "${DIR_LOGS}/plot_metadata.log"
else
  hr_log "  Skipping metadata overlays (missing per_chr table)"
fi

# ── 7. Statistics ────────────────────────────────────────────────────
hr_log "Running run_stats.R..."
Rscript "${SCRIPTDIR}/run_stats.R" \
  "${MASTER}" \
  "${DIR_STATS}" \
  ${ANCESTRY_LABELS:+"${ANCESTRY_LABELS}"} \
  2>&1 | tee "${DIR_LOGS}/run_stats.log"

# ── 8. Report ────────────────────────────────────────────────────────
hr_log "Running 06_write_report.py..."
python3 "${SCRIPTDIR}/06_write_report.py" \
  --tables-dir "${DIR_TABLES}" \
  --stats-dir "${DIR_STATS}" \
  --out-dir "${DIR_REPORT}" \
  --ngsf-reps "${NGSFHMM_REPS}" \
  2>&1 | tee "${DIR_LOGS}/write_report.log"

hr_log "=== Step 05 complete ==="
hr_log "Core plots:     ${DIR_PLOTS_CORE}"
hr_log "Metadata plots: ${DIR_PLOTS_META}"
hr_log "Statistics:     ${DIR_STATS}"
hr_log "Report:         ${DIR_REPORT}"
