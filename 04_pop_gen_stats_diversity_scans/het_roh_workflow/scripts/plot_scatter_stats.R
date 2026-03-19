#!/usr/bin/env Rscript
# =============================================================================
# plot_scatter_stats.R — Scatterplots of het vs ROH/FROH/depth
# =============================================================================
# Usage:
#   Rscript plot_scatter_stats.R <master_summary.tsv> <out_dir> \
#     [sample_qc_table.tsv] [ancestry_labels.tsv]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: plot_scatter_stats.R <master.tsv> <out_dir> [qc_table] [ancestry]")

master_file <- args[1]
out_dir     <- args[2]
qc_file     <- if (length(args) >= 3 && file.exists(args[3])) args[3] else NULL
anc_file    <- if (length(args) >= 4 && file.exists(args[4])) args[4] else NULL

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# ── Read ──────────────────────────────────────────────────────────────────
d <- fread(master_file)
d[, het_genomewide := as.numeric(het_genomewide)]
d[, froh := as.numeric(froh)]
d[, roh_total_bp := as.numeric(roh_total_bp)]
d[, longest_roh := as.numeric(longest_roh)]

# Merge depth if available
if (!is.null(qc_file)) {
  qc <- fread(qc_file)
  depth_col <- intersect(names(qc), c("mean_depth", "MeanDepth", "depth", "avgDP"))
  sample_col <- intersect(names(qc), c("sample", "Sample", "SAMPLE"))
  if (length(depth_col) > 0 && length(sample_col) > 0) {
    setnames(qc, sample_col[1], "sample")
    setnames(qc, depth_col[1], "mean_depth")
    d <- merge(d, qc[, .(sample, mean_depth)], by = "sample", all.x = TRUE)
  }
}

# Merge ancestry if available
if (!is.null(anc_file)) {
  anc <- fread(anc_file)
  if (ncol(anc) >= 2) {
    setnames(anc, 1:2, c("sample", "ancestry"))
    d <- merge(d, anc[, .(sample, ancestry)], by = "sample", all.x = TRUE)
  }
}
if (!"ancestry" %in% names(d)) d[, ancestry := "all"]

# ── Helper ────────────────────────────────────────────────────────────────
make_scatter <- function(dt, x, y, xlab, ylab, title, fname, color = "ancestry") {
  p <- ggplot(dt, aes_string(x = x, y = y, color = color)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    labs(title = title, x = xlab, y = ylab) +
    theme_pub
  ggsave(file.path(out_dir, paste0(fname, ".png")), p, width = 7, height = 5, dpi = 400)
  ggsave(file.path(out_dir, paste0(fname, ".pdf")), p, width = 7, height = 5)
}

# ── Scatterplots ──────────────────────────────────────────────────────────
make_scatter(d, "froh", "het_genomewide",
  expression(F[ROH]), "Genome-wide heterozygosity",
  expression("Heterozygosity vs"~F[ROH]),
  "scatter_het_vs_froh")

make_scatter(d, "roh_total_bp", "het_genomewide",
  "Total ROH (bp)", "Genome-wide heterozygosity",
  "Heterozygosity vs total ROH",
  "scatter_het_vs_roh_total")

make_scatter(d, "longest_roh", "het_genomewide",
  "Longest ROH (bp)", "Genome-wide heterozygosity",
  "Heterozygosity vs longest ROH",
  "scatter_het_vs_longest_roh")

# Depth scatterplots if available
if ("mean_depth" %in% names(d)) {
  d[, mean_depth := as.numeric(mean_depth)]
  make_scatter(d[!is.na(mean_depth)], "mean_depth", "het_genomewide",
    "Mean depth", "Genome-wide heterozygosity",
    "Depth vs heterozygosity",
    "scatter_depth_vs_het")

  make_scatter(d[!is.na(mean_depth)], "mean_depth", "froh",
    "Mean depth", expression(F[ROH]),
    expression("Depth vs"~F[ROH]),
    "scatter_depth_vs_froh")

  make_scatter(d[!is.na(mean_depth)], "mean_depth", "roh_total_bp",
    "Mean depth", "Total ROH (bp)",
    "Depth vs total ROH",
    "scatter_depth_vs_roh_total")
}

# ── Het in vs out ROH ────────────────────────────────────────────────────
if (all(c("theta_proxy_in_roh", "theta_proxy_out_roh") %in% names(d))) {
  d[, theta_proxy_in_roh := as.numeric(theta_proxy_in_roh)]
  d[, theta_proxy_out_roh := as.numeric(theta_proxy_out_roh)]
  dsub <- d[!is.na(theta_proxy_in_roh) & !is.na(theta_proxy_out_roh)]
  if (nrow(dsub) > 0) {
    p_inout <- ggplot(dsub, aes(x = theta_proxy_out_roh, y = theta_proxy_in_roh, color = ancestry)) +
      geom_point(size = 2, alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      labs(
        title = "Diversity proxy: inside vs outside ROH",
        x = "Theta proxy (outside ROH)", y = "Theta proxy (inside ROH)"
      ) +
      theme_pub
    ggsave(file.path(out_dir, "scatter_theta_in_vs_out_roh.png"), p_inout, width = 7, height = 5, dpi = 400)
    ggsave(file.path(out_dir, "scatter_theta_in_vs_out_roh.pdf"), p_inout, width = 7, height = 5)
  }
}

cat("Scatter plots written to:", out_dir, "\n")
