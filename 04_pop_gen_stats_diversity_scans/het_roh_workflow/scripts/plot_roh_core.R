#!/usr/bin/env Rscript
# =============================================================================
# plot_roh_core.R — Core ROH/FROH genome-wide plots
# =============================================================================
# Usage:
#   Rscript plot_roh_core.R <master_summary.tsv> <roh_bins_long.tsv> <out_dir>
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: plot_roh_core.R <master_summary.tsv> <roh_bins_long.tsv> <out_dir>")
}

master_file <- args[1]
bins_file   <- args[2]
out_dir     <- args[3]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# ── Read data ──────────────────────────────────────────────────────────────
d <- fread(master_file)
d[, het_genomewide := as.numeric(het_genomewide)]
d[, froh := as.numeric(froh)]

# ── 1. FROH boxplot & violin ─────────────────────────────────────────────
p1 <- ggplot(d, aes(x = "All", y = froh)) +
  geom_boxplot(outlier.shape = NA, fill = "#E07A5F", alpha = 0.4) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  labs(title = expression(F[ROH]~"distribution"), x = "", y = expression(F[ROH])) +
  theme_pub
ggsave(file.path(out_dir, "froh_boxplot.png"), p1, width = 5, height = 5, dpi = 400)
ggsave(file.path(out_dir, "froh_boxplot.pdf"), p1, width = 5, height = 5)

p1v <- ggplot(d, aes(x = "All", y = froh)) +
  geom_violin(fill = "#E07A5F", alpha = 0.3, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.1, size = 1.2, alpha = 0.6) +
  labs(title = expression(F[ROH]~"distribution"), x = "", y = expression(F[ROH])) +
  theme_pub
ggsave(file.path(out_dir, "froh_violin.png"), p1v, width = 5, height = 5, dpi = 400)
ggsave(file.path(out_dir, "froh_violin.pdf"), p1v, width = 5, height = 5)

# ── 2. Longest ROH boxplot ──────────────────────────────────────────────
p2 <- ggplot(d, aes(x = "All", y = longest_roh / 1e6)) +
  geom_boxplot(outlier.shape = NA, fill = "#81B29A", alpha = 0.4) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  labs(title = "Longest ROH tract", x = "", y = "Longest ROH (Mb)") +
  theme_pub
ggsave(file.path(out_dir, "longest_roh_boxplot.png"), p2, width = 5, height = 5, dpi = 400)
ggsave(file.path(out_dir, "longest_roh_boxplot.pdf"), p2, width = 5, height = 5)

# ── 3. ROH count boxplot ────────────────────────────────────────────────
p3 <- ggplot(d, aes(x = "All", y = n_tracts)) +
  geom_boxplot(outlier.shape = NA, fill = "#3D405B", alpha = 0.3) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6, color = "#3D405B") +
  labs(title = "Number of ROH tracts per sample", x = "", y = "Number of ROH") +
  theme_pub
ggsave(file.path(out_dir, "roh_count_boxplot.png"), p3, width = 5, height = 5, dpi = 400)
ggsave(file.path(out_dir, "roh_count_boxplot.pdf"), p3, width = 5, height = 5)

# ── 4. Total ROH bp boxplot ─────────────────────────────────────────────
p4 <- ggplot(d, aes(x = "All", y = roh_total_bp / 1e6)) +
  geom_boxplot(outlier.shape = NA, fill = "#F2CC8F", alpha = 0.5) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  labs(title = "Total ROH length per sample", x = "", y = "Total ROH (Mb)") +
  theme_pub
ggsave(file.path(out_dir, "roh_total_bp_boxplot.png"), p4, width = 5, height = 5, dpi = 400)
ggsave(file.path(out_dir, "roh_total_bp_boxplot.pdf"), p4, width = 5, height = 5)

# ── 5. ROH bins stacked bar ─────────────────────────────────────────────
if (file.exists(bins_file)) {
  bins <- fread(bins_file)
  # Expect: sample, order_index, ancestry_label, bin, roh_bp, pct_genome_in_bin
  if ("bin" %in% names(bins) && "pct_genome_in_bin" %in% names(bins)) {
    bins[, pct_genome_in_bin := as.numeric(pct_genome_in_bin)]
    # Order samples by total FROH
    sample_order <- bins[, .(total = sum(pct_genome_in_bin, na.rm = TRUE)), by = sample]
    sample_order <- sample_order[order(total)]
    bins[, sample := factor(sample, levels = sample_order$sample)]

    bin_colors <- c("1-2Mb" = "#A8DADC", "2-4Mb" = "#457B9D",
                    "4-8Mb" = "#1D3557", "8-16Mb" = "#E63946", ">16Mb" = "#6D071A")

    p5 <- ggplot(bins, aes(x = sample, y = pct_genome_in_bin, fill = bin)) +
      geom_col(width = 0.9) +
      scale_fill_manual(values = bin_colors, name = "ROH size class") +
      labs(
        title = "ROH length-class contributions per sample",
        x = "Sample (ranked by total FROH)", y = "% genome in ROH"
      ) +
      theme_pub +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))

    ggsave(file.path(out_dir, "roh_bins_stacked.png"), p5, width = 14, height = 5, dpi = 400)
    ggsave(file.path(out_dir, "roh_bins_stacked.pdf"), p5, width = 14, height = 5)

    # ── ROH length-class distribution (violin per bin) ───────────────────
    p6 <- ggplot(bins[pct_genome_in_bin > 0], aes(x = bin, y = pct_genome_in_bin)) +
      geom_boxplot(fill = "grey80") +
      geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
      labs(
        title = "Distribution of ROH size classes across samples",
        x = "ROH size class", y = "% genome in class"
      ) +
      theme_pub

    ggsave(file.path(out_dir, "roh_bins_distribution.png"), p6, width = 7, height = 5, dpi = 400)
    ggsave(file.path(out_dir, "roh_bins_distribution.pdf"), p6, width = 7, height = 5)
  }
}

cat("ROH core plots written to:", out_dir, "\n")
