#!/usr/bin/env Rscript
# =============================================================================
# plot_heterozygosity_core.R — Core genome-wide heterozygosity plots
# =============================================================================
# Usage:
#   Rscript plot_heterozygosity_core.R \
#     <genomewide_heterozygosity.tsv> \
#     <out_dir> \
#     [optional: ancestry_labels.tsv]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: plot_heterozygosity_core.R <het.tsv> <out_dir> [ancestry.tsv]")
}

het_file <- args[1]
out_dir  <- args[2]
anc_file <- if (length(args) >= 3 && file.exists(args[3])) args[3] else NULL

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Read data ──────────────────────────────────────────────────────────────
het <- fread(het_file)
# Expect: sample, het_genomewide (or het_saf1)
het_col <- intersect(names(het), c("het_genomewide", "het_saf1"))[1]
if (is.na(het_col)) stop("Cannot find het_genomewide or het_saf1 column")
setnames(het, het_col, "het")
het <- het[!is.na(het) & het != "NA"]
het[, het := as.numeric(het)]

# Optional ancestry
if (!is.null(anc_file)) {
  anc <- fread(anc_file)
  if (ncol(anc) >= 2) {
    setnames(anc, 1:2, c("sample", "ancestry"))
    het <- merge(het, anc[, .(sample, ancestry)], by = "sample", all.x = TRUE)
  }
}
if (!"ancestry" %in% names(het)) het[, ancestry := "all"]

# ── Theme ──────────────────────────────────────────────────────────────────
theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom"
  )

# ── 1. Boxplot ─────────────────────────────────────────────────────────────
p1 <- ggplot(het, aes(x = ancestry, y = het)) +
  geom_boxplot(outlier.shape = NA, fill = "steelblue", alpha = 0.4) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  labs(
    title = "Genome-wide heterozygosity (ANGSD single-sample SFS)",
    x = "", y = "Heterozygosity"
  ) +
  theme_pub

ggsave(file.path(out_dir, "het_genomewide_boxplot.png"), p1, width = 7, height = 5, dpi = 400)
ggsave(file.path(out_dir, "het_genomewide_boxplot.pdf"), p1, width = 7, height = 5)

# ── 2. Violin ─────────────────────────────────────────────────────────────
p2 <- ggplot(het, aes(x = ancestry, y = het)) +
  geom_violin(fill = "steelblue", alpha = 0.3, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.1, size = 1.2, alpha = 0.6) +
  labs(
    title = "Genome-wide heterozygosity distribution",
    x = "", y = "Heterozygosity"
  ) +
  theme_pub

ggsave(file.path(out_dir, "het_genomewide_violin.png"), p2, width = 7, height = 5, dpi = 400)
ggsave(file.path(out_dir, "het_genomewide_violin.pdf"), p2, width = 7, height = 5)

# ── 3. Histogram ──────────────────────────────────────────────────────────
p3 <- ggplot(het, aes(x = het)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_vline(aes(xintercept = median(het, na.rm = TRUE)),
             linetype = "dashed", color = "red") +
  labs(
    title = "Genome-wide heterozygosity histogram",
    x = "Heterozygosity", y = "Count"
  ) +
  theme_pub

ggsave(file.path(out_dir, "het_genomewide_histogram.png"), p3, width = 7, height = 4, dpi = 400)
ggsave(file.path(out_dir, "het_genomewide_histogram.pdf"), p3, width = 7, height = 4)

# ── 4. Ranked bar plot ────────────────────────────────────────────────────
het_sorted <- het[order(het)]
het_sorted[, rank := .I]

p4 <- ggplot(het_sorted, aes(x = rank, y = het, fill = ancestry)) +
  geom_col(width = 1) +
  labs(
    title = "Per-sample heterozygosity (ranked)",
    x = "Sample rank", y = "Heterozygosity", fill = "Ancestry"
  ) +
  theme_pub

ggsave(file.path(out_dir, "het_genomewide_ranked.png"), p4, width = 10, height = 4, dpi = 400)
ggsave(file.path(out_dir, "het_genomewide_ranked.pdf"), p4, width = 10, height = 4)

cat("Heterozygosity core plots written to:", out_dir, "\n")
