#!/usr/bin/env Rscript
# =============================================================================
# plot_roh_by_chromosome.R — Per-chromosome ROH/FROH/het plots + heatmaps
# =============================================================================
# Usage:
#   Rscript plot_roh_by_chromosome.R \
#     <per_chr_roh_summary.tsv> \
#     <out_dir> \
#     [chrom_sizes.tsv] [ancestry_labels.tsv]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: plot_roh_by_chromosome.R <per_chr_roh.tsv> <out_dir> [chrom_sizes] [ancestry]")

chr_file  <- args[1]
out_dir   <- args[2]
sizes_file <- if (length(args) >= 3 && file.exists(args[3])) args[3] else NULL
anc_file   <- if (length(args) >= 4 && file.exists(args[4])) args[4] else NULL

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# ── Read ──────────────────────────────────────────────────────────────────
d <- fread(chr_file)
# Expect: sample, chrom, roh_bp, n_tracts, longest_tract, callable_bp_chr, froh_chr

# Order chromosomes naturally
chrom_order <- unique(d$chrom)
# Natural sort: extract numeric part if possible
chrom_nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chrom_order)))
chrom_order <- chrom_order[order(ifelse(is.na(chrom_nums), 999, chrom_nums), chrom_order)]
d[, chrom := factor(chrom, levels = chrom_order)]

# Optional ancestry
if (!is.null(anc_file)) {
  anc <- fread(anc_file)
  if (ncol(anc) >= 2) {
    setnames(anc, 1:2, c("sample", "ancestry"))
    d <- merge(d, anc[, .(sample, ancestry)], by = "sample", all.x = TRUE)
  }
}
if (!"ancestry" %in% names(d)) d[, ancestry := "all"]

# ── 1. Per-chromosome FROH boxplot ───────────────────────────────────────
p1 <- ggplot(d, aes(x = chrom, y = froh_chr)) +
  geom_boxplot(outlier.size = 0.5, fill = "#E07A5F", alpha = 0.3) +
  labs(
    title = expression("Per-chromosome"~F[ROH]),
    x = "Chromosome", y = expression(F[ROH])
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
ggsave(file.path(out_dir, "per_chr_froh_boxplot.png"), p1, width = 12, height = 5, dpi = 400)
ggsave(file.path(out_dir, "per_chr_froh_boxplot.pdf"), p1, width = 12, height = 5)

# ── 2. Per-chromosome ROH burden (total bp) ─────────────────────────────
p2 <- ggplot(d, aes(x = chrom, y = roh_bp / 1e6)) +
  geom_boxplot(outlier.size = 0.5, fill = "#457B9D", alpha = 0.3) +
  labs(
    title = "Per-chromosome ROH burden",
    x = "Chromosome", y = "ROH burden (Mb)"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
ggsave(file.path(out_dir, "per_chr_roh_burden_boxplot.png"), p2, width = 12, height = 5, dpi = 400)
ggsave(file.path(out_dir, "per_chr_roh_burden_boxplot.pdf"), p2, width = 12, height = 5)

# ── 3. Heatmap: FROH by sample × chromosome ─────────────────────────────
# Order samples by total FROH
sample_total <- d[, .(total_froh = sum(roh_bp, na.rm = TRUE)), by = sample]
sample_total <- sample_total[order(total_froh)]
d[, sample := factor(sample, levels = sample_total$sample)]

p3 <- ggplot(d, aes(x = chrom, y = sample, fill = froh_chr)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", name = expression(F[ROH])) +
  labs(
    title = expression("Sample × Chromosome"~F[ROH]~"heatmap"),
    x = "Chromosome", y = "Sample"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 3)
  )
ggsave(file.path(out_dir, "heatmap_froh_sample_chr.png"), p3, width = 12, height = 14, dpi = 400)
ggsave(file.path(out_dir, "heatmap_froh_sample_chr.pdf"), p3, width = 12, height = 14)

# ── 4. Heatmap: ROH burden by sample × chromosome ───────────────────────
p4 <- ggplot(d, aes(x = chrom, y = sample, fill = roh_bp / 1e6)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = "ROH (Mb)") +
  labs(
    title = "Sample × Chromosome ROH burden heatmap",
    x = "Chromosome", y = "Sample"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 3)
  )
ggsave(file.path(out_dir, "heatmap_roh_burden_sample_chr.png"), p4, width = 12, height = 14, dpi = 400)
ggsave(file.path(out_dir, "heatmap_roh_burden_sample_chr.pdf"), p4, width = 12, height = 14)

# ── 5. ROH frequency by chromosome (fraction of samples with ROH) ───────
roh_freq <- d[, .(
  n_with_roh = sum(roh_bp > 0),
  n_total = .N,
  frac_with_roh = sum(roh_bp > 0) / .N
), by = chrom]

p5 <- ggplot(roh_freq, aes(x = chrom, y = frac_with_roh)) +
  geom_col(fill = "#3D405B", alpha = 0.7) +
  labs(
    title = "Fraction of samples with ROH per chromosome",
    x = "Chromosome", y = "Fraction with ROH"
  ) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
ggsave(file.path(out_dir, "per_chr_roh_frequency.png"), p5, width = 12, height = 5, dpi = 400)
ggsave(file.path(out_dir, "per_chr_roh_frequency.pdf"), p5, width = 12, height = 5)

cat("Chromosome-level plots written to:", out_dir, "\n")
