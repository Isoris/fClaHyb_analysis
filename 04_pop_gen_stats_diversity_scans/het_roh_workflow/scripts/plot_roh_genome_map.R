#!/usr/bin/env Rscript
# =============================================================================
# plot_roh_genome_map.R — Panel-A-style ROH genome map + stacked bar sidebar
# =============================================================================
# Usage:
#   Rscript plot_roh_genome_map.R \
#     <per_sample_genome_roh_segments.tsv> \
#     <per_sample_roh_bins_wide_fixedBins.tsv> \
#     <per_sample_roh_bins_wide_adaptiveBins.tsv> \
#     <chrom_sizes.tsv> \
#     <out_dir> \
#     [order_file] [ancestry_labels]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Need: segments.tsv fixedBins.tsv adaptiveBins.tsv chrom_sizes.tsv out_dir")

seg_file   <- args[1]
fixed_file <- args[2]
adapt_file <- args[3]
sizes_file <- args[4]
out_dir    <- args[5]
order_file <- if (length(args) >= 6 && file.exists(args[6])) args[6] else NULL
anc_file   <- if (length(args) >= 7 && file.exists(args[7])) args[7] else NULL

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "bottom")

# ── Load data ────────────────────────────────────────────────────────────
seg <- fread(seg_file)
sizes <- fread(sizes_file, header = FALSE, col.names = c("chrom", "len"))
sizes[, cum_start := cumsum(c(0, len[-.N]))]

# Merge cumulative positions
seg <- merge(seg, sizes[, .(chrom, cum_start)], by = "chrom")
seg[, gw_start := cum_start + start]
seg[, gw_end := cum_start + end]

# Sample ordering
if (!is.null(order_file)) {
  ord <- fread(order_file, header = FALSE)$V1
  seg[, sample := factor(sample, levels = rev(ord))]
} else {
  sample_order <- seg[, .(total = sum(roh_length_bp)), by = sample][order(total)]$sample
  seg[, sample := factor(sample, levels = rev(sample_order))]
}

# Bin colors
bin_colors_fixed <- c(
  "1-2Mb" = "#A8DADC", "2-4Mb" = "#457B9D", "4-8Mb" = "#1D3557",
  "8-16Mb" = "#E63946", "gt16Mb" = "#6D071A", "lt1Mb" = "#CCCCCC"
)

# ── Helper: genome map ──────────────────────────────────────────────────
make_genome_map <- function(seg_dt, bin_colors, suffix) {
  # Chromosome boundaries
  chr_bounds <- sizes[, .(xmin = cum_start, xmax = cum_start + len, chrom = chrom)]

  p <- ggplot() +
    # Chromosome alternating background
    geom_rect(data = chr_bounds[seq(1, .N, 2)],
              aes(xmin = xmin / 1e6, xmax = xmax / 1e6, ymin = -Inf, ymax = Inf),
              fill = "grey95", alpha = 0.5) +
    # ROH segments
    geom_rect(data = seg_dt,
              aes(xmin = gw_start / 1e6, xmax = gw_end / 1e6,
                  ymin = as.numeric(sample) - 0.4, ymax = as.numeric(sample) + 0.4,
                  fill = roh_bin_fixed)) +
    scale_fill_manual(values = bin_colors, name = "ROH size class") +
    scale_y_continuous(breaks = seq_along(levels(seg_dt$sample)),
                       labels = levels(seg_dt$sample)) +
    labs(title = paste0("ROH genome map (", suffix, ")"),
         x = "Genomic position (Mb)", y = "") +
    theme_pub +
    theme(axis.text.y = element_text(size = 3))

  h <- max(5, length(levels(seg_dt$sample)) * 0.08 + 2)
  ggsave(file.path(out_dir, paste0("roh_genome_map_", suffix, ".png")),
         p, width = 16, height = h, dpi = 400)
  ggsave(file.path(out_dir, paste0("roh_genome_map_", suffix, ".pdf")),
         p, width = 16, height = h)
}

make_genome_map(seg, bin_colors_fixed, "fixedBins")

# ── Stacked bar (SROH by bin) ──────────────────────────────────────────
make_stacked_bar <- function(bins_file, suffix) {
  bins <- fread(bins_file)
  # Melt to long format
  bp_cols <- grep("^roh_bp_", names(bins), value = TRUE)
  pct_cols <- grep("^pct_genome_", names(bins), value = TRUE)

  long <- melt(bins, id.vars = "sample", measure.vars = pct_cols,
               variable.name = "bin_pct", value.name = "pct")
  long[, bin := gsub("^pct_genome_", "", bin_pct)]
  long[, pct := as.numeric(pct)]

  # Order samples by total FROH
  froh_order <- bins[, .(froh = as.numeric(froh)), by = sample][order(froh)]$sample
  long[, sample := factor(sample, levels = froh_order)]

  p <- ggplot(long, aes(x = sample, y = pct, fill = bin)) +
    geom_col(width = 0.9) +
    coord_flip() +
    labs(title = paste0("ROH size-class contribution (", suffix, ")"),
         x = "", y = "% genome in ROH", fill = "ROH bin") +
    theme_pub +
    theme(axis.text.y = element_text(size = 3))

  h <- max(5, length(unique(long$sample)) * 0.08 + 2)
  ggsave(file.path(out_dir, paste0("roh_stacked_bar_", suffix, ".png")),
         p, width = 8, height = h, dpi = 400)
  ggsave(file.path(out_dir, paste0("roh_stacked_bar_", suffix, ".pdf")),
         p, width = 8, height = h)
}

make_stacked_bar(fixed_file, "fixedBins")
make_stacked_bar(adapt_file, "adaptiveBins")

# ── SROH vs NROH scatter ────────────────────────────────────────────────
roh_summary <- seg[, .(SROH = sum(roh_length_bp), NROH = .N), by = sample]
ct <- cor.test(roh_summary$SROH, roh_summary$NROH, method = "pearson")
label_text <- sprintf("r = %.3f, p = %.2e", ct$estimate, ct$p.value)

p_scatter <- ggplot(roh_summary, aes(x = SROH / 1e6, y = NROH)) +
  geom_point(size = 2, alpha = 0.7, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.5) +
  annotate("text", x = Inf, y = Inf, label = label_text,
           hjust = 1.1, vjust = 1.5, size = 4) +
  labs(title = expression(S[ROH]~"vs"~N[ROH]),
       x = expression(S[ROH]~"(Mb)"), y = expression(N[ROH])) +
  theme_pub

ggsave(file.path(out_dir, "scatter_SROH_vs_NROH.png"), p_scatter, width = 7, height = 5, dpi = 400)
ggsave(file.path(out_dir, "scatter_SROH_vs_NROH.pdf"), p_scatter, width = 7, height = 5)

cat("ROH genome map + scatter plots written to:", out_dir, "\n")
