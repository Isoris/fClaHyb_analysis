#!/usr/bin/env Rscript
# =============================================================================
# plot_breeding_heatmaps.R — Sample × chr/segment theta heatmaps (multiscale)
# =============================================================================
# Usage:
#   Rscript plot_breeding_heatmaps.R <breeding_dir> <out_dir> [order_file]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: plot_breeding_heatmaps.R <breeding_dir> <out_dir> [order_file]")

breed_dir  <- args[1]
out_dir    <- args[2]
order_file <- if (length(args) >= 3 && file.exists(args[3])) args[3] else NULL

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), legend.position = "bottom")

# Sample order
sample_order <- NULL
if (!is.null(order_file)) {
  sample_order <- rev(readLines(order_file))
  sample_order <- sample_order[nchar(trimws(sample_order)) > 0]
}

nat_chrom_order <- function(chroms) {
  nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chroms)))
  chroms[order(ifelse(is.na(nums), 999, nums), chroms)]
}

# ── A. Sample × chromosome theta heatmaps ──────────────────────────────
chr_files <- list.files(breed_dir, pattern = "^per_chr_theta_summary_.*\\.tsv$", full.names = TRUE)

for (fp in chr_files) {
  scale_label <- sub(".*per_chr_theta_summary_(.*)\\.tsv$", "\\1", basename(fp))
  d <- fread(fp)
  if (nrow(d) == 0) next

  d[, mean_theta := as.numeric(mean_theta)]
  chroms <- nat_chrom_order(unique(d$chrom))
  d[, chrom := factor(chrom, levels = chroms)]

  if (!is.null(sample_order)) {
    d[, sample := factor(sample, levels = sample_order)]
  } else {
    sord <- d[, .(m = mean(mean_theta, na.rm = TRUE)), by = sample][order(m)]$sample
    d[, sample := factor(sample, levels = rev(sord))]
  }

  p <- ggplot(d[!is.na(mean_theta)], aes(x = chrom, y = sample, fill = mean_theta)) +
    geom_tile() +
    scale_fill_viridis_c(option = "viridis", name = "Mean theta/site") +
    labs(title = paste0("Sample × Chromosome mean theta (", scale_label, ")"),
         subtitle = "Diversity proxy, NOT literal Hobs",
         x = "Chromosome", y = "") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 3))

  h <- max(5, length(unique(d$sample)) * 0.08 + 2)
  ggsave(file.path(out_dir, paste0("heatmap_theta_chr_", scale_label, ".png")),
         p, width = 12, height = h, dpi = 400)
  ggsave(file.path(out_dir, paste0("heatmap_theta_chr_", scale_label, ".pdf")),
         p, width = 12, height = h)
}

# ── B. Sample × segment theta heatmaps ─────────────────────────────────
seg_files <- list.files(breed_dir, pattern = "^per_sample_segment_breeding_features_.*\\.tsv$",
                        full.names = TRUE)

for (fp in seg_files) {
  label <- sub(".*per_sample_segment_breeding_features_(.*)\\.tsv$", "\\1", basename(fp))

  d <- fread(fp, select = c("sample", "chrom", "segment_id", "start", "end",
                             "mean_theta_segment", "froh_segment", "low_div_flag"))
  if (nrow(d) == 0) next

  d[, mean_theta_segment := as.numeric(mean_theta_segment)]
  d[, froh_segment := as.numeric(froh_segment)]

  # Order segments by genome position
  chroms <- nat_chrom_order(unique(d$chrom))
  d[, chrom := factor(chrom, levels = chroms)]
  d <- d[order(chrom, start)]
  d[, seg_idx := .I]

  if (!is.null(sample_order)) {
    d[, sample := factor(sample, levels = sample_order)]
  }

  # Theta heatmap
  p1 <- ggplot(d[!is.na(mean_theta_segment)],
               aes(x = seg_idx, y = sample, fill = mean_theta_segment)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_viridis_c(option = "viridis", name = "Mean theta") +
    labs(title = paste0("Sample × segment theta (", label, ")"),
         x = "Genome segment (ordered)", y = "") +
    theme_pub +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 2))

  h <- max(5, length(unique(d$sample)) * 0.06 + 2)
  w <- max(10, length(unique(d$seg_idx)) * 0.01 + 4)
  ggsave(file.path(out_dir, paste0("heatmap_theta_segments_", label, ".png")),
         p1, width = min(w, 30), height = min(h, 25), dpi = 300)

  # FROH heatmap
  p2 <- ggplot(d, aes(x = seg_idx, y = sample, fill = froh_segment)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_viridis_c(option = "inferno", name = expression(F[ROH])) +
    labs(title = paste0("Sample × segment FROH (", label, ")"),
         x = "Genome segment (ordered)", y = "") +
    theme_pub +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 2))

  ggsave(file.path(out_dir, paste0("heatmap_froh_segments_", label, ".png")),
         p2, width = min(w, 30), height = min(h, 25), dpi = 300)
}

cat("Breeding heatmaps written to:", out_dir, "\n")
