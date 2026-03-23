#!/usr/bin/env Rscript

# =============================================================================
# STEP13_make_combined_panels_ABC.R
#
# Combines STEP12 outputs into a single manuscript-style multi-panel figure:
#   (a) plain regional PCA
#   (b) regional PCA colored by regional heterozygosity, with representatives
#   (c) chromosome-wide tP diversity curves by inferred group
#
# Usage:
#   Rscript STEP13_make_combined_panels_ABC.R \
#     <step12_prefix> <candidate_id> <outprefix> [width=8] [height=13]
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(paste(
    "Usage: Rscript STEP13_make_combined_panels_ABC.R",
    "<step12_prefix> <candidate_id> <outprefix>",
    "[width=8] [height=13]"
  ))
}

step12_prefix <- args[1]
candidate_id  <- args[2]
outprefix     <- args[3]
fig_width     <- if (length(args) >= 4) as.numeric(args[4]) else 8
fig_height    <- if (length(args) >= 5) as.numeric(args[5]) else 13

# ── Read STEP12 outputs ───────────────────────────────────────────────────
pcs_file   <- paste0(step12_prefix, ".candidate_", candidate_id, ".regional_pca_samples.tsv.gz")
plotc_file <- paste0(step12_prefix, ".candidate_", candidate_id, ".plotC_summary.tsv.gz")
rep_file   <- paste0(step12_prefix, ".candidate_", candidate_id, ".representatives.tsv.gz")

for (f in c(pcs_file, plotc_file, rep_file)) {
  if (!file.exists(f)) stop("Missing input file: ", f)
}

pcs      <- fread(pcs_file)
plotc_dt <- fread(plotc_file)
reps     <- fread(rep_file)

# ── Validate columns ──────────────────────────────────────────────────────
req_pcs <- c("sample", "PC1", "PC2", "ordered_group", "regional_het", "group_label")
miss_pcs <- setdiff(req_pcs, names(pcs))
if (length(miss_pcs) > 0) {
  stop("regional_pca_samples file missing required columns: ", paste(miss_pcs, collapse = ", "))
}

req_plotc <- c("group_label", "ordered_group", "WinCenter", "tP_mean")
miss_plotc <- setdiff(req_plotc, names(plotc_dt))
if (length(miss_plotc) > 0) {
  stop("plotC_summary file missing required columns: ", paste(miss_plotc, collapse = ", "))
}

# ── Recover region coordinates from reps table ───────────────────────────
# STEP12 carries candidate_start_bp / candidate_end_bp / candidate_chrom through reps
region_start <- NA_real_
region_end   <- NA_real_
region_chrom <- ""

if ("candidate_start_bp" %in% names(reps) && "candidate_end_bp" %in% names(reps)) {
  region_start <- unique(reps$candidate_start_bp)[1]
  region_end   <- unique(reps$candidate_end_bp)[1]
}
if ("candidate_chrom" %in% names(reps)) {
  region_chrom <- unique(reps$candidate_chrom)[1]
}

# ── Ensure group_label is a factor in correct order ──────────────────────
group_levels <- unique(pcs[order(ordered_group)]$group_label)
pcs[, group_label := factor(group_label, levels = group_levels)]
plotc_dt[, group_label := factor(group_label, levels = group_levels)]
reps[, group_label := factor(group_label, levels = group_levels)]

# ── Group colors ─────────────────────────────────────────────────────────
# Use a colorblind-friendly palette
n_groups <- length(group_levels)
if (n_groups == 3) {
  group_colors <- c("Homo_1" = "#0072B2", "Het" = "#D55E00", "Homo_2" = "#009E73")
  # If labels aren't exactly these, use generic
  if (!all(group_levels %in% names(group_colors))) {
    group_colors <- setNames(c("#0072B2", "#D55E00", "#009E73"), group_levels)
  }
} else {
  group_colors <- setNames(
    scales::hue_pal()(n_groups),
    group_levels
  )
}

# ── Panel (a): plain regional PCA ────────────────────────────────────────
p_a <- ggplot(pcs, aes(x = PC1, y = PC2)) +
  geom_point(shape = 1, alpha = 0.6, size = 2) +
  theme_bw(base_size = 11) +
  labs(
    title = "(a)  Regional PCA",
    subtitle = if (nchar(region_chrom) > 0) {
      paste0(region_chrom, ": ",
             format(region_start, big.mark = ","), " – ",
             format(region_end, big.mark = ","))
    } else NULL,
    x = "PC1",
    y = "PC2"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 13))

# ── Panel (b): PCA colored by regional het, with group labels ────────────
p_b <- ggplot(pcs, aes(x = PC1, y = PC2, color = regional_het)) +
  geom_point(aes(shape = group_label), size = 2.5, alpha = 0.8) +
  geom_point(
    data = reps,
    aes(x = PC1, y = PC2),
    inherit.aes = FALSE,
    shape = 1, size = 6, stroke = 1.2, color = "black"
  ) +
  geom_text(
    data = reps,
    aes(x = PC1, y = PC2, label = group_label),
    inherit.aes = FALSE,
    nudge_y = diff(range(pcs$PC2, na.rm = TRUE)) * 0.07,
    color = "black", size = 3.2, fontface = "bold"
  ) +
  theme_bw(base_size = 11) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "(b)  Regional PCA — heterozygosity",
    subtitle = "Circled = representative nearest group center",
    x = "PC1",
    y = "PC2",
    color = "Dosage\nhet",
    shape = "Group"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 13))

# ── Panel (c): chromosome-wide tP by group ──────────────────────────────
p_c <- ggplot(plotc_dt, aes(x = WinCenter / 1e6, y = tP_mean, color = group_label)) +
  {
    if (is.finite(region_start) && is.finite(region_end)) {
      annotate(
        "rect",
        xmin = region_start / 1e6, xmax = region_end / 1e6,
        ymin = -Inf, ymax = Inf,
        alpha = 0.15, fill = "grey60"
      )
    }
  } +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = group_colors) +
  theme_bw(base_size = 11) +
  labs(
    title = paste0("(c)  ", region_chrom,
                   " — pairwise θ (tP) by inferred group"),
    subtitle = paste0("Shaded = candidate ", candidate_id),
    x = "Chromosome position (Mb)",
    y = "Mean tP (pairwise θ)",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 13),
    legend.position = "bottom"
  )

# ── Combine panels ───────────────────────────────────────────────────────
combined <- p_a / p_b / p_c + plot_layout(heights = c(1, 1, 0.85))

pdf_out <- paste0(outprefix, ".candidate_", candidate_id, ".combined_ABC.pdf")
png_out <- paste0(outprefix, ".candidate_", candidate_id, ".combined_ABC.png")

ggsave(pdf_out, combined, width = fig_width, height = fig_height)
ggsave(png_out, combined, width = fig_width, height = fig_height, dpi = 400)

message("[DONE] Wrote:")
message("  ", pdf_out)
message("  ", png_out)
