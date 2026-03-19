#!/usr/bin/env Rscript
# =============================================================================
# plot_roh_metadata_overlays.R — Ancestry/family/relatedness grouped plots
# =============================================================================
# Produces grouped/colorized versions of core het/ROH plots when metadata
# is available. Compatible with existing metadata ecosystem:
#   - Q matrix ancestry labels
#   - covtree ordering
#   - ngsRelate relatedness pairs
#   - pruned81 representative subset
#
# Usage:
#   Rscript plot_roh_metadata_overlays.R \
#     <master_summary.tsv> \
#     <per_chr_roh_summary.tsv> \
#     <roh_bins_long.tsv> \
#     <out_dir> \
#     [--ancestry ancestry_labels.tsv] \
#     [--order covtree_order.txt] \
#     [--pruned81 pruned81_samples.txt] \
#     [--kinship ngsRelate_pairs.tsv]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ── Parse arguments ──────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plot_roh_metadata_overlays.R <master.tsv> <per_chr.tsv> <bins.tsv> <out_dir> [--ancestry ...] [--order ...] [--pruned81 ...] [--kinship ...]")
}

master_file <- args[1]
chr_file    <- args[2]
bins_file   <- args[3]
out_dir     <- args[4]

# Parse optional named arguments
anc_file    <- NULL
order_file  <- NULL
pruned_file <- NULL
kin_file    <- NULL

i <- 5
while (i <= length(args)) {
  if (args[i] == "--ancestry" && i < length(args)) {
    anc_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--order" && i < length(args)) {
    order_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--pruned81" && i < length(args)) {
    pruned_file <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--kinship" && i < length(args)) {
    kin_file <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# ── Ancestry palette ─────────────────────────────────────────────────────
# Colorblind-friendly palette, extendable
anc_palette <- c(
  "anc1" = "#E69F00", "anc2" = "#56B4E9", "anc3" = "#009E73",
  "anc4" = "#F0E442", "anc5" = "#0072B2", "anc6" = "#D55E00",
  "Q1" = "#E69F00", "Q2" = "#56B4E9", "Q3" = "#009E73",
  "Q4" = "#F0E442", "Q5" = "#0072B2", "Q6" = "#D55E00",
  "all" = "grey50", "NA" = "grey80"
)

# ── Read data ────────────────────────────────────────────────────────────
d <- fread(master_file)
d[, het_genomewide := as.numeric(het_genomewide)]
d[, froh := as.numeric(froh)]
d[, roh_total_bp := as.numeric(roh_total_bp)]

chr_d <- fread(chr_file)

# ── Merge ancestry ───────────────────────────────────────────────────────
has_ancestry <- FALSE
if (!is.null(anc_file) && file.exists(anc_file)) {
  anc <- fread(anc_file)
  if (ncol(anc) >= 2) {
    setnames(anc, 1:2, c("sample", "ancestry"))
    d <- merge(d, anc[, .(sample, ancestry)], by = "sample", all.x = TRUE)
    chr_d <- merge(chr_d, anc[, .(sample, ancestry)], by = "sample", all.x = TRUE)
    has_ancestry <- TRUE
  }
}
if (!"ancestry" %in% names(d)) {
  d[, ancestry := "all"]
  chr_d[, ancestry := "all"]
}

# ── Merge ordering ───────────────────────────────────────────────────────
has_order <- FALSE
if (!is.null(order_file) && file.exists(order_file)) {
  ord <- fread(order_file, header = FALSE)$V1
  d[, order_idx := match(sample, ord)]
  d[is.na(order_idx), order_idx := 9999L]
  chr_d[, order_idx := match(sample, ord)]
  chr_d[is.na(order_idx), order_idx := 9999L]
  has_order <- TRUE
}

# ── Pruned 81 subset ────────────────────────────────────────────────────
has_pruned <- FALSE
pruned_samples <- NULL
if (!is.null(pruned_file) && file.exists(pruned_file)) {
  pruned_samples <- fread(pruned_file, header = FALSE)$V1
  d[, is_pruned81 := sample %in% pruned_samples]
  chr_d[, is_pruned81 := sample %in% pruned_samples]
  has_pruned <- TRUE
}

# ══════════════════════════════════════════════════════════════════════════
# ANCESTRY-GROUPED PLOTS
# ══════════════════════════════════════════════════════════════════════════
if (has_ancestry && length(unique(d$ancestry[!is.na(d$ancestry)])) > 1) {
  cat("Generating ancestry-grouped plots...\n")

  # ── Boxplots by ancestry ───────────────────────────────────────────────
  for (vname in c("het_genomewide", "froh")) {
    ylab <- if (vname == "het_genomewide") "Genome-wide heterozygosity" else expression(F[ROH])
    ttl  <- if (vname == "het_genomewide") "Heterozygosity by ancestry" else expression(F[ROH]~"by ancestry")

    p <- ggplot(d[!is.na(get(vname))], aes_string(x = "ancestry", y = vname, fill = "ancestry")) +
      geom_boxplot(outlier.shape = NA, alpha = 0.4) +
      geom_jitter(aes(color = ancestry), width = 0.15, size = 1.5, alpha = 0.6) +
      scale_fill_manual(values = anc_palette, guide = "none") +
      scale_color_manual(values = anc_palette, guide = "none") +
      labs(title = ttl, x = "Ancestry cluster", y = ylab) +
      theme_pub

    ggsave(file.path(out_dir, paste0(vname, "_by_ancestry_boxplot.png")), p, width = 8, height = 5, dpi = 400)
    ggsave(file.path(out_dir, paste0(vname, "_by_ancestry_boxplot.pdf")), p, width = 8, height = 5)
  }

  # ── Violin plots by ancestry ──────────────────────────────────────────
  for (vname in c("het_genomewide", "froh", "roh_total_bp")) {
    ylab <- switch(vname,
      "het_genomewide" = "Heterozygosity",
      "froh" = expression(F[ROH]),
      "roh_total_bp" = "Total ROH (bp)")

    p <- ggplot(d[!is.na(get(vname))], aes_string(x = "ancestry", y = vname, fill = "ancestry")) +
      geom_violin(alpha = 0.3, draw_quantiles = c(0.25, 0.5, 0.75)) +
      geom_jitter(aes(color = ancestry), width = 0.1, size = 1.2, alpha = 0.5) +
      scale_fill_manual(values = anc_palette, guide = "none") +
      scale_color_manual(values = anc_palette, guide = "none") +
      labs(title = paste(vname, "by ancestry"), x = "Ancestry cluster", y = ylab) +
      theme_pub

    ggsave(file.path(out_dir, paste0(vname, "_by_ancestry_violin.png")), p, width = 8, height = 5, dpi = 400)
    ggsave(file.path(out_dir, paste0(vname, "_by_ancestry_violin.pdf")), p, width = 8, height = 5)
  }

  # ── Scatter: het vs FROH colored by ancestry ──────────────────────────
  p_sc <- ggplot(d[!is.na(het_genomewide) & !is.na(froh)],
                 aes(x = froh, y = het_genomewide, color = ancestry)) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_manual(values = anc_palette) +
    labs(
      title = expression("Heterozygosity vs"~F[ROH]~"by ancestry"),
      x = expression(F[ROH]), y = "Genome-wide heterozygosity",
      color = "Ancestry"
    ) +
    theme_pub

  ggsave(file.path(out_dir, "scatter_het_vs_froh_ancestry.png"), p_sc, width = 8, height = 5, dpi = 400)
  ggsave(file.path(out_dir, "scatter_het_vs_froh_ancestry.pdf"), p_sc, width = 8, height = 5)

  # ── Chromosome heatmap ordered by ancestry ────────────────────────────
  chrom_order <- unique(chr_d$chrom)
  chrom_nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chrom_order)))
  chrom_order <- chrom_order[order(ifelse(is.na(chrom_nums), 999, chrom_nums), chrom_order)]
  chr_d[, chrom := factor(chrom, levels = chrom_order)]

  # Order samples by ancestry, then by total FROH within ancestry
  sample_ord <- d[, .(total_roh = sum(roh_total_bp, na.rm = TRUE)), by = .(sample, ancestry)]
  sample_ord <- sample_ord[order(ancestry, total_roh)]
  chr_d[, sample := factor(sample, levels = sample_ord$sample)]

  p_hm <- ggplot(chr_d, aes(x = chrom, y = sample, fill = froh_chr)) +
    geom_tile() +
    scale_fill_viridis_c(option = "inferno", name = expression(F[ROH])) +
    labs(
      title = expression("Sample × Chromosome"~F[ROH]~"(ordered by ancestry)"),
      x = "Chromosome", y = "Sample"
    ) +
    theme_pub +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 3)
    )

  ggsave(file.path(out_dir, "heatmap_froh_by_ancestry.png"), p_hm, width = 12, height = 14, dpi = 400)
  ggsave(file.path(out_dir, "heatmap_froh_by_ancestry.pdf"), p_hm, width = 12, height = 14)

  # ── Per-chromosome FROH by ancestry ───────────────────────────────────
  p_chr_anc <- ggplot(chr_d, aes(x = chrom, y = froh_chr, fill = ancestry)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.4, position = position_dodge(0.8)) +
    scale_fill_manual(values = anc_palette) +
    labs(
      title = expression("Per-chromosome"~F[ROH]~"by ancestry"),
      x = "Chromosome", y = expression(F[ROH]), fill = "Ancestry"
    ) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  ggsave(file.path(out_dir, "per_chr_froh_by_ancestry.png"), p_chr_anc, width = 14, height = 5, dpi = 400)
  ggsave(file.path(out_dir, "per_chr_froh_by_ancestry.pdf"), p_chr_anc, width = 14, height = 5)
}

# ══════════════════════════════════════════════════════════════════════════
# PRUNED-81 SUBSET PLOTS (cleaner representative views)
# ══════════════════════════════════════════════════════════════════════════
if (has_pruned && !is.null(pruned_samples) && length(pruned_samples) > 0) {
  cat("Generating pruned-81 subset plots...\n")
  d81 <- d[is_pruned81 == TRUE]
  chr81 <- chr_d[is_pruned81 == TRUE]

  if (nrow(d81) > 5) {
    # Heatmap for pruned subset
    chrom_order <- unique(chr81$chrom)
    chrom_nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chrom_order)))
    chrom_order <- chrom_order[order(ifelse(is.na(chrom_nums), 999, chrom_nums), chrom_order)]
    chr81[, chrom := factor(chrom, levels = chrom_order)]

    sample_ord81 <- d81[, .(total_roh = sum(roh_total_bp, na.rm = TRUE)), by = sample]
    sample_ord81 <- sample_ord81[order(total_roh)]
    chr81[, sample := factor(sample, levels = sample_ord81$sample)]

    p81 <- ggplot(chr81, aes(x = chrom, y = sample, fill = froh_chr)) +
      geom_tile() +
      scale_fill_viridis_c(option = "inferno", name = expression(F[ROH])) +
      labs(
        title = expression("Pruned-81 representative subset: Sample × Chromosome"~F[ROH]),
        x = "Chromosome", y = "Sample"
      ) +
      theme_pub +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 5)
      )

    ggsave(file.path(out_dir, "heatmap_froh_pruned81.png"), p81, width = 12, height = 10, dpi = 400)
    ggsave(file.path(out_dir, "heatmap_froh_pruned81.pdf"), p81, width = 12, height = 10)

    # Boxplots for pruned subset
    p81b <- ggplot(d81, aes(x = ancestry, y = froh, fill = ancestry)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.4) +
      geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
      scale_fill_manual(values = anc_palette, guide = "none") +
      labs(
        title = expression(F[ROH]~"(pruned-81 unrelated subset)"),
        x = "Ancestry", y = expression(F[ROH])
      ) +
      theme_pub

    ggsave(file.path(out_dir, "froh_pruned81_boxplot.png"), p81b, width = 7, height = 5, dpi = 400)
    ggsave(file.path(out_dir, "froh_pruned81_boxplot.pdf"), p81b, width = 7, height = 5)
  }
}

# ══════════════════════════════════════════════════════════════════════════
# ROH BINS STACKED BAR (ancestry-colored)
# ══════════════════════════════════════════════════════════════════════════
if (file.exists(bins_file) && has_ancestry) {
  bins <- fread(bins_file)
  if (all(c("bin", "pct_genome_in_bin", "ancestry_label") %in% names(bins))) {
    bins[, pct_genome_in_bin := as.numeric(pct_genome_in_bin)]

    # Order by ancestry then total
    sample_total <- bins[, .(total = sum(pct_genome_in_bin, na.rm = TRUE)), by = .(sample, ancestry_label)]
    sample_total <- sample_total[order(ancestry_label, total)]
    bins[, sample := factor(sample, levels = sample_total$sample)]

    bin_colors <- c("1-2Mb" = "#A8DADC", "2-4Mb" = "#457B9D",
                    "4-8Mb" = "#1D3557", "8-16Mb" = "#E63946", ">16Mb" = "#6D071A")

    p_bins <- ggplot(bins, aes(x = sample, y = pct_genome_in_bin, fill = bin)) +
      geom_col(width = 0.9) +
      scale_fill_manual(values = bin_colors, name = "ROH size class") +
      facet_grid(. ~ ancestry_label, scales = "free_x", space = "free_x") +
      labs(
        title = "ROH length-class contributions by ancestry",
        x = "Sample", y = "% genome in ROH"
      ) +
      theme_pub +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))

    ggsave(file.path(out_dir, "roh_bins_stacked_by_ancestry.png"), p_bins, width = 16, height = 5, dpi = 400)
    ggsave(file.path(out_dir, "roh_bins_stacked_by_ancestry.pdf"), p_bins, width = 16, height = 5)
  }
}

cat("Metadata overlay plots written to:", out_dir, "\n")
