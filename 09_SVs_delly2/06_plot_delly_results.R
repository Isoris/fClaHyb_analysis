#!/usr/bin/env Rscript
# =============================================================================
# 06_plot_delly_results.R
# Publication-style DEL catalog summary plots
# Modeled on Chinese alligator Fig. 2 (Wan et al.)
#
# Panels:
#   A — Genome-wide DEL density heatmap (1-Mb windows, per sample)
#       with cumulative barplot on top and per-sample total on right
#   B — Per-sample DEL count barplot (sorted)
#   C — PCA on binary DEL genotype matrix
#   D — Pairwise shared DEL heatmap with hierarchical clustering
#   E — DEL size distribution (log-scale histogram)
#   F — Per-chromosome DEL count barplot
#
# Usage:
#   Rscript 06_plot_delly_results.R \
#     --summary_dir DIR --final_dir DIR --plot_dir DIR \
#     --ref_fai PATH --samples_unrelated PATH
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
  library(stats)
})

# ── Parse arguments ─────────────────────────────────────────────────────────
option_list <- list(
  make_option("--summary_dir", type="character"),
  make_option("--final_dir",   type="character"),
  make_option("--plot_dir",    type="character"),
  make_option("--ref_fai",     type="character"),
  make_option("--samples_unrelated", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$plot_dir, showWarnings=FALSE, recursive=TRUE)

# ── Theme ───────────────────────────────────────────────────────────────────
theme_pub <- theme_classic(base_size=10) +
  theme(
    plot.title = element_text(face="bold", size=12),
    axis.text = element_text(color="black"),
    legend.position = "right",
    panel.grid = element_blank()
  )

# ── Colors ──────────────────────────────────────────────────────────────────
# Matching the teal/orange palette from the alligator figure
col_del   <- "#E69F00"  # orange for DEL
col_main  <- "#2E8B8B"  # teal
col_light <- "#7ECACA"
col_dark  <- "#1A5252"

cat("=== Loading data ===\n")

# ── Load data ───────────────────────────────────────────────────────────────
per_sample   <- fread(file.path(opt$summary_dir, "per_sample_DEL_counts.tsv"))
svlen_dist   <- fread(file.path(opt$summary_dir, "DEL_svlen_distribution.tsv"))
per_chr      <- fread(file.path(opt$summary_dir, "per_chromosome_DEL_counts.tsv"))
window_counts <- fread(file.path(opt$summary_dir, "DEL_window_counts_1Mb.tsv"))
pairwise     <- fread(file.path(opt$summary_dir, "pairwise_shared_DEL.tsv"))
binary_gt    <- fread(file.path(opt$summary_dir, "DEL_binary_genotype_matrix.tsv"))

# Reference FAI for chromosome lengths
fai <- fread(opt$ref_fai, header=FALSE, col.names=c("chrom","length","offset","line_bases","line_bytes"))

# Unrelated sample list (for coloring)
unrelated_samples <- character(0)
if (file.exists(opt$samples_unrelated)) {
  unrelated_samples <- scan(opt$samples_unrelated, what="character", quiet=TRUE)
}

cat("  Samples:", nrow(per_sample), "\n")
cat("  DEL sites:", nrow(binary_gt), "\n")
cat("  Chromosomes:", nrow(fai), "\n")

# =============================================================================
# PANEL A: Genome-wide DEL density heatmap per sample
# =============================================================================
cat("=== Panel A: Genome-wide heatmap ===\n")

# We need per-sample per-window counts from the GT matrix
# This requires re-parsing the GT matrix with position info
gt_matrix <- fread(file.path(opt$final_dir, "catalog_226.DEL.GT_matrix.tsv"))
samples <- colnames(gt_matrix)[6:ncol(gt_matrix)]

# Parse chromosome and position
gt_matrix[, window_1Mb := floor(POS / 1e6) * 1e6]

# For each sample, count DELs per 1-Mb window
cat("  Building per-sample window matrix...\n")

# Build cumulative chromosome offsets for genome-wide x-axis
chr_order <- fai$chrom
fai_dt <- data.table(chrom=fai$chrom, length=fai$length)
fai_dt[, cum_offset := cumsum(c(0, length[-.N]))]
fai_dt[, cum_end := cum_offset + length]

# Per-sample per-window matrix
window_sample_list <- list()
for (s in samples) {
  # Get DELs carried by this sample
  carried <- gt_matrix[get(s) %in% c("0/1","1/1","0|1","1|0","1|1"), 
                        .(CHROM, window_1Mb)]
  if (nrow(carried) > 0) {
    cts <- carried[, .N, by=.(CHROM, window_1Mb)]
    cts[, sample := s]
    window_sample_list[[s]] <- cts
  }
}
window_sample_dt <- rbindlist(window_sample_list)
setnames(window_sample_dt, "N", "n_DELs")

# Merge with cumulative offset
window_sample_dt <- merge(window_sample_dt, fai_dt[, .(chrom, cum_offset)], 
                          by.x="CHROM", by.y="chrom", all.x=TRUE)
window_sample_dt[, cum_pos := cum_offset + window_1Mb + 5e5]  # midpoint

# Sort samples by total DEL count
sample_order <- per_sample[order(-n_DELs), sample]

# Build wide matrix for ComplexHeatmap
# Create all windows
all_windows <- fai_dt[, .(window_1Mb = seq(0, length - 1, by=1e6)), by=chrom]
all_windows <- merge(all_windows, fai_dt[, .(chrom, cum_offset)], by="chrom")
all_windows[, window_id := paste0(chrom, ":", window_1Mb)]

# Wide matrix: samples (rows) x windows (cols)
wide_dt <- dcast(window_sample_dt, sample ~ paste0(CHROM, ":", window_1Mb),
                 value.var="n_DELs", fill=0)
# Align to sample order
wide_dt <- wide_dt[match(sample_order, sample)]

mat_a <- as.matrix(wide_dt[, -1, with=FALSE])
rownames(mat_a) <- wide_dt$sample
# Align columns to genomic order
col_order <- all_windows[order(cum_offset, window_1Mb), window_id]
col_order <- col_order[col_order %in% colnames(mat_a)]
mat_a <- mat_a[, col_order, drop=FALSE]

# Cumulative barplot (sum across samples per window)
cum_counts <- colSums(mat_a)

# Chromosome boundary positions (for vertical lines)
chr_boundaries <- all_windows[, .(first_idx = .I[1]), by=chrom]

cat("  Drawing Panel A heatmap...\n")

# Use ComplexHeatmap for the genome-wide view
col_fun_a <- colorRamp2(c(0, 5, 10, 15), c("#440154", "#31688e", "#35b779", "#fde725"))

pdf(file.path(opt$plot_dir, "panel_A_genome_heatmap.pdf"), width=18, height=10)

# Top annotation: cumulative barplot
top_anno <- HeatmapAnnotation(
  cumulative = anno_barplot(cum_counts, height=unit(2, "cm"),
                            gp=gpar(fill=col_dark, col=NA),
                            axis_param=list(gp=gpar(fontsize=6))),
  annotation_name_side="left"
)

# Right annotation: per-sample total
right_anno <- rowAnnotation(
  total = anno_barplot(per_sample[match(sample_order, sample), n_DELs],
                       width=unit(2, "cm"),
                       gp=gpar(fill=col_main, col=NA),
                       axis_param=list(gp=gpar(fontsize=6))),
  annotation_name_side="top"
)

ht <- Heatmap(mat_a,
  name="# DEL",
  col=col_fun_a,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_row_names=(nrow(mat_a) <= 80),
  show_column_names=FALSE,
  row_names_gp=gpar(fontsize=5),
  top_annotation=top_anno,
  right_annotation=right_anno,
  column_title="Genomic position (1-Mb windows)",
  row_title="Samples",
  use_raster=TRUE,
  raster_quality=3,
  heatmap_legend_param=list(title="# DEL", title_gp=gpar(fontsize=8))
)
draw(ht)
dev.off()

cat("  Panel A saved.\n")

# =============================================================================
# PANEL B: Per-sample DEL count barplot (sorted, colored by source)
# =============================================================================
cat("=== Panel B: Per-sample barplot ===\n")

per_sample[, source := ifelse(sample %in% unrelated_samples, "Unrelated (81)", "Related")]
per_sample[, sample := factor(sample, levels=sample_order)]

p_b <- ggplot(per_sample, aes(x=sample, y=n_DELs, fill=source)) +
  geom_col(width=0.8) +
  scale_fill_manual(values=c("Unrelated (81)"=col_main, "Related"=col_light),
                    name="Subset") +
  labs(x="Sample", y="DEL count", title="B") +
  theme_pub +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=3),
        plot.title = element_text(face="bold", size=14))

ggsave(file.path(opt$plot_dir, "panel_B_per_sample_DEL.pdf"), p_b,
       width=16, height=5, dpi=400)
ggsave(file.path(opt$plot_dir, "panel_B_per_sample_DEL.png"), p_b,
       width=16, height=5, dpi=400)

cat("  Panel B saved.\n")

# =============================================================================
# PANEL C: PCA on binary DEL genotype matrix
# =============================================================================
cat("=== Panel C: PCA ===\n")

# Transpose: rows=samples, cols=DELs
gt_for_pca <- as.matrix(binary_gt[, -1, with=FALSE])
rownames(gt_for_pca) <- NULL  # DELs as rows
# Transpose so samples are rows
pca_input <- t(gt_for_pca)
rownames(pca_input) <- colnames(binary_gt)[-1]

# Remove zero-variance columns
col_var <- apply(pca_input, 2, var)
pca_input <- pca_input[, col_var > 0, drop=FALSE]

cat("  PCA input:", nrow(pca_input), "samples x", ncol(pca_input), "DELs\n")

pca_res <- prcomp(pca_input, center=TRUE, scale.=FALSE)
pca_df <- data.table(
  sample = rownames(pca_input),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2]
)
pca_df[, source := ifelse(sample %in% unrelated_samples, "Unrelated (81)", "Related")]

var_pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

p_c <- ggplot(pca_df, aes(x=PC1, y=PC2, color=source)) +
  geom_point(size=2.5, alpha=0.8) +
  scale_color_manual(values=c("Unrelated (81)"=col_main, "Related"=col_light),
                     name="Subset") +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  geom_vline(xintercept=0, linetype="dashed", color="grey50") +
  labs(x=paste0("PC1 (", var_pct[1], "%)"),
       y=paste0("PC2 (", var_pct[2], "%)"),
       title="C") +
  theme_pub +
  theme(plot.title = element_text(face="bold", size=14))

ggsave(file.path(opt$plot_dir, "panel_C_PCA_DEL.pdf"), p_c,
       width=7, height=6, dpi=400)
ggsave(file.path(opt$plot_dir, "panel_C_PCA_DEL.png"), p_c,
       width=7, height=6, dpi=400)

cat("  Panel C saved.\n")

# =============================================================================
# PANEL D: Pairwise shared DEL heatmap
# =============================================================================
cat("=== Panel D: Pairwise shared heatmap ===\n")

# Build symmetric matrix
pw_wide <- dcast(pairwise, sample_i ~ sample_j, value.var="n_shared", fill=0)
pw_samples <- as.character(pw_wide$sample_i)
pw_mat <- as.matrix(pw_wide[, -1, with=FALSE])
rownames(pw_mat) <- pw_samples
# Fill symmetric
for (i in seq_along(pw_samples)) {
  for (j in seq_along(pw_samples)) {
    if (pw_mat[i,j] == 0 && pw_mat[j,i] > 0) pw_mat[i,j] <- pw_mat[j,i]
  }
}

# Hierarchical clustering
hc <- hclust(as.dist(max(pw_mat) - pw_mat), method="ward.D2")

# Color function
col_fun_d <- colorRamp2(
  quantile(pw_mat[upper.tri(pw_mat)], c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE),
  c("#440154", "#31688e", "#21918c", "#5ec962", "#fde725")
)

# Source annotation
source_colors <- c("Unrelated (81)"=col_main, "Related"=col_light)
sample_source <- ifelse(pw_samples %in% unrelated_samples, "Unrelated (81)", "Related")

row_anno <- rowAnnotation(
  Source = sample_source,
  col = list(Source = source_colors),
  show_annotation_name = FALSE
)

# Top annotation: per-sample DEL count
top_sv <- per_sample[match(pw_samples, sample), n_DELs]
top_anno_d <- HeatmapAnnotation(
  `# DELs` = anno_barplot(top_sv, height=unit(1.5, "cm"),
                           gp=gpar(fill=col_dark, col=NA)),
  annotation_name_side = "left"
)

pdf(file.path(opt$plot_dir, "panel_D_pairwise_shared.pdf"), width=12, height=11)

ht_d <- Heatmap(pw_mat,
  name="# Shared\nDEL",
  col=col_fun_d,
  cluster_rows=hc,
  cluster_columns=hc,
  show_row_names=(nrow(pw_mat) <= 80),
  show_column_names=(ncol(pw_mat) <= 80),
  row_names_gp=gpar(fontsize=4),
  column_names_gp=gpar(fontsize=4),
  left_annotation=row_anno,
  top_annotation=top_anno_d,
  column_title="D  Pairwise shared DELs",
  column_title_gp=gpar(fontsize=12, fontface="bold"),
  use_raster=TRUE,
  raster_quality=3
)
draw(ht_d)
dev.off()

cat("  Panel D saved.\n")

# =============================================================================
# PANEL E: DEL size distribution
# =============================================================================
cat("=== Panel E: DEL size distribution ===\n")

svlen_dist[, svlen_kb := svlen_bp / 1000]

p_e <- ggplot(svlen_dist, aes(x=svlen_bp)) +
  geom_histogram(bins=80, fill=col_main, color="white", linewidth=0.2) +
  scale_x_log10(labels=scales::comma) +
  annotation_logticks(sides="b") +
  labs(x="DEL size (bp, log scale)", y="Count", title="E") +
  theme_pub +
  theme(plot.title = element_text(face="bold", size=14))

ggsave(file.path(opt$plot_dir, "panel_E_DEL_size_distribution.pdf"), p_e,
       width=7, height=5, dpi=400)
ggsave(file.path(opt$plot_dir, "panel_E_DEL_size_distribution.png"), p_e,
       width=7, height=5, dpi=400)

cat("  Panel E saved.\n")

# =============================================================================
# PANEL F: Per-chromosome DEL counts
# =============================================================================
cat("=== Panel F: Per-chromosome ===\n")

# Order chromosomes naturally
per_chr[, chrom_num := as.integer(gsub("\\D+", "", chrom))]
per_chr <- per_chr[order(chrom_num)]
per_chr[, chrom_label := gsub("C_gar_", "", chrom)]
per_chr[, chrom_label := factor(chrom_label, levels=chrom_label)]

p_f <- ggplot(per_chr, aes(x=chrom_label, y=n_DELs)) +
  geom_col(fill=col_main, width=0.7) +
  labs(x="Chromosome", y="DEL count", title="F") +
  theme_pub +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(face="bold", size=14))

ggsave(file.path(opt$plot_dir, "panel_F_per_chromosome.pdf"), p_f,
       width=10, height=5, dpi=400)
ggsave(file.path(opt$plot_dir, "panel_F_per_chromosome.png"), p_f,
       width=10, height=5, dpi=400)

cat("  Panel F saved.\n")

# =============================================================================
# COMPOSITE: Combine panels B, C, E, F into one figure
# =============================================================================
cat("=== Composite figure ===\n")

# 2x2 layout: B | C / E | F
composite <- plot_grid(
  p_b + theme(legend.position="none"),
  p_c + theme(legend.position="none"),
  p_e,
  p_f,
  labels=c("B","C","E","F"),
  ncol=2, nrow=2,
  rel_widths=c(1.5, 1),
  rel_heights=c(1, 1)
)

ggsave(file.path(opt$plot_dir, "composite_BCEF.pdf"), composite,
       width=18, height=10, dpi=400)
ggsave(file.path(opt$plot_dir, "composite_BCEF.png"), composite,
       width=18, height=10, dpi=400)

cat("  Composite saved.\n")

# =============================================================================
# Summary stats for figure legend
# =============================================================================
cat("\n=== Figure summary stats ===\n")
cat("  Total DELs:", nrow(binary_gt), "\n")
cat("  Samples:", length(samples), "\n")
cat("  PC1 variance:", var_pct[1], "%\n")
cat("  PC2 variance:", var_pct[2], "%\n")
cat("  Median DELs/sample:", median(per_sample$n_DELs), "\n")
cat("  Median DEL size:", median(svlen_dist$svlen_bp), "bp\n")

# Write figure legend text
legend_file <- file.path(opt$plot_dir, "figure_legend.txt")
writeLines(c(
  "Figure X. Structural deletion variants in catfish hatchery cohort.",
  "",
  paste0("(A) Genome-wide distribution of deletions across ", length(samples),
         " samples. Heatmap shows DEL count per non-overlapping 1-Mb window. ",
         "Top barplot: cumulative DEL count. Right barplot: per-sample total."),
  paste0("(B) Per-sample DEL counts, sorted by total. ",
         "Teal = unrelated subset (n=", length(unrelated_samples), ")."),
  paste0("(C) PCA based on binary DEL genotype matrix. ",
         "PC1 explains ", var_pct[1], "% and PC2 explains ", var_pct[2], "% of variance."),
  "(D) Pairwise shared DEL heatmap with hierarchical clustering (Ward's method).",
  paste0("(E) DEL size distribution (log scale). Median size: ",
         round(median(svlen_dist$svlen_bp)), " bp."),
  "(F) DEL counts per chromosome/linkage group."
), legend_file)

cat("\nAll plots written to:", opt$plot_dir, "\n")
cat("Done.\n")
