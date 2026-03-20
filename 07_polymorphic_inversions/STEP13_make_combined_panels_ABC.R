#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    paste(
      "Usage:",
      "Rscript STEP13_make_combined_panels_ABC.R",
      "<step12_prefix>",
      "<candidate_id>",
      "<outprefix>"
    )
  )
}

step12_prefix <- args[1]
candidate_id <- args[2]
outprefix <- args[3]

pcs_file <- paste0(step12_prefix, ".candidate_", candidate_id, ".regional_pca_samples.tsv.gz")
plotc_file <- paste0(step12_prefix, ".candidate_", candidate_id, ".plotC_summary.tsv.gz")
rep_file <- paste0(step12_prefix, ".candidate_", candidate_id, ".representatives.tsv.gz")

for (f in c(pcs_file, plotc_file, rep_file)) {
  if (!file.exists(f)) stop("Missing input file: ", f)
}

pcs <- fread(pcs_file)
plotc_dt <- fread(plotc_file)
reps <- fread(rep_file)

req_pcs <- c("sample", "PC1", "PC2", "ordered_group", "regional_het")
miss_pcs <- setdiff(req_pcs, names(pcs))
if (length(miss_pcs) > 0) {
  stop("regional_pca_samples file missing required columns: ", paste(miss_pcs, collapse = ", "))
}

req_plotc <- c("group_label", "ordered_group", "mid_bp", "het_mean")
miss_plotc <- setdiff(req_plotc, names(plotc_dt))
if (length(miss_plotc) > 0) {
  stop("plotC_summary file missing required columns: ", paste(miss_plotc, collapse = ", "))
}

# Try to recover region info from pcs if present in attributes? not available.
# So we infer shaded region only if region columns were carried through reps or pcs.
region_start <- if ("candidate_start_bp" %in% names(reps)) unique(reps$candidate_start_bp)[1] else NA_real_
region_end   <- if ("candidate_end_bp" %in% names(reps)) unique(reps$candidate_end_bp)[1] else NA_real_

# Panel labels like the example
label_map <- data.table(
  ordered_group = sort(unique(pcs$ordered_group)),
  group_tag = c("3a", "3b", "5a")[seq_along(sort(unique(pcs$ordered_group)))]
)

pcs <- merge(pcs, label_map, by = "ordered_group", all.x = TRUE)
reps <- merge(reps, label_map, by = "ordered_group", all.x = TRUE)
plotc_dt <- merge(plotc_dt, label_map, by = "ordered_group", all.x = TRUE)

# Panel A
p_a <- ggplot(pcs, aes(x = PC1, y = PC2)) +
  geom_point(shape = 1, alpha = 0.6, size = 2) +
  theme_bw(base_size = 12) +
  labs(
    title = "(a)",
    x = "PC1",
    y = "PC2"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0)
  )

# Panel B
p_b <- ggplot(pcs, aes(x = PC1, y = PC2, color = regional_het, label = ordered_group)) +
  geom_text(size = 3) +
  geom_point(
    data = reps,
    aes(x = PC1, y = PC2),
    inherit.aes = FALSE,
    shape = 1,
    size = 6,
    stroke = 1.2,
    color = "black"
  ) +
  geom_text(
    data = reps,
    aes(x = PC1, y = PC2, label = group_tag),
    inherit.aes = FALSE,
    nudge_y = 1.5,
    color = "black",
    size = 4
  ) +
  theme_bw(base_size = 12) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "(b)",
    x = "PC1",
    y = "PC2",
    color = "Heterozygosity"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0)
  )

# Panel C
p_c <- ggplot(plotc_dt, aes(x = mid_bp / 1e6, y = het_mean, color = group_tag)) +
  {
    if (is.finite(region_start) && is.finite(region_end)) {
      annotate(
        "rect",
        xmin = region_start / 1e6,
        xmax = region_end / 1e6,
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.15,
        fill = "grey60"
      )
    }
  } +
  geom_line(linewidth = 1) +
  geom_text(
    data = plotc_dt[, .SD[which.max(mid_bp)], by = group_tag],
    aes(label = group_tag),
    inherit.aes = FALSE,
    hjust = -0.1,
    size = 4,
    color = "black"
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "(c)",
    x = "Chromosome position (Mb)",
    y = "Heterozygosity (%)",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "none"
  )

# Layout: A over B over C
combined <- p_a / p_b / p_c + plot_layout(heights = c(1, 1, 0.8))

pdf_out <- paste0(outprefix, ".candidate_", candidate_id, ".combined_ABC.pdf")
png_out <- paste0(outprefix, ".candidate_", candidate_id, ".combined_ABC.png")

ggsave(pdf_out, combined, width = 7.5, height = 12)
ggsave(png_out, combined, width = 7.5, height = 12, dpi = 300)

message("[DONE] Wrote:")
message("  ", pdf_out)
message("  ", png_out)
