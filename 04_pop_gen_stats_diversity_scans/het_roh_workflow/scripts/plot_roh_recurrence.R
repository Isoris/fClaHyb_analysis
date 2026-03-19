#!/usr/bin/env Rscript
# =============================================================================
# plot_roh_recurrence.R — ROH recurrence across windows per chromosome
# =============================================================================
# Usage:
#   Rscript plot_roh_recurrence.R <recurrence_dir> <out_dir>
#   Expects: recurrence_roh_windows_{100kb,250kb,500kb,1Mb}.tsv
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: plot_roh_recurrence.R <table_dir> <out_dir>")

table_dir <- args[1]
out_dir   <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))

labels <- c("100kb", "250kb", "500kb", "1Mb")

for (lbl in labels) {
  fname <- file.path(table_dir, paste0("recurrence_roh_windows_", lbl, ".tsv"))
  if (!file.exists(fname)) {
    cat("Skipping", lbl, "(file not found)\n")
    next
  }

  d <- fread(fname)
  # Use pct_with_roh_all column
  pct_col <- grep("^pct_with_roh_all$", names(d), value = TRUE)
  if (length(pct_col) == 0) pct_col <- grep("^pct_with_roh", names(d), value = TRUE)[1]
  if (is.na(pct_col)) next

  d[, pct := as.numeric(get(pct_col))]
  d[, mid := (start + end) / 2]

  # Natural chromosome order
  chroms <- unique(d$chrom)
  chrom_nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chroms)))
  chroms <- chroms[order(ifelse(is.na(chrom_nums), 999, chrom_nums), chroms)]
  d[, chrom := factor(chrom, levels = chroms)]

  p <- ggplot(d, aes(x = mid / 1e6, y = pct)) +
    geom_area(fill = "steelblue", alpha = 0.4) +
    geom_line(color = "steelblue", linewidth = 0.3) +
    facet_wrap(~ chrom, scales = "free_x", ncol = 4) +
    labs(
      title = paste0("ROH recurrence across samples (", lbl, " windows)"),
      x = "Position (Mb)", y = "% samples with ROH"
    ) +
    theme_pub +
    theme(axis.text.x = element_text(size = 6))

  ggsave(file.path(out_dir, paste0("roh_recurrence_", lbl, ".png")),
         p, width = 14, height = 10, dpi = 400)
  ggsave(file.path(out_dir, paste0("roh_recurrence_", lbl, ".pdf")),
         p, width = 14, height = 10)
}

cat("Recurrence plots written to:", out_dir, "\n")
