#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# =========================================================
# USER SETTINGS
# =========================================================
input_dir <- "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/03-mosdepth_LG12/windows_50bp"
out_dir   <- file.path(input_dir, "analysis_R")

# coverage thresholds interpreted as fraction of each 50 bp window
window_bp <- 50
low_cut   <- 0.10   # <= 10% of bases in window reach threshold
high_cut  <- 0.80   # >= 80% of bases in window reach threshold

# Candidate P/A-like windows:
# at least this many samples low, and at least this many samples high
min_n_low8  <- 5
min_n_high8 <- 20

# Merge adjacent candidate windows into larger blocks if gap <= this many bp
merge_gap_bp <- 50

# number of samples to show in facet depth plots
n_facet_samples <- 18

# if TRUE, also read .regions.bed.gz and make per-sample depth facet plots
make_depth_facets <- TRUE

# =========================================================
# SETUP
# =========================================================
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[INFO] input_dir = ", input_dir)
message("[INFO] out_dir   = ", out_dir)

thr_files <- Sys.glob(file.path(input_dir, "*.thresholds.bed.gz"))
if (length(thr_files) == 0) {
  stop("[ERROR] No *.thresholds.bed.gz files found in: ", input_dir)
}

reg_files <- Sys.glob(file.path(input_dir, "*.regions.bed.gz"))

# =========================================================
# HELPERS
# =========================================================
extract_sample_name <- function(path, suffix_regex) {
  x <- basename(path)
  sub(suffix_regex, "", x)
}

read_threshold_file <- function(f, window_bp = 50) {
  samp <- extract_sample_name(f, "\\.C_gar_LG12\\.w50\\.thresholds\\.bed\\.gz$")
  dt <- fread(cmd = paste("zcat", shQuote(f)), sep = "\t", header = TRUE)

  # clean header names
  setnames(dt, gsub("^#", "", names(dt)))

  # expected columns: chrom,start,end,region,1X,4X,8X,12X,20X
  needed <- c("chrom", "start", "end", "region", "1X", "4X", "8X", "12X", "20X")
  missing_cols <- setdiff(needed, names(dt))
  if (length(missing_cols) > 0) {
    stop("[ERROR] Missing columns in ", f, ": ", paste(missing_cols, collapse = ", "))
  }

  dt[, sample := samp]
  dt[, mid := (start + end) / 2]

  # proportions of bases in each 50 bp window reaching threshold
  dt[, prop1  := `1X`  / window_bp]
  dt[, prop4  := `4X`  / window_bp]
  dt[, prop8  := `8X`  / window_bp]
  dt[, prop12 := `12X` / window_bp]
  dt[, prop20 := `20X` / window_bp]

  dt[]
}

read_region_file <- function(f) {
  samp <- extract_sample_name(f, "\\.C_gar_LG12\\.w50\\.regions\\.bed\\.gz$")
  dt <- fread(cmd = paste("zcat", shQuote(f)), sep = "\t", header = FALSE)

  if (ncol(dt) < 4) {
    stop("[ERROR] Unexpected .regions.bed.gz format in ", f)
  }

  # columns are usually: chrom start end depth
  setnames(dt, c("chrom", "start", "end", "depth")[seq_len(ncol(dt))])
  dt[, sample := samp]
  dt[, mid := (start + end) / 2]
  dt[]
}

merge_candidate_windows <- function(dt, gap_bp = 50) {
  # dt must be sorted and contain: chrom, start, end, n_low8, n_high8, mean_prop8, median_prop8
  if (nrow(dt) == 0) return(copy(dt)[, grp := integer()])

  dt <- copy(dt)
  setorder(dt, chrom, start, end)

  grp <- integer(nrow(dt))
  grp[1] <- 1L

  for (i in 2:nrow(dt)) {
    same_chr <- dt$chrom[i] == dt$chrom[i - 1]
    close_enough <- (dt$start[i] - dt$end[i - 1]) <= gap_bp
    grp[i] <- if (same_chr && close_enough) grp[i - 1] else grp[i - 1] + 1L
  }

  dt[, grp := grp]
  dt
}

safe_write <- function(dt, path) {
  fwrite(dt, path, sep = "\t", quote = FALSE, na = "NA")
  message("[INFO] wrote: ", path)
}

# =========================================================
# READ THRESHOLDS
# =========================================================
message("[INFO] Reading threshold files...")
thr_list <- lapply(thr_files, read_threshold_file, window_bp = window_bp)
thr <- rbindlist(thr_list, use.names = TRUE, fill = TRUE)

n_samples <- uniqueN(thr$sample)
message("[INFO] Samples detected: ", n_samples)
message("[INFO] Total rows: ", nrow(thr))

# expected windows should be same across samples
window_counts <- thr[, .N, by = sample][order(sample)]
safe_write(window_counts, file.path(out_dir, "LG12_per_sample_window_counts.tsv"))

# =========================================================
# PER-WINDOW SUMMARY ACROSS ALL SAMPLES
# =========================================================
message("[INFO] Building per-window summary...")
win <- thr[, .(
  n_samples     = .N,
  n_low1        = sum(prop1  <= low_cut, na.rm = TRUE),
  n_low4        = sum(prop4  <= low_cut, na.rm = TRUE),
  n_low8        = sum(prop8  <= low_cut, na.rm = TRUE),
  n_low12       = sum(prop12 <= low_cut, na.rm = TRUE),
  n_low20       = sum(prop20 <= low_cut, na.rm = TRUE),

  n_high1       = sum(prop1  >= high_cut, na.rm = TRUE),
  n_high4       = sum(prop4  >= high_cut, na.rm = TRUE),
  n_high8       = sum(prop8  >= high_cut, na.rm = TRUE),
  n_high12      = sum(prop12 >= high_cut, na.rm = TRUE),
  n_high20      = sum(prop20 >= high_cut, na.rm = TRUE),

  mean_prop1    = mean(prop1,  na.rm = TRUE),
  mean_prop4    = mean(prop4,  na.rm = TRUE),
  mean_prop8    = mean(prop8,  na.rm = TRUE),
  mean_prop12   = mean(prop12, na.rm = TRUE),
  mean_prop20   = mean(prop20, na.rm = TRUE),

  median_prop1  = median(prop1,  na.rm = TRUE),
  median_prop4  = median(prop4,  na.rm = TRUE),
  median_prop8  = median(prop8,  na.rm = TRUE),
  median_prop12 = median(prop12, na.rm = TRUE),
  median_prop20 = median(prop20, na.rm = TRUE)
), by = .(chrom, start, end, mid)]

setorder(win, chrom, start, end)
safe_write(win, file.path(out_dir, "LG12_window_summary_across_samples.tsv"))

# =========================================================
# CANDIDATE P/A-LIKE WINDOWS
# =========================================================
message("[INFO] Calling candidate P/A-like windows...")
cand <- win[
  n_low8  >= min_n_low8 &
  n_high8 >= min_n_high8
]

safe_write(cand, file.path(out_dir, "LG12_candidate_PAV_windows.tsv"))

# =========================================================
# MERGE CANDIDATE WINDOWS INTO BLOCKS
# =========================================================
message("[INFO] Merging candidate windows into blocks...")
cand_grp <- merge_candidate_windows(cand, gap_bp = merge_gap_bp)

if (nrow(cand_grp) > 0) {
  blocks <- cand_grp[, .(
    chrom         = first(chrom),
    start         = min(start),
    end           = max(end),
    block_bp      = max(end) - min(start),
    n_windows     = .N,
    max_n_low8    = max(n_low8, na.rm = TRUE),
    mean_n_low8   = mean(n_low8, na.rm = TRUE),
    max_n_high8   = max(n_high8, na.rm = TRUE),
    mean_n_high8  = mean(n_high8, na.rm = TRUE),
    mean_prop8    = mean(mean_prop8, na.rm = TRUE),
    median_prop8  = median(median_prop8, na.rm = TRUE)
  ), by = grp][order(chrom, start, end)]
} else {
  blocks <- data.table(
    grp = integer(),
    chrom = character(),
    start = integer(),
    end = integer(),
    block_bp = integer(),
    n_windows = integer(),
    max_n_low8 = integer(),
    mean_n_low8 = numeric(),
    max_n_high8 = integer(),
    mean_n_high8 = numeric(),
    mean_prop8 = numeric(),
    median_prop8 = numeric()
  )
}

safe_write(blocks, file.path(out_dir, "LG12_candidate_PAV_blocks.tsv"))

# =========================================================
# PER-SAMPLE BURDEN
# =========================================================
message("[INFO] Computing per-sample low-coverage burden...")
sample_burden <- thr[, .(
  n_windows      = .N,
  frac_low1      = mean(prop1  <= low_cut, na.rm = TRUE),
  frac_low4      = mean(prop4  <= low_cut, na.rm = TRUE),
  frac_low8      = mean(prop8  <= low_cut, na.rm = TRUE),
  frac_low12     = mean(prop12 <= low_cut, na.rm = TRUE),
  frac_low20     = mean(prop20 <= low_cut, na.rm = TRUE),
  mean_prop8     = mean(prop8, na.rm = TRUE),
  median_prop8   = median(prop8, na.rm = TRUE)
), by = sample][order(-frac_low8, sample)]

safe_write(sample_burden, file.path(out_dir, "LG12_per_sample_low_coverage_burden.tsv"))

# =========================================================
# OPTIONAL: READ REGIONS FOR DEPTH FACETS
# =========================================================
reg <- NULL
if (make_depth_facets) {
  if (length(reg_files) == 0) {
    warning("[WARN] No *.regions.bed.gz files found; skipping facet depth plots.")
    make_depth_facets <- FALSE
  } else {
    message("[INFO] Reading region depth files...")
    reg_list <- lapply(reg_files, read_region_file)
    reg <- rbindlist(reg_list, use.names = TRUE, fill = TRUE)
    safe_write(
      reg[, .(chrom, start, end, mid, sample, depth)],
      file.path(out_dir, "LG12_regions_depth_long.tsv")
    )
  }
}

# =========================================================
# PLOT 1: SUMMARY TRACK
# =========================================================
message("[INFO] Plotting summary track...")
p_track1 <- ggplot(win, aes(x = mid, y = n_low8)) +
  geom_line(linewidth = 0.3) +
  labs(
    x = "LG12 position (bp)",
    y = "Number of samples with low 8X support",
    title = "LG12 recurrent low-coverage windows across samples"
  ) +
  theme_minimal(base_size = 11)

p_track2 <- ggplot(win, aes(x = mid, y = n_high8)) +
  geom_line(linewidth = 0.3) +
  labs(
    x = "LG12 position (bp)",
    y = "Number of samples with high 8X support",
    title = "LG12 well-covered windows across samples"
  ) +
  theme_minimal(base_size = 11)

track_pdf <- file.path(out_dir, "LG12_summary_tracks_nlow8_nhigh8.pdf")
ggsave(track_pdf, p_track1 / p_track2, width = 13, height = 7)
message("[INFO] wrote: ", track_pdf)

# =========================================================
# PLOT 2: HEATMAP OF PROP8
# =========================================================
message("[INFO] Plotting heatmap...")
sample_order <- sample_burden$sample
thr[, sample := factor(sample, levels = sample_order)]

p_heat <- ggplot(thr, aes(x = mid, y = sample, fill = prop8)) +
  geom_raster() +
  scale_fill_viridis_c(
    name = "Prop bases\n>=8X",
    limits = c(0, 1)
  ) +
  labs(
    x = "LG12 position (bp)",
    y = "Sample",
    title = "LG12 sample-by-window coverage heatmap (8X threshold)"
  ) +
  theme_minimal(base_size = 9)

heat_pdf <- file.path(out_dir, "LG12_prop8_heatmap.pdf")
ggsave(heat_pdf, p_heat, width = 14, height = max(8, 0.08 * n_samples))
message("[INFO] wrote: ", heat_pdf)

# =========================================================
# PLOT 3: FACET DEPTH PLOTS (SELECTED SAMPLES)
# =========================================================
if (make_depth_facets && !is.null(reg)) {
  message("[INFO] Plotting per-sample facet depth profiles...")
  top_samples <- head(sample_burden$sample, n_facet_samples)
  reg_sub <- reg[sample %in% top_samples]
  reg_sub[, sample := factor(sample, levels = top_samples)]

  p_facet <- ggplot(reg_sub, aes(x = mid, y = depth)) +
    geom_line(linewidth = 0.2) +
    facet_wrap(~ sample, ncol = 3, scales = "free_y") +
    labs(
      x = "LG12 position (bp)",
      y = "Median depth per 50 bp window",
      title = "LG12 depth profiles for selected samples"
    ) +
    theme_minimal(base_size = 9)

  facet_pdf <- file.path(out_dir, "LG12_depth_facets_selected_samples.pdf")
  ggsave(facet_pdf, p_facet, width = 14, height = 10)
  message("[INFO] wrote: ", facet_pdf)
}

# =========================================================
# PLOT 4: HISTOGRAM OF PER-SAMPLE LOW BURDEN
# =========================================================
message("[INFO] Plotting sample burden histogram...")
p_burden <- ggplot(sample_burden, aes(x = frac_low8)) +
  geom_histogram(bins = 40) +
  labs(
    x = "Fraction of LG12 windows with low 8X support",
    y = "Number of samples",
    title = "Distribution of per-sample LG12 low-coverage burden"
  ) +
  theme_minimal(base_size = 11)

burden_pdf <- file.path(out_dir, "LG12_sample_low8_burden_histogram.pdf")
ggsave(burden_pdf, p_burden, width = 8, height = 5)
message("[INFO] wrote: ", burden_pdf)

# =========================================================
# TEXT SUMMARY
# =========================================================
message("[INFO] Writing summary report...")
summary_txt <- file.path(out_dir, "LG12_analysis_summary.txt")
sink(summary_txt)

cat("=== LG12 mosdepth threshold analysis summary ===\n")
cat("Input directory: ", input_dir, "\n", sep = "")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Samples analyzed: ", n_samples, "\n", sep = "")
cat("Window size (bp): ", window_bp, "\n", sep = "")
cat("Low cutoff:  prop <= ", low_cut, "\n", sep = "")
cat("High cutoff: prop >= ", high_cut, "\n", sep = "")
cat("Candidate P/A thresholds: n_low8 >= ", min_n_low8,
    " and n_high8 >= ", min_n_high8, "\n", sep = "")
cat("Merge gap (bp): ", merge_gap_bp, "\n\n", sep = "")

cat("--- Per-window summary ---\n")
cat("Total windows: ", nrow(win), "\n", sep = "")
cat("Candidate windows: ", nrow(cand), "\n", sep = "")
cat("Candidate merged blocks: ", nrow(blocks), "\n\n", sep = "")

if (nrow(blocks) > 0) {
  cat("--- Top 10 candidate blocks by size ---\n")
  print(head(blocks[order(-block_bp)], 10))
  cat("\n")
}

cat("--- Per-sample low burden (top 20 by frac_low8) ---\n")
print(head(sample_burden, 20))
cat("\n")

cat("--- Files written ---\n")
cat("LG12_window_summary_across_samples.tsv\n")
cat("LG12_candidate_PAV_windows.tsv\n")
cat("LG12_candidate_PAV_blocks.tsv\n")
cat("LG12_per_sample_low_coverage_burden.tsv\n")
cat("LG12_summary_tracks_nlow8_nhigh8.pdf\n")
cat("LG12_prop8_heatmap.pdf\n")
if (make_depth_facets) cat("LG12_depth_facets_selected_samples.pdf\n")
cat("LG12_sample_low8_burden_histogram.pdf\n")

sink()
message("[INFO] wrote: ", summary_txt)

message("[DONE] All analyses completed.")
