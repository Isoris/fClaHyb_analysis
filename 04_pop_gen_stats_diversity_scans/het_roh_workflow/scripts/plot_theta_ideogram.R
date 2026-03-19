#!/usr/bin/env Rscript
# =============================================================================
# plot_theta_ideogram.R — Local theta/diversity ideogram-like tracks
# =============================================================================
# Produces genome-wide ideogram views of local theta (diversity proxy).
# Automatically detects all thetaStat scales present in the directory:
#   sample.win5000.step1000.pestPG
#   sample.win10000.step2000.pestPG
#   sample.win50000.step10000.pestPG
# etc.
#
# NOTE: These are LOCAL DIVERSITY TRACKS, not literal per-site Hobs.
#
# Usage:
#   Rscript plot_theta_ideogram.R \
#     <theta_dir> <chrom_sizes.tsv> <out_dir> [sample_list.txt] [max_samples=10]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: plot_theta_ideogram.R <theta_dir> <chrom_sizes.tsv> <out_dir> [sample_list] [max_samples]")
}

theta_dir   <- args[1]
sizes_file  <- args[2]
out_dir     <- args[3]
sample_file <- if (length(args) >= 4 && file.exists(args[4])) args[4] else NULL
max_samples <- if (length(args) >= 5) as.integer(args[5]) else 10

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95")
  )

# ── Chromosome sizes ─────────────────────────────────────────────────────
sizes <- fread(sizes_file, header = FALSE, col.names = c("chrom", "len"))
chrom_order <- sizes$chrom
sizes[, cum_start := cumsum(c(0, len[-.N]))]

# ── Find theta window files ──────────────────────────────────────────────
theta_files <- list.files(theta_dir, pattern = "\\.pestPG$", full.names = TRUE)
if (length(theta_files) == 0) {
  cat("No .pestPG files found in", theta_dir, "\n")
  quit(save = "no")
}

# Optional sample restriction
keep_samples <- NULL
if (!is.null(sample_file)) {
  keep_samples <- fread(sample_file, header = FALSE)$V1
  keep_samples <- unique(as.character(keep_samples))
}

# ── Helpers ──────────────────────────────────────────────────────────────
get_sample <- function(f) {
  bn <- basename(f)
  sub("\\.win[0-9]+\\.step[0-9]+\\.pestPG$", "", bn)
}

get_scale <- function(f) {
  bn <- basename(f)
  m <- regexec("\\.win([0-9]+)\\.step([0-9]+)\\.pestPG$", bn)
  reg <- regmatches(bn, m)[[1]]
  if (length(reg) != 3) return(NULL)
  list(
    win = as.integer(reg[2]),
    step = as.integer(reg[3]),
    label = paste0("win", reg[2], "_step", reg[3])
  )
}

safe_read_theta <- function(f) {
  # thetaStat writes a header line starting with "#"
  # fread(comment.char = "#") would drop the header entirely, so we read manually
  lines <- readLines(f, warn = FALSE)
  if (length(lines) < 2) return(NULL)

  header_line <- lines[1]
  data_lines  <- lines[-1]
  if (length(data_lines) == 0) return(NULL)

  # remove leading # from header
  header_line <- sub("^#", "", header_line)

  tf <- tempfile(fileext = ".tsv")
  writeLines(c(header_line, data_lines), tf)

  out <- tryCatch(
    fread(tf, sep = "\t", header = TRUE, data.table = TRUE),
    error = function(e) NULL
  )
  unlink(tf)
  out
}

detect_theta_cols <- function(t) {
  cn <- names(t)

  chr_col <- cn[cn %in% c("Chr", "chr", "Chromo", "chrom")][1]
  center_col <- cn[cn %in% c("WinCenter", "midPos", "wincenter")][1]
  tp_col <- cn[cn %in% c("tP", "thetaPi", "ThetaPi", "pi")][1]
  nsites_col <- cn[cn %in% c("nSites", "nsite", "Nsites")][1]

  # fallback for standard thetaStat output
  # expected names often include:
  # "(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)"
  # "Chr" "WinCenter" "tW" "tP" ... "nSites"
  if (is.na(chr_col) && "Chr" %in% cn) chr_col <- "Chr"
  if (is.na(center_col) && "WinCenter" %in% cn) center_col <- "WinCenter"
  if (is.na(tp_col) && "tP" %in% cn) tp_col <- "tP"
  if (is.na(nsites_col) && "nSites" %in% cn) nsites_col <- "nSites"

  list(
    chr_col = chr_col,
    center_col = center_col,
    tp_col = tp_col,
    nsites_col = nsites_col
  )
}

nat_sort_chr <- function(x) {
  nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", x)))
  x[order(ifelse(is.na(nums), 999999, nums), x)]
}

# ── Build file table with sample and scale ──────────────────────────────
file_dt <- data.table(file = theta_files)
file_dt[, sample := vapply(file, get_sample, character(1))]

scale_info <- lapply(theta_files, get_scale)
file_dt[, win := vapply(scale_info, function(z) if (is.null(z)) NA_integer_ else z$win, integer(1))]
file_dt[, step := vapply(scale_info, function(z) if (is.null(z)) NA_integer_ else z$step, integer(1))]
file_dt[, scale_label := vapply(scale_info, function(z) if (is.null(z)) NA_character_ else z$label, character(1))]

file_dt <- file_dt[!is.na(win) & !is.na(step)]

if (!is.null(keep_samples)) {
  file_dt <- file_dt[sample %in% keep_samples]
}

if (nrow(file_dt) == 0) {
  cat("No matching theta files after filtering.\n")
  quit(save = "no")
}

# limit samples globally
all_samples <- unique(file_dt$sample)
if (length(all_samples) > max_samples) {
  cat("Plotting first", max_samples, "of", length(all_samples), "samples\n")
  keep <- all_samples[1:max_samples]
  file_dt <- file_dt[sample %in% keep]
}

# ── Loop over scales ─────────────────────────────────────────────────────
scale_dt <- unique(file_dt[, .(win, step, scale_label)])
setorder(scale_dt, win, step)

for (k in seq_len(nrow(scale_dt))) {
  WIN <- scale_dt$win[k]
  STEP <- scale_dt$step[k]
  SCALE_LABEL <- scale_dt$scale_label[k]

  cat("Processing scale:", SCALE_LABEL, "\n")

  sub_files <- file_dt[win == WIN & step == STEP]

  all_theta <- list()

  for (i in seq_len(nrow(sub_files))) {
    samp <- sub_files$sample[i]
    f    <- sub_files$file[i]

    t <- safe_read_theta(f)
    if (is.null(t) || nrow(t) == 0) next

    cols <- detect_theta_cols(t)
    chr_col    <- cols$chr_col
    center_col <- cols$center_col
    tp_col     <- cols$tp_col
    nsites_col <- cols$nsites_col

    if (any(is.na(c(chr_col, center_col, tp_col)))) {
      cat("Skipping file due to missing expected columns:", basename(f), "\n")
      next
    }

    keep_cols <- c(chr_col, center_col, tp_col)
    if (!is.na(nsites_col)) keep_cols <- c(keep_cols, nsites_col)

    sub_t <- t[, ..keep_cols]
    setnames(sub_t, old = c(chr_col, center_col, tp_col), new = c("chrom", "center", "tP"))
    if (!is.na(nsites_col)) {
      setnames(sub_t, old = nsites_col, new = "nSites")
    }

    sub_t[, sample := samp]
    sub_t[, win := WIN]
    sub_t[, step := STEP]
    sub_t[, scale_label := SCALE_LABEL]

    sub_t[, center := as.numeric(center)]
    sub_t[, tP := as.numeric(tP)]

    if ("nSites" %in% names(sub_t)) {
      sub_t[, nSites := as.numeric(nSites)]
      sub_t[nSites > 0, theta_per_site := tP / nSites]
      sub_t[nSites <= 0 | is.na(nSites), theta_per_site := NA_real_]
    } else {
      sub_t[, nSites := NA_real_]
      sub_t[, theta_per_site := tP]
    }

    sub_t <- sub_t[!is.na(chrom) & !is.na(center) & !is.na(theta_per_site)]
    if (nrow(sub_t) == 0) next

    all_theta[[paste0(samp, "_", SCALE_LABEL)]] <- sub_t
  }

  if (length(all_theta) == 0) {
    cat("No theta data loaded for scale:", SCALE_LABEL, "\n")
    next
  }

  dt <- rbindlist(all_theta, fill = TRUE)
  dt[, chrom := as.character(chrom)]

  # reorder chromosomes naturally if needed
  chr_ord <- nat_sort_chr(unique(dt$chrom))
  # keep genome size order preference if present
  chr_ord2 <- chrom_order[chrom_order %in% chr_ord]
  chr_ord  <- c(chr_ord2, setdiff(chr_ord, chr_ord2))

  dt <- merge(dt, sizes[, .(chrom, cum_start)], by = "chrom", all.x = TRUE)
  dt <- dt[!is.na(cum_start)]
  dt[, gw_pos := cum_start + center]
  dt[, chrom := factor(chrom, levels = chr_ord)]

  # ── Per-sample ideogram ───────────────────────────────────────────────
  for (samp in unique(dt$sample)) {
    sub <- dt[sample == samp]
    if (nrow(sub) < 10) next

    p <- ggplot(sub, aes(x = gw_pos / 1e6, y = theta_per_site)) +
      geom_point(size = 0.3, alpha = 0.4, color = "steelblue") +
      geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "red", linewidth = 0.5) +
      labs(
        title = paste0(samp, " — Local theta diversity proxy (", SCALE_LABEL, ")"),
        subtitle = "NOT literal Hobs",
        x = "Genomic position (Mb, genome-wide)",
        y = "Theta/site (diversity proxy)"
      ) +
      theme_pub

    ggsave(
      file.path(out_dir, paste0(samp, "_theta_genomewide_", SCALE_LABEL, ".png")),
      p, width = 14, height = 4, dpi = 300
    )
    ggsave(
      file.path(out_dir, paste0(samp, "_theta_genomewide_", SCALE_LABEL, ".pdf")),
      p, width = 14, height = 4
    )

    p2 <- ggplot(sub, aes(x = center / 1e6, y = theta_per_site)) +
      geom_point(size = 0.2, alpha = 0.3, color = "steelblue") +
      facet_wrap(~ chrom, scales = "free_x", ncol = 4) +
      labs(
        title = paste0(samp, " — Per-chromosome local diversity proxy (", SCALE_LABEL, ")"),
        subtitle = "NOT literal Hobs",
        x = "Position (Mb)",
        y = "Theta/site"
      ) +
      theme_pub +
      theme(axis.text.x = element_text(size = 6))

    ggsave(
      file.path(out_dir, paste0(samp, "_theta_per_chr_", SCALE_LABEL, ".png")),
      p2, width = 14, height = 10, dpi = 300
    )
    ggsave(
      file.path(out_dir, paste0(samp, "_theta_per_chr_", SCALE_LABEL, ".pdf")),
      p2, width = 14, height = 10
    )
  }

  # ── Multi-sample overlay (mean across samples per window) ────────────
  mean_dt <- dt[, .(
    mean_theta = mean(theta_per_site, na.rm = TRUE),
    sd_theta = sd(theta_per_site, na.rm = TRUE),
    n_samples = .N
  ), by = .(chrom, center)]

  mean_dt <- merge(mean_dt, sizes[, .(chrom, cum_start)], by = "chrom", all.x = TRUE)
  mean_dt <- mean_dt[!is.na(cum_start)]
  mean_dt[, gw_pos := cum_start + center]
  mean_dt[, chrom := factor(chrom, levels = chr_ord)]

  p_mean <- ggplot(mean_dt, aes(x = gw_pos / 1e6, y = mean_theta)) +
    geom_point(size = 0.3, alpha = 0.4, color = "darkgreen") +
    geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "red", linewidth = 0.5) +
    labs(
      title = paste0("Mean local theta diversity proxy across all samples (", SCALE_LABEL, ")"),
      subtitle = "NOT literal per-site Hobs",
      x = "Genomic position (Mb)",
      y = "Mean theta/site"
    ) +
    theme_pub

  ggsave(
    file.path(out_dir, paste0("mean_theta_genomewide_", SCALE_LABEL, ".png")),
    p_mean, width = 14, height = 4, dpi = 300
  )
  ggsave(
    file.path(out_dir, paste0("mean_theta_genomewide_", SCALE_LABEL, ".pdf")),
    p_mean, width = 14, height = 4
  )

  # Optional: faceted mean per chromosome
  p_mean_chr <- ggplot(mean_dt, aes(x = center / 1e6, y = mean_theta)) +
    geom_point(size = 0.25, alpha = 0.4, color = "darkgreen") +
    geom_smooth(method = "loess", span = 0.05, se = FALSE, color = "red", linewidth = 0.5) +
    facet_wrap(~ chrom, scales = "free_x", ncol = 4) +
    labs(
      title = paste0("Mean per-chromosome theta diversity proxy (", SCALE_LABEL, ")"),
      subtitle = "NOT literal per-site Hobs",
      x = "Position (Mb)",
      y = "Mean theta/site"
    ) +
    theme_pub +
    theme(axis.text.x = element_text(size = 6))

  ggsave(
    file.path(out_dir, paste0("mean_theta_per_chr_", SCALE_LABEL, ".png")),
    p_mean_chr, width = 14, height = 10, dpi = 300
  )
  ggsave(
    file.path(out_dir, paste0("mean_theta_per_chr_", SCALE_LABEL, ".pdf")),
    p_mean_chr, width = 14, height = 10
  )
}

cat("Theta ideogram plots written to:", out_dir, "\n")
