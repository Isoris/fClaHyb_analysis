#!/usr/bin/env Rscript

# =============================================================================
# STEP11_overlap_candidate_regions_with_theta_and_het.R
#
# Overlap candidate local-PCA outlier regions (from STEP10) with per-window
# theta/diversity data (from STEP10b or similar).
#
# The theta input should be a POPULATION-LEVEL summary (mean tP across samples
# per window), not per-sample data. For per-sample overlays, use STEP12 directly.
#
# This step produces:
#   - candidate-level summary (mean/median/min/max of theta in overlapping windows)
#   - detailed overlap table (every theta window × candidate intersection)
#
# Usage:
#   Rscript STEP11_overlap_candidate_regions_with_theta_and_het.R \
#     <candidate_regions.tsv.gz> \
#     <theta_windows.tsv.gz> \
#     <outprefix> \
#     [theta_chrom_col=chrom] \
#     [theta_start_col=WinStart] \
#     [theta_end_col=WinStop] \
#     [theta_value_col=tP]
#
# The theta_windows file can be:
#   - Output of STEP10b (population mean per window)
#   - Any TSV with columns: chrom, start, end, and one or more numeric value columns
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(paste(
    "Usage: Rscript STEP11_overlap_candidate_regions_with_theta_and_het.R",
    "<candidate_regions.tsv.gz> <theta_windows.tsv.gz> <outprefix>",
    "[theta_chrom_col=chrom] [theta_start_col=WinStart] [theta_end_col=WinStop]",
    "[theta_value_col=tP]"
  ))
}

cand_file       <- args[1]
thet_file       <- args[2]
outprefix       <- args[3]
theta_chrom_col <- if (length(args) >= 4) args[4] else "chrom"
theta_start_col <- if (length(args) >= 5) args[5] else "WinStart"
theta_end_col   <- if (length(args) >= 6) args[6] else "WinStop"
theta_value_col <- if (length(args) >= 7) args[7] else "tP"

# ── Read inputs ────────────────────────────────────────────────────────────
cand <- fread(cand_file)
thet <- fread(thet_file)

message("[INFO] Candidate regions: ", nrow(cand))
message("[INFO] Theta windows:     ", nrow(thet))

# Validate candidate columns
req_cand <- c("candidate_id", "chrom", "start_bp", "end_bp")
miss_cand <- setdiff(req_cand, names(cand))
if (length(miss_cand) > 0) stop("Candidate file missing required columns: ", paste(miss_cand, collapse = ", "))

# Validate theta columns
req_thet <- c(theta_chrom_col, theta_start_col, theta_end_col)
miss_thet <- setdiff(req_thet, names(thet))
if (length(miss_thet) > 0) stop("Theta file missing required columns: ", paste(miss_thet, collapse = ", "))

# ── Standardize column names ──────────────────────────────────────────────
# Rename theta columns to standard names for foverlaps
if (theta_chrom_col != "chrom")    setnames(thet, theta_chrom_col, "chrom")
if (theta_start_col != "start_bp") setnames(thet, theta_start_col, "start_bp")
if (theta_end_col   != "end_bp")   setnames(thet, theta_end_col,   "end_bp")

cand[, start_bp := as.numeric(start_bp)]
cand[, end_bp   := as.numeric(end_bp)]
thet[, start_bp := as.numeric(start_bp)]
thet[, end_bp   := as.numeric(end_bp)]

# Add theta window ID if not present
if (!("theta_window_id" %in% names(thet))) {
  thet[, theta_window_id := .I]
}

# ── foverlaps ──────────────────────────────────────────────────────────────
# y = candidate regions (must have key set on chrom, start_bp, end_bp)
# x = theta windows (queried against y)
#
# In the output:
#   start_bp, end_bp = from y (candidate regions)
#   i.start_bp, i.end_bp = from x (theta windows)
setkey(cand, chrom, start_bp, end_bp)
setkey(thet, chrom, start_bp, end_bp)

ov <- foverlaps(
  x = thet,
  y = cand,
  by.x = c("chrom", "start_bp", "end_bp"),
  by.y = c("chrom", "start_bp", "end_bp"),
  type = "any",
  nomatch = 0L
)

if (nrow(ov) == 0) {
  warning("No overlaps found between candidate regions and theta windows")
}

message("[INFO] Overlap rows: ", nrow(ov))

# Compute overlap extent
# In foverlaps output: start_bp/end_bp = candidate (y), i.start_bp/i.end_bp = theta (x)
ov[, overlap_start := pmax(start_bp, i.start_bp)]
ov[, overlap_end   := pmin(end_bp,   i.end_bp)]
ov[, overlap_bp    := overlap_end - overlap_start]
ov <- ov[overlap_bp > 0]

# ── Summarize per candidate ───────────────────────────────────────────────
# Find all numeric columns from theta (excluding coordinate columns)
skip_cols <- c("chrom", "start_bp", "end_bp", "theta_window_id",
               "i.start_bp", "i.end_bp", "overlap_start", "overlap_end", "overlap_bp",
               # Also skip candidate columns that got merged in
               "candidate_id", "mds_axis", "n_windows",
               "first_global_window_id", "last_global_window_id")
numeric_cols <- setdiff(names(thet), c("chrom", "start_bp", "end_bp", "theta_window_id"))
numeric_cols <- numeric_cols[sapply(thet[, ..numeric_cols], is.numeric)]

summary_list <- vector("list", nrow(cand))

for (i in seq_len(nrow(cand))) {
  rid <- cand$candidate_id[i]
  sub <- ov[candidate_id == rid]

  if (nrow(sub) == 0) {
    base <- copy(cand[i])
    base[, n_theta_windows := 0L]
    base[, overlap_bp_total := 0]
    summary_list[[i]] <- base
    next
  }

  base <- copy(cand[i])
  base[, n_theta_windows := uniqueN(sub$theta_window_id)]
  base[, overlap_bp_total := sum(sub$overlap_bp, na.rm = TRUE)]

  for (cc in numeric_cols) {
    if (cc %in% names(sub)) {
      vals <- sub[[cc]]
      vals <- vals[is.finite(vals)]
      if (length(vals) > 0) {
        base[[paste0(cc, "_mean")]]   <- mean(vals)
        base[[paste0(cc, "_median")]] <- median(vals)
        base[[paste0(cc, "_min")]]    <- min(vals)
        base[[paste0(cc, "_max")]]    <- max(vals)
      }
    }
  }

  summary_list[[i]] <- base
}

cand_summary <- rbindlist(summary_list, fill = TRUE)

# ── Build detail table ────────────────────────────────────────────────────
# Rename i.start_bp / i.end_bp to theta_start_bp / theta_end_bp for clarity
detail_dt <- copy(ov)
if ("i.start_bp" %in% names(detail_dt)) setnames(detail_dt, "i.start_bp", "theta_start_bp")
if ("i.end_bp"   %in% names(detail_dt)) setnames(detail_dt, "i.end_bp",   "theta_end_bp")

# Also rename candidate start/end for clarity
setnames(detail_dt, "start_bp", "candidate_start_bp")
setnames(detail_dt, "end_bp",   "candidate_end_bp")

# ── Write outputs ──────────────────────────────────────────────────────────
summary_out <- paste0(outprefix, ".candidate_theta_summary.tsv.gz")
detail_out  <- paste0(outprefix, ".candidate_theta_detail.tsv.gz")
rds_out     <- paste0(outprefix, ".candidate_theta_overlap.rds")

fwrite(cand_summary, summary_out, sep = "\t")
fwrite(detail_dt, detail_out, sep = "\t")
saveRDS(
  list(
    candidate_regions  = cand,
    theta_windows      = thet,
    overlap_detail     = detail_dt,
    candidate_summary  = cand_summary
  ),
  rds_out
)

message("[DONE] Wrote:")
message("  ", summary_out)
message("  ", detail_out)
message("  ", rds_out)
