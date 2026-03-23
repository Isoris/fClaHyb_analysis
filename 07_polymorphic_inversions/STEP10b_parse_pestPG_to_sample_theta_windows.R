#!/usr/bin/env Rscript

# =============================================================================
# STEP10b_parse_pestPG_to_sample_theta_windows.R
#
# Bridge script: parse per-sample .pestPG files from the het_roh workflow
# into a single unified table of per-sample windowed tP values suitable
# for STEP11 (overlap) and STEP12 (plot C).
#
# The .pestPG header looks like:
#   #(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop) \
#     Chr WinCenter tW tP tF tH tL Tajima fuf fud fayh zeng nSites
#
# We extract: Chr, WinStart, WinStop, WinCenter, tP, nSites
# and add the sample name from the filename.
#
# Usage:
#   Rscript STEP10b_parse_pestPG_to_sample_theta_windows.R \
#     <pestPG_dir> <samples_file> <outprefix> \
#     [scale_pattern="win50000.step10000"] [chr_filter=""]
#
# Arguments:
#   pestPG_dir      — directory containing *.pestPG files (e.g. het_roh/02_heterozygosity/03_theta/multiscale/)
#   samples_file    — samples list (one per line, or samples.ind with tab-separated fields, col1=sample)
#   outprefix       — output prefix
#   scale_pattern   — filename pattern to select the right scale (default: win50000.step10000)
#   chr_filter      — optional: restrict to this chromosome (e.g. "C_gar_LG12")
#
# Outputs:
#   <outprefix>.sample_tP_windows.tsv.gz  — long-format: sample, chrom, WinStart, WinStop, WinCenter, tP, nSites
#   <outprefix>.sample_tP_windows.rds
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(paste(
    "Usage: Rscript STEP10b_parse_pestPG_to_sample_theta_windows.R",
    "<pestPG_dir> <samples_file> <outprefix>",
    "[scale_pattern=win50000.step10000] [chr_filter=]"
  ))
}

pestPG_dir    <- args[1]
samples_file  <- args[2]
outprefix     <- args[3]
scale_pattern <- if (length(args) >= 4) args[4] else "win50000.step10000"
chr_filter    <- if (length(args) >= 5 && nchar(args[5]) > 0) args[5] else ""

# ── Read sample list ───────────────────────────────────────────────────────
if (!file.exists(samples_file)) stop("Samples file not found: ", samples_file)

# Handle both plain list and .ind format (tab-separated, col1=sample)
samp_raw <- fread(samples_file, header = FALSE, sep = "\t", fill = TRUE)
sample_names <- samp_raw[[1]]
sample_names <- sample_names[nchar(sample_names) > 0]

message("[INFO] Samples loaded: ", length(sample_names))

# ── Find .pestPG files ────────────────────────────────────────────────────
all_files <- list.files(pestPG_dir, pattern = "\\.pestPG$", full.names = TRUE, recursive = TRUE)
if (length(all_files) == 0) stop("No .pestPG files found in: ", pestPG_dir)

# Filter by scale pattern
if (nchar(scale_pattern) > 0) {
  all_files <- all_files[grepl(scale_pattern, basename(all_files), fixed = TRUE)]
}
if (length(all_files) == 0) stop("No .pestPG files matching pattern '", scale_pattern, "' found")

message("[INFO] Found ", length(all_files), " .pestPG files matching '", scale_pattern, "'")

# ── Parse one .pestPG file ────────────────────────────────────────────────
# The header is special: first field is a combined parenthetical mess
# that we need to parse carefully.
# Format: #(idxS,idxE)(fpdS,fpdE)(WinS,WinE)\tChr\tWinCenter\ttW\ttP\t...
parse_one_pestPG <- function(filepath, sample_id) {
  lines <- readLines(filepath, warn = FALSE)
  if (length(lines) < 2) return(NULL)

  # Find header line (starts with #)
  hdr_idx <- grep("^#", lines)
  if (length(hdr_idx) == 0) return(NULL)
  hdr_line <- lines[hdr_idx[1]]

  # Data lines: everything after header, non-empty
  data_lines <- lines[(hdr_idx[1] + 1):length(lines)]
  data_lines <- data_lines[nchar(trimws(data_lines)) > 0]
  if (length(data_lines) == 0) return(NULL)

  # Parse each data line
  # The first field contains parenthetical groups: (a,b)(c,d)(e,f)
  # followed by tab-separated: Chr WinCenter tW tP tF tH tL Tajima fuf fud fayh zeng nSites
  results <- rbindlist(lapply(data_lines, function(line) {
    # Split off the parenthetical prefix from the tab-separated rest
    # The parenthetical part ends after the last ")"
    # Everything after the last ")" is tab-separated fields

    # Strategy: find the position of the last ")" then split
    paren_end <- regexpr("\\)\t", line)
    if (paren_end < 0) return(NULL)

    paren_part <- substr(line, 1, paren_end)
    rest_part  <- substr(line, paren_end + 2, nchar(line))

    # Extract (WinStart,WinStop) from the third parenthetical group
    # Pattern: (num,num)(num,num)(WinStart,WinStop)
    paren_matches <- gregexpr("\\(([^)]+)\\)", paren_part)
    paren_contents <- regmatches(paren_part, paren_matches)[[1]]

    if (length(paren_contents) < 3) return(NULL)

    # Third group = (WinStart,WinStop)
    win_pair <- gsub("[()]", "", paren_contents[3])
    win_parts <- strsplit(win_pair, ",")[[1]]
    if (length(win_parts) != 2) return(NULL)

    WinStart <- as.numeric(win_parts[1])
    WinStop  <- as.numeric(win_parts[2])

    # Parse the rest (tab-separated)
    rest_fields <- strsplit(rest_part, "\t")[[1]]
    # Expected: Chr WinCenter tW tP tF tH tL Tajima fuf fud fayh zeng nSites
    if (length(rest_fields) < 5) return(NULL)

    data.table(
      chrom     = rest_fields[1],
      WinStart  = WinStart,
      WinStop   = WinStop,
      WinCenter = as.numeric(rest_fields[2]),
      tP        = as.numeric(rest_fields[4]),
      nSites    = as.integer(rest_fields[length(rest_fields)])
    )
  }), fill = TRUE)

  if (nrow(results) == 0) return(NULL)

  results[, sample := sample_id]
  results
}

# ── Match files to samples ────────────────────────────────────────────────
# Strategy: for each sample, look for a .pestPG file that contains the sample name
# Common naming: <sample>.win50000.step10000.pestPG or catfish.<sample>.win...
all_results <- list()

for (sid in sample_names) {
  # Find files matching this sample
  matching <- all_files[grepl(sid, basename(all_files), fixed = TRUE)]

  if (length(matching) == 0) {
    message("[WARN] No .pestPG file found for sample: ", sid)
    next
  }

  if (length(matching) > 1) {
    # Prefer exact match over partial
    exact <- matching[grepl(paste0("^", sid, "\\."), basename(matching))]
    if (length(exact) >= 1) matching <- exact[1] else matching <- matching[1]
  } else {
    matching <- matching[1]
  }

  parsed <- parse_one_pestPG(matching, sid)
  if (!is.null(parsed) && nrow(parsed) > 0) {
    all_results[[length(all_results) + 1]] <- parsed
  } else {
    message("[WARN] Could not parse .pestPG for sample: ", sid, " (file: ", matching, ")")
  }
}

if (length(all_results) == 0) stop("No .pestPG data successfully parsed for any sample")

merged <- rbindlist(all_results, fill = TRUE)

# Remove rows with NA tP or non-finite
merged <- merged[is.finite(tP)]

message("[INFO] Parsed ", nrow(merged), " window-sample records across ", uniqueN(merged$sample), " samples")

# ── Optional chromosome filter ─────────────────────────────────────────────
if (nchar(chr_filter) > 0) {
  merged <- merged[chrom == chr_filter]
  message("[INFO] After chr filter '", chr_filter, "': ", nrow(merged), " records")
}

# ── Write outputs ──────────────────────────────────────────────────────────
tsv_out <- paste0(outprefix, ".sample_tP_windows.tsv.gz")
rds_out <- paste0(outprefix, ".sample_tP_windows.rds")

fwrite(merged, tsv_out, sep = "\t")
saveRDS(merged, rds_out)

message("[DONE] Wrote:")
message("  ", tsv_out, "  (", nrow(merged), " rows, ", uniqueN(merged$sample), " samples, ", uniqueN(merged$chrom), " chromosomes)")
message("  ", rds_out)
