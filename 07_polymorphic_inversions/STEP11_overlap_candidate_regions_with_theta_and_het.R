#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript STEP11_overlap_candidate_regions_with_theta_and_het.R <candidate_regions.tsv.gz> <theta_het_windows.tsv.gz> <outprefix> <region_id_col> [chrom_col=chrom] [start_col=start_bp] [end_col=end_bp]")
}

cand_file <- args[1]
thet_file <- args[2]
outprefix <- args[3]
region_id_col <- args[4]

chrom_col <- if (length(args) >= 5) args[5] else "chrom"
start_col <- if (length(args) >= 6) args[6] else "start_bp"
end_col   <- if (length(args) >= 7) args[7] else "end_bp"

cand <- fread(cand_file)
thet <- fread(thet_file)

req_cand <- c(region_id_col, "chrom", "start_bp", "end_bp")
miss_cand <- setdiff(req_cand, names(cand))
if (length(miss_cand) > 0) stop("Candidate file missing required columns: ", paste(miss_cand, collapse = ", "))

req_thet <- c(chrom_col, start_col, end_col)
miss_thet <- setdiff(req_thet, names(thet))
if (length(miss_thet) > 0) stop("Theta/het file missing required columns: ", paste(miss_thet, collapse = ", "))

setnames(thet, old = chrom_col, new = "chrom")
setnames(thet, old = start_col, new = "start_bp")
setnames(thet, old = end_col, new = "end_bp")

cand[, start_bp := as.numeric(start_bp)]
cand[, end_bp := as.numeric(end_bp)]
thet[, start_bp := as.numeric(start_bp)]
thet[, end_bp := as.numeric(end_bp)]

if (!("theta_window_id" %in% names(thet))) {
  thet[, theta_window_id := .I]
}

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

if (nrow(ov) == 0) warning("No overlaps found")

ov[, overlap_start := pmax(start_bp, i.start_bp)]
ov[, overlap_end   := pmin(end_bp,   i.end_bp)]
ov[, overlap_bp := overlap_end - overlap_start + 1]
ov <- ov[overlap_bp > 0]

numeric_cols <- setdiff(names(thet), c("chrom", "start_bp", "end_bp", "theta_window_id"))
numeric_cols <- numeric_cols[sapply(thet[, ..numeric_cols], is.numeric)]

summary_list <- vector("list", nrow(cand))

for (i in seq_len(nrow(cand))) {
  rid <- cand[[region_id_col]][i]
  sub <- ov[get(region_id_col) == rid]

  if (nrow(sub) == 0) {
    base <- copy(cand[i])
    base[, n_theta_windows := 0L]
    base[, overlap_bp_total := 0]
    summary_list[[i]] <- base
    next
  }

  base <- unique(sub[, ..region_id_col])
  base <- merge(base, cand, by = region_id_col, all.x = TRUE)

  base[, n_theta_windows := uniqueN(sub$theta_window_id)]
  base[, overlap_bp_total := sum(sub$overlap_bp, na.rm = TRUE)]

  for (cc in numeric_cols) {
    vals <- sub[[cc]]
    base[[paste0(cc, "_mean")]] <- mean(vals, na.rm = TRUE)
    base[[paste0(cc, "_median")]] <- median(vals, na.rm = TRUE)
    base[[paste0(cc, "_min")]] <- min(vals, na.rm = TRUE)
    base[[paste0(cc, "_max")]] <- max(vals, na.rm = TRUE)
  }

  summary_list[[i]] <- base
}

cand_summary <- rbindlist(summary_list, fill = TRUE)

detail_cols <- c(
  region_id_col, "chrom",
  "i.start_bp", "i.end_bp",
  "theta_window_id",
  "start_bp", "end_bp",
  "overlap_bp",
  setdiff(names(thet), c("chrom", "start_bp", "end_bp"))
)
detail_cols <- unique(detail_cols[detail_cols %in% names(ov)])
detail_dt <- ov[, ..detail_cols]

if ("i.start_bp" %in% names(detail_dt)) setnames(detail_dt, "i.start_bp", "candidate_start_bp")
if ("i.end_bp" %in% names(detail_dt)) setnames(detail_dt, "i.end_bp", "candidate_end_bp")

summary_out <- paste0(outprefix, ".candidate_theta_het_summary.tsv.gz")
detail_out  <- paste0(outprefix, ".candidate_theta_het_detail.tsv.gz")
rds_out     <- paste0(outprefix, ".candidate_theta_het_overlap.rds")

fwrite(cand_summary, summary_out, sep = "\t")
fwrite(detail_dt, detail_out, sep = "\t")
saveRDS(
  list(
    candidate_regions = cand,
    theta_het = thet,
    overlap_detail = detail_dt,
    candidate_summary = cand_summary
  ),
  rds_out
)

message("[DONE] Wrote:")
message("  ", summary_out)
message("  ", detail_out)
message("  ", rds_out)
