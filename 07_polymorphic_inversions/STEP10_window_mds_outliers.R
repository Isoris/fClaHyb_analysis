#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript STEP10_window_mds_outliers.R <step09_rds_dir> <outprefix> [npc=2] [mds_dims=20] [z_thresh=3] [gap_bp=1000000] [min_windows=3]")
}

rds_dir <- args[1]
outprefix <- args[2]
npc <- if (length(args) >= 3) as.integer(args[3]) else 2L
mds_dims <- if (length(args) >= 4) as.integer(args[4]) else 20L
z_thresh <- if (length(args) >= 5) as.numeric(args[5]) else 3
gap_bp <- if (length(args) >= 6) as.numeric(args[6]) else 1000000
min_windows <- if (length(args) >= 7) as.integer(args[7]) else 3L

rds_files <- list.files(rds_dir, pattern = "\\.window_pca\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop("No STEP09 .window_pca.rds files found")

dist_sq_from_pcs <- function(values1, vectors1, values2, vectors2) {
  Xt <- crossprod(vectors2, vectors1)
  sum(values1^2) + sum(values2^2) - 2 * sum(values1 * colSums(Xt * (values2 * Xt)))
}

pc_dist_from_step09 <- function(pca_df, sample_names, npc, normalize = "L1") {
  values <- as.matrix(pca_df[, paste0("lam_", seq_len(npc)), drop = FALSE])
  vec_cols <- unlist(lapply(seq_len(npc), function(pc) paste0("PC_", pc, "_", sample_names)))
  vectors <- as.matrix(pca_df[, vec_cols, drop = FALSE])

  if (normalize == "L1") {
    rs <- rowSums(abs(values))
    rs[rs == 0] <- NA
    values <- values / rs
  }

  n <- nrow(values)
  out <- matrix(NA_real_, n, n)
  emat <- function(u) matrix(u, ncol = npc)

  for (i in seq_len(n)) {
    vi <- values[i, ]
    Ui <- emat(vectors[i, ])
    for (j in i:n) {
      vj <- values[j, ]
      Uj <- emat(vectors[j, ])
      if (anyNA(vi) || anyNA(vj) || anyNA(Ui) || anyNA(Uj)) {
        d2 <- NA_real_
      } else {
        d2 <- dist_sq_from_pcs(vi, Ui, vj, Uj)
        if (!is.finite(d2) || d2 < 0) d2 <- 0
      }
      out[i, j] <- sqrt(d2)
      out[j, i] <- out[i, j]
    }
  }
  out
}

cluster_outliers_bp <- function(dt, flag_col, gap_bp, min_windows) {
  idx <- which(dt[[flag_col]])
  if (length(idx) == 0) return(NULL)
  clusters <- list()
  cur <- idx[1]

  if (length(idx) > 1) {
    for (ii in idx[-1]) {
      prev_end <- dt$end_bp[max(cur)]
      this_start <- dt$start_bp[ii]
      if ((this_start - prev_end) <= gap_bp) {
        cur <- c(cur, ii)
      } else {
        if (length(cur) >= min_windows) clusters[[length(clusters) + 1]] <- cur
        cur <- ii
      }
    }
  }

  if (length(cur) >= min_windows) clusters[[length(clusters) + 1]] <- cur
  clusters
}

meta_all <- list()
pca_all <- list()
sample_names_ref <- NULL

for (f in rds_files) {
  obj <- readRDS(f)
  if (is.null(sample_names_ref)) {
    sample_names_ref <- obj$sample_names
  } else if (!identical(sample_names_ref, obj$sample_names)) {
    stop("Sample names differ across STEP09 files")
  }
  meta_all[[length(meta_all) + 1]] <- as.data.table(obj$window_meta)
  pca_all[[length(pca_all) + 1]] <- as.data.table(obj$pca)
}

meta_dt <- rbindlist(meta_all, fill = TRUE)
pca_dt <- rbindlist(pca_all, fill = TRUE)
meta_dt[, global_window_id := .I]
pca_dt[, global_window_id := .I]

setkey(meta_dt, global_window_id)
setkey(pca_dt, global_window_id)
dt <- merge(meta_dt, pca_dt, by = c("global_window_id", "window_id"))

dmat <- pc_dist_from_step09(as.data.frame(dt), sample_names_ref, npc = npc, normalize = "L1")
keep <- which(apply(dmat, 1, function(x) all(is.finite(x))))
if (length(keep) < 3) stop("Too few windows with finite distances")

dt_keep <- dt[keep]
dmat_keep <- dmat[keep, keep, drop = FALSE]

k_mds <- min(mds_dims, nrow(dmat_keep) - 1L)
mds <- cmdscale(as.dist(dmat_keep), k = k_mds, eig = TRUE)

mds_dt <- as.data.table(mds$points)
setnames(mds_dt, paste0("MDS", seq_len(ncol(mds_dt))))
mds_dt[, global_window_id := dt_keep$global_window_id]

for (i in seq_len(ncol(mds$points))) {
  coln <- paste0("MDS", i)
  zcol <- paste0("MDS", i, "_z")
  vv <- mds_dt[[coln]]
  mds_dt[[zcol]] <- (vv - mean(vv, na.rm = TRUE)) / sd(vv, na.rm = TRUE)
  mds_dt[[paste0("MDS", i, "_outlier")]] <- abs(mds_dt[[zcol]]) >= z_thresh
}

setkey(mds_dt, global_window_id)
setkey(dt_keep, global_window_id)
out_dt <- merge(dt_keep, mds_dt, by = "global_window_id")

cand_list <- list()
cand_id <- 1L

for (i in seq_len(k_mds)) {
  flag_col <- paste0("MDS", i, "_outlier")
  for (chr in unique(out_dt$chrom)) {
    sub <- out_dt[chrom == chr][order(start_bp)]
    clusters <- cluster_outliers_bp(sub, flag_col, gap_bp = gap_bp, min_windows = min_windows)
    if (!is.null(clusters)) {
      for (cl in clusters) {
        xx <- sub[cl]
        cand_list[[length(cand_list) + 1]] <- data.table(
          candidate_id = cand_id,
          mds_axis = i,
          chrom = chr,
          start_bp = min(xx$start_bp),
          end_bp = max(xx$end_bp),
          n_windows = nrow(xx),
          first_global_window_id = min(xx$global_window_id),
          last_global_window_id = max(xx$global_window_id)
        )
        cand_id <- cand_id + 1L
      }
    }
  }
}

cand_dt <- if (length(cand_list) > 0) rbindlist(cand_list) else data.table(
  candidate_id = integer(),
  mds_axis = integer(),
  chrom = character(),
  start_bp = numeric(),
  end_bp = numeric(),
  n_windows = integer(),
  first_global_window_id = integer(),
  last_global_window_id = integer()
)

mds_out <- paste0(outprefix, ".window_mds.tsv.gz")
cand_out <- paste0(outprefix, ".candidate_regions.tsv.gz")
dist_rds_out <- paste0(outprefix, ".pcdist.rds")
mds_rds_out <- paste0(outprefix, ".mds.rds")

fwrite(out_dt, mds_out, sep = "\t")
fwrite(cand_dt, cand_out, sep = "\t")
saveRDS(dmat_keep, dist_rds_out)
saveRDS(list(dt = out_dt, candidate_regions = cand_dt, cmdscale = mds), mds_rds_out)

message("[DONE] Wrote:")
message("  ", mds_out)
message("  ", cand_out)
message("  ", dist_rds_out)
message("  ", mds_rds_out)
