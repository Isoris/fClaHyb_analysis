#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(RSpectra)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript STEP09_local_pca_windows_by_chr.R <sites.tsv.gz> <dosage.tsv.gz> <outprefix> [npc=2] [winsize=100]")
}

sites_file <- args[1]
dosage_file <- args[2]
outprefix <- args[3]
npc <- if (length(args) >= 4) as.integer(args[4]) else 2L
winsize <- if (length(args) >= 5) as.integer(args[5]) else 100L

sites <- fread(sites_file)
dos <- fread(dosage_file)

req_sites <- c("marker", "chrom", "pos", "allele1", "allele2")
miss_sites <- setdiff(req_sites, names(sites))
if (length(miss_sites) > 0) stop("Missing required columns in sites file: ", paste(miss_sites, collapse = ", "))
if (!("marker" %in% names(dos))) stop("Dosage file must contain 'marker' column")

if (!identical(sites$marker, dos$marker)) {
  setkeyv(sites, "marker")
  setkeyv(dos, "marker")
  dos <- dos[sites$marker]
  if (!identical(sites$marker, dos$marker)) stop("Sites and dosage markers do not match")
}

sample_names <- setdiff(names(dos), "marker")
chroms <- unique(sites$chrom)
if (length(chroms) != 1) stop("STEP09 expects one chromosome per input file")
chrom <- chroms[1]

X <- as.matrix(dos[, ..sample_names])
storage.mode(X) <- "double"

n_snps <- nrow(X)
n_windows <- floor(n_snps / winsize)
if (n_windows < 1) stop("Not enough SNPs for one full window")

cov_pca_one_window <- function(x, k) {
  rm <- rowMeans(x, na.rm = TRUE)
  xc <- x - rm
  covmat <- suppressWarnings(cov(xc, use = "pairwise.complete.obs"))

  if (anyNA(covmat) || nrow(covmat) < k) {
    return(list(total = NA_real_, eigval = rep(NA_real_, k), eigvec = matrix(NA_real_, nrow = ncol(x), ncol = k)))
  }

  total <- sum(covmat^2)

  if (k == ncol(covmat)) {
    ee <- eigen(covmat, symmetric = TRUE)
    vals <- ee$values[seq_len(k)]
    vecs <- ee$vectors[, seq_len(k), drop = FALSE]
  } else {
    ee <- tryCatch(
      RSpectra::eigs_sym(covmat, k = k, which = "LM"),
      error = function(e) {
        ee2 <- eigen(covmat, symmetric = TRUE)
        list(values = ee2$values[seq_len(k)], vectors = ee2$vectors[, seq_len(k), drop = FALSE])
      }
    )
    vals <- ee$values
    vecs <- ee$vectors
  }

  list(total = total, eigval = vals, eigvec = vecs)
}

meta_list <- vector("list", n_windows)
pca_rows <- vector("list", n_windows)

for (w in seq_len(n_windows)) {
  i1 <- (w - 1L) * winsize + 1L
  i2 <- w * winsize
  xw <- X[i1:i2, , drop = FALSE]
  sw <- sites[i1:i2]
  pca <- cov_pca_one_window(xw, npc)

  meta_list[[w]] <- data.table(
    window_id = w,
    chrom = chrom,
    snp_start_index = i1,
    snp_end_index = i2,
    n_snps = nrow(xw),
    start_bp = sw$pos[1],
    end_bp = sw$pos[nrow(sw)],
    start_marker = sw$marker[1],
    end_marker = sw$marker[nrow(sw)]
  )

  rowvec <- c(total = pca$total, setNames(pca$eigval, paste0("lam_", seq_len(npc))))
  for (pc in seq_len(npc)) {
    rowvec <- c(rowvec, setNames(pca$eigvec[, pc], paste0("PC_", pc, "_", sample_names)))
  }
  pca_rows[[w]] <- rowvec
}

window_meta <- rbindlist(meta_list)
pca_mat <- do.call(rbind, pca_rows)
pca_mat <- as.data.frame(pca_mat)
pca_mat$window_id <- seq_len(n_windows)
pca_mat <- pca_mat[, c("window_id", setdiff(names(pca_mat), "window_id"))]

meta_out <- paste0(outprefix, ".window_meta.tsv.gz")
pca_rds_out <- paste0(outprefix, ".window_pca.rds")
pca_tsv_out <- paste0(outprefix, ".window_pca.tsv.gz")

fwrite(window_meta, meta_out, sep = "\t")
saveRDS(list(chrom = chrom, sample_names = sample_names, npc = npc, winsize = winsize, window_meta = window_meta, pca = pca_mat), pca_rds_out)
fwrite(as.data.table(pca_mat), pca_tsv_out, sep = "\t")

message("[DONE] Wrote:")
message("  ", meta_out)
message("  ", pca_rds_out)
message("  ", pca_tsv_out)#
