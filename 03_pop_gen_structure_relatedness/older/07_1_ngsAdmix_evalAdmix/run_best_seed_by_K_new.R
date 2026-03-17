#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# 04_best_seed_by_K.R
#
# Purpose:
#   1) Summarize NGSadmix + evalAdmix metrics across K and seeds
#   2) Select best seed per K
#   3) Copy/link the selected best files to standardized *_best names
#   4) Write stable cluster color palettes per K
#   5) Write canonical sample order table
#   6) Write per-sample dominant ancestry/color table from best qopt per K
#
# Outputs:
#   - all_seed_metrics_by_K.tsv
#   - best_seed_by_K.tsv
#   - best_seed_copied_files.tsv
#   - cluster_palette_by_K.tsv
#   - sample_order_reference.tsv
#   - sample_main_ancestry_by_K.tsv
#
# Notes:
#   - sample order is assumed to follow the canonical sample file / manifest
#   - qopt rows must match canonical sample count exactly
# =============================================================================

# =============================================================================
# PATHS / SETTINGS (HPC)
# =============================================================================

base_dir <- "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/05_ngsadmix_global"

run_dir  <- file.path(base_dir, "runs_thin500")
eval_dir <- file.path(base_dir, "evaladmix_thin500")

# Canonical sample order input:
# Prefer sample_file if available; otherwise manifest_file with sample_id column.
sample_file   <- "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"
manifest_file <- "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/pa_roary/pav/bam_manifest.tsv"

# If TRUE, overwrite existing *_best files
overwrite_existing <- TRUE

# If TRUE, try symlink first; if FALSE, copy directly
use_symlink <- FALSE

# K range and seeds to evaluate
K_values    <- 2:12
seed_values <- 1:3

# Prefix style:
# thin500_K08_seed1.log
# thin500_K08_seed1.corres.txt
thin_label <- "500"

# Stable palette name
palette_name <- "catfish_ngsadmix_v1"

# Output files
out_all              <- file.path(base_dir, "all_seed_metrics_by_K.tsv")
out_best             <- file.path(base_dir, "best_seed_by_K.tsv")
out_copied           <- file.path(base_dir, "best_seed_copied_files.tsv")
out_palette          <- file.path(base_dir, "cluster_palette_by_K.tsv")
out_sample_order     <- file.path(base_dir, "sample_order_reference.tsv")
out_sample_ancestry  <- file.path(base_dir, "sample_main_ancestry_by_K.tsv")

# =============================================================================
# HELPERS
# =============================================================================

extract_ll_from_log <- function(logfile) {
  if (!file.exists(logfile)) return(NA_real_)

  x <- readLines(logfile, warn = FALSE)
  hit <- grep("(like|LL|logl|log-likelihood)", x, ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) return(NA_real_)

  last <- tail(hit, 1)
  nums <- regmatches(last, gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", last))[[1]]
  if (length(nums) == 0) return(NA_real_)

  as.numeric(tail(nums, 1))
}

read_evaladmix_metrics <- function(corfile) {
  if (!file.exists(corfile)) {
    return(list(
      mean_abs_resid = NA_real_,
      max_abs_resid  = NA_real_,
      mean_raw_resid = NA_real_
    ))
  }

  M <- tryCatch(
    as.matrix(read.table(corfile, header = FALSE, sep = "", check.names = FALSE)),
    error = function(e) NULL
  )

  if (is.null(M)) {
    warning("Could not read evalAdmix matrix: ", corfile)
    return(list(
      mean_abs_resid = NA_real_,
      max_abs_resid  = NA_real_,
      mean_raw_resid = NA_real_
    ))
  }

  if (nrow(M) != ncol(M)) {
    warning("Skipping non-square evalAdmix matrix: ", corfile)
    return(list(
      mean_abs_resid = NA_real_,
      max_abs_resid  = NA_real_,
      mean_raw_resid = NA_real_
    ))
  }

  diag(M) <- NA_real_

  list(
    mean_abs_resid = mean(abs(M), na.rm = TRUE),
    max_abs_resid  = max(abs(M), na.rm = TRUE),
    mean_raw_resid = mean(M, na.rm = TRUE)
  )
}

build_prefix <- function(K, seed, thin_label = "500") {
  sprintf("thin%s_K%02d_seed%d", thin_label, K, seed)
}

copy_or_link <- function(from, to, overwrite = TRUE, use_symlink = FALSE) {
  if (!file.exists(from)) {
    warning("Missing source file: ", from)
    return(FALSE)
  }

  if (file.exists(to)) {
    if (!overwrite) {
      message("Exists, skipping: ", to)
      return(TRUE)
    }
    ok_rm <- file.remove(to)
    if (!ok_rm) {
      warning("Could not remove existing file: ", to)
      return(FALSE)
    }
  }

  if (use_symlink) {
    ok <- tryCatch(file.symlink(from, to), error = function(e) FALSE)
    if (isTRUE(ok)) return(TRUE)
  }

  ok <- file.copy(from, to, overwrite = overwrite)
  if (!ok) warning("Failed to copy: ", from, " -> ", to)
  ok
}

read_samples_canonical <- function(sample_file = NULL, manifest_file = NULL) {
  if (!is.null(sample_file) && file.exists(sample_file)) {
    x <- fread(sample_file, header = FALSE)
    samples <- as.character(x[[1]])
    if (length(samples) == 0) stop("sample_file is empty: ", sample_file)

    return(data.table(
      sample_index = seq_along(samples),
      sample = samples,
      source_type = "sample_file",
      source_file = sample_file
    ))
  }

  if (!is.null(manifest_file) && file.exists(manifest_file)) {
    x <- fread(manifest_file)
    if (!"sample_id" %in% names(x)) {
      stop("manifest_file must contain sample_id column: ", manifest_file)
    }
    samples <- as.character(x$sample_id)
    if (length(samples) == 0) stop("manifest_file contains zero samples: ", manifest_file)

    return(data.table(
      sample_index = seq_along(samples),
      sample = samples,
      source_type = "manifest_file",
      source_file = manifest_file
    ))
  }

  stop("No valid sample_file or manifest_file found")
}

read_qopt_matrix <- function(qopt_file, K) {
  if (!file.exists(qopt_file)) {
    stop("Missing qopt file: ", qopt_file)
  }

  lines <- readLines(qopt_file, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  if (length(lines) == 0) {
    stop("qopt file is empty: ", qopt_file)
  }

  split_rows <- strsplit(lines, "[[:space:]]+")
  lens <- vapply(split_rows, length, integer(1))

  if (length(unique(lens)) != 1) {
    stop("Inconsistent qopt row widths in: ", qopt_file,
         " (widths=", paste(unique(lens), collapse = ","), ")")
  }

  qmat <- do.call(rbind, lapply(split_rows, as.numeric))
  mode(qmat) <- "numeric"

  if (ncol(qmat) != K) {
    stop(sprintf("qopt columns (%d) != K (%d) in %s", ncol(qmat), K, qopt_file))
  }

  qmat
}

make_cluster_palette <- function(K, palette_name = "catfish_ngsadmix_v1") {
  base_cols <- c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
    "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
    "#9C755F", "#BAB0AC", "#86BCB6", "#D37295"
  )

  cols <- if (K <= length(base_cols)) {
    base_cols[seq_len(K)]
  } else {
    grDevices::colorRampPalette(base_cols)(K)
  }

  data.table(
    K = K,
    cluster_index = seq_len(K),
    cluster_label = paste0("K", seq_len(K)),
    color_hex = cols,
    palette_name = palette_name
  )
}

# =============================================================================
# READ CANONICAL SAMPLE ORDER
# =============================================================================

sample_dt <- read_samples_canonical(
  sample_file   = sample_file,
  manifest_file = manifest_file
)

n_samples_canonical <- nrow(sample_dt)
sample_source_type  <- unique(sample_dt$source_type)
sample_source_file  <- unique(sample_dt$source_file)

if (length(sample_source_type) != 1L || length(sample_source_file) != 1L) {
  stop("Unexpected multiple sample sources detected")
}

fwrite(sample_dt, out_sample_order, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote sample-order reference:\n  ", out_sample_order, "\n", sep = "")
cat("Canonical sample count: ", n_samples_canonical, "\n", sep = "")

# =============================================================================
# STEP 1: SUMMARIZE ALL RUNS
# =============================================================================

all_rows <- rbindlist(lapply(K_values, function(K) {

  rbindlist(lapply(seed_values, function(seed) {

    prefix <- build_prefix(K, seed, thin_label = thin_label)

    log_file <- file.path(run_dir,  paste0(prefix, ".log"))
    cor_file <- file.path(eval_dir, paste0(prefix, ".corres.txt"))
    q_file   <- file.path(run_dir,  paste0(prefix, ".qopt"))
    f_file   <- file.path(run_dir,  paste0(prefix, ".fopt.gz"))

    ll <- extract_ll_from_log(log_file)
    em <- read_evaladmix_metrics(cor_file)

    qopt_nrows <- if (file.exists(q_file)) {
      length(readLines(q_file, warn = FALSE))
    } else {
      NA_integer_
    }

    qopt_ncols <- if (file.exists(q_file)) {
      tryCatch({
        qmat_tmp <- read_qopt_matrix(q_file, K)
        ncol(qmat_tmp)
      }, error = function(e) NA_integer_)
    } else {
      NA_integer_
    }

    data.table(
      K = K,
      seed = seed,
      seed_tag = paste0("seed", seed),
      prefix = prefix,

      qopt_exists = file.exists(q_file),
      fopt_exists = file.exists(f_file),
      log_exists  = file.exists(log_file),
      cor_exists  = file.exists(cor_file),

      qopt_nrows = qopt_nrows,
      qopt_ncols = qopt_ncols,
      expected_n_samples = n_samples_canonical,
      qopt_rows_match_expected = !is.na(qopt_nrows) & (qopt_nrows == n_samples_canonical),

      loglik = ll,
      mean_abs_resid = em$mean_abs_resid,
      max_abs_resid  = em$max_abs_resid,
      mean_raw_resid = em$mean_raw_resid,

      q_file   = q_file,
      f_file   = f_file,
      log_file = log_file,
      cor_file = cor_file
    )
  }))
}), use.names = TRUE)

all_rows[, usable_for_ranking := qopt_exists & fopt_exists & log_exists & qopt_rows_match_expected]

rank_dt <- copy(all_rows)

rank_dt[is.na(loglik),          loglik_rank   := -Inf]
rank_dt[!is.na(loglik),         loglik_rank   := loglik]

rank_dt[is.na(mean_abs_resid),  mean_abs_rank := Inf]
rank_dt[!is.na(mean_abs_resid), mean_abs_rank := mean_abs_resid]

rank_dt[is.na(max_abs_resid),   max_abs_rank  := Inf]
rank_dt[!is.na(max_abs_resid),  max_abs_rank  := max_abs_resid]

setorder(rank_dt, K, -loglik_rank, mean_abs_rank, max_abs_rank, seed)

best_by_K <- rank_dt[usable_for_ranking == TRUE, .SD[1], by = K]

if (nrow(best_by_K) == 0) {
  stop("No usable best seeds found. Check qopt/fopt/log existence and sample-count match.")
}

best_by_K[, selection_rule := "highest loglik, then lowest mean_abs_resid, then lowest max_abs_resid"]
best_by_K[, n_samples_expected := n_samples_canonical]
best_by_K[, sample_source_type := sample_source_type]
best_by_K[, sample_source_file := sample_source_file]
best_by_K[, palette_name := palette_name]

fwrite(rank_dt[, .(
  K, seed, seed_tag, prefix,
  qopt_exists, fopt_exists, log_exists, cor_exists,
  qopt_nrows, qopt_ncols, expected_n_samples, qopt_rows_match_expected,
  usable_for_ranking,
  loglik, mean_abs_resid, max_abs_resid, mean_raw_resid,
  q_file, f_file, log_file, cor_file
)], out_all, sep = "\t", quote = FALSE, na = "NA")

fwrite(best_by_K[, .(
  K,
  seed,
  seed_tag,
  prefix,
  loglik,
  mean_abs_resid,
  max_abs_resid,
  mean_raw_resid,
  qopt_nrows,
  qopt_ncols,
  n_samples_expected,
  sample_source_type,
  sample_source_file,
  palette_name,
  selection_rule,
  q_file,
  f_file,
  log_file,
  cor_file
)], out_best, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote full metrics table:\n  ", out_all, "\n", sep = "")
cat("Wrote best-seed summary:\n  ", out_best, "\n", sep = "")

cat("\nBest seed by K:\n")
print(best_by_K[, .(
  K, seed_tag, loglik, mean_abs_resid, max_abs_resid, qopt_nrows
)])

# =============================================================================
# STEP 2: COPY / LINK BEST FILES
# =============================================================================

summary_rows <- rbindlist(lapply(seq_len(nrow(best_by_K)), function(i) {

  K <- best_by_K$K[i]
  prefix <- best_by_K$prefix[i]

  src_qopt <- file.path(run_dir,  paste0(prefix, ".qopt"))
  src_fopt <- file.path(run_dir,  paste0(prefix, ".fopt.gz"))
  src_log  <- file.path(run_dir,  paste0(prefix, ".log"))
  src_cor  <- file.path(eval_dir, paste0(prefix, ".corres.txt"))

  dst_base_run  <- sprintf("thin%s_K%02d_best", thin_label, K)
  dst_base_eval <- sprintf("thin%s_K%02d_best", thin_label, K)

  dst_qopt <- file.path(run_dir,  paste0(dst_base_run,  ".qopt"))
  dst_fopt <- file.path(run_dir,  paste0(dst_base_run,  ".fopt.gz"))
  dst_log  <- file.path(run_dir,  paste0(dst_base_run,  ".log"))
  dst_cor  <- file.path(eval_dir, paste0(dst_base_eval, ".corres.txt"))

  ok_qopt <- copy_or_link(src_qopt, dst_qopt, overwrite = overwrite_existing, use_symlink = use_symlink)
  ok_fopt <- copy_or_link(src_fopt, dst_fopt, overwrite = overwrite_existing, use_symlink = use_symlink)
  ok_log  <- copy_or_link(src_log,  dst_log,  overwrite = overwrite_existing, use_symlink = use_symlink)
  ok_cor  <- copy_or_link(src_cor,  dst_cor,  overwrite = overwrite_existing, use_symlink = use_symlink)

  data.table(
    K = K,
    selected_prefix = prefix,
    best_qopt = dst_qopt,
    best_fopt = dst_fopt,
    best_log  = dst_log,
    best_cor  = dst_cor,
    qopt_ok = ok_qopt,
    fopt_ok = ok_fopt,
    log_ok  = ok_log,
    cor_ok  = ok_cor
  )
}), use.names = TRUE)

fwrite(summary_rows, out_copied, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote copied-file summary:\n  ", out_copied, "\n", sep = "")
cat("\nCreated best files for K values:\n")
print(summary_rows[, .(K, selected_prefix, qopt_ok, fopt_ok, log_ok, cor_ok)])

# =============================================================================
# STEP 3: WRITE CLUSTER PALETTE TABLE
# =============================================================================

palette_dt <- rbindlist(
  lapply(sort(unique(best_by_K$K)), function(K) make_cluster_palette(K, palette_name = palette_name)),
  use.names = TRUE
)

fwrite(palette_dt, out_palette, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote palette table:\n  ", out_palette, "\n", sep = "")

# =============================================================================
# STEP 4: WRITE SAMPLE MAIN ANCESTRY TABLE FROM BEST QOPT PER K
# =============================================================================

sample_anc_list <- lapply(seq_len(nrow(best_by_K)), function(i) {
  K_use <- best_by_K$K[i]
  qf    <- best_by_K$q_file[i]

  if (!file.exists(qf)) {
    warning("Skipping missing qopt for K=", K_use, ": ", qf)
    return(NULL)
  }

  qmat <- read_qopt_matrix(qf, K_use)

  if (nrow(qmat) != n_samples_canonical) {
    stop(sprintf(
      "K=%d: qopt rows (%d) do not match canonical sample count (%d)",
      K_use, nrow(qmat), n_samples_canonical
    ))
  }

  qdt <- copy(sample_dt[, .(sample_index, sample)])

  for (j in seq_len(K_use)) {
    qdt[[paste0("Q", j)]] <- qmat[, j]
  }

  qcols <- paste0("Q", seq_len(K_use))

  qdt[, qmax := apply(.SD, 1, max), .SDcols = qcols]
  qdt[, cluster_index := max.col(as.matrix(.SD), ties.method = "first"), .SDcols = qcols]
  qdt[, cluster_label := paste0("K", cluster_index)]
  qdt[, K := K_use]
  qdt[, q_file := qf]
  qdt[, n_samples_expected := n_samples_canonical]
  qdt[, sample_source_type := sample_source_type]
  qdt[, sample_source_file := sample_source_file]
  qdt[, palette_name := palette_name]

  pal_sub <- palette_dt[K == K_use, .(cluster_index, cluster_label, color_hex)]

  qdt <- merge(
    qdt,
    pal_sub,
    by = c("cluster_index", "cluster_label"),
    all.x = TRUE,
    sort = FALSE
  )

  setcolorder(qdt, c(
    "K",
    "sample_index",
    "sample",
    "q_file",
    "cluster_index",
    "cluster_label",
    "qmax",
    "color_hex",
    qcols,
    "n_samples_expected",
    "sample_source_type",
    "sample_source_file",
    "palette_name"
  ))

  qdt
})

sample_anc_dt <- rbindlist(sample_anc_list, use.names = TRUE, fill = TRUE)

fwrite(sample_anc_dt, out_sample_ancestry, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote sample ancestry/color table:\n  ", out_sample_ancestry, "\n", sep = "")

# =============================================================================
# FINAL MESSAGE
# =============================================================================

cat("\n============================================================\n")
cat("NGSadmix best-seed metadata build complete\n")
cat("Base dir: ", base_dir, "\n", sep = "")
cat("Canonical sample count: ", n_samples_canonical, "\n", sep = "")
cat("Outputs:\n")
cat("  - ", out_all, "\n", sep = "")
cat("  - ", out_best, "\n", sep = "")
cat("  - ", out_copied, "\n", sep = "")
cat("  - ", out_palette, "\n", sep = "")
cat("  - ", out_sample_order, "\n", sep = "")
cat("  - ", out_sample_ancestry, "\n", sep = "")
cat("============================================================\n")
