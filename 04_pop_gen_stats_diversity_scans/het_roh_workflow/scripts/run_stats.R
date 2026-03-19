#!/usr/bin/env Rscript
# =============================================================================
# run_stats.R — Simple robust statistics for het/ROH summaries
# v2 note:
#   - Keeps original sample-level statistics unchanged
#   - Adds chromosome-level statistics from per_chr_roh_summary.tsv
#   - Optionally adds per-chromosome theta statistics from per_chr_theta_summary_*.tsv
#   - Designed to support breeding-oriented chromosome-aware interpretation
# =============================================================================
# Usage:
#   Rscript run_stats.R <master_summary.tsv> <stats_dir> [ancestry_labels.tsv]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: run_stats.R <master.tsv> <stats_dir> [ancestry]")

master_file <- args[1]
stats_dir   <- args[2]
anc_file    <- if (length(args) >= 3 && file.exists(args[3])) args[3] else NULL

dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
safe_fwrite <- function(x, file) {
  fwrite(x, file, sep = "\t", quote = FALSE, na = "NA")
}

safe_as_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

detect_tables_dir <- function(master_path) {
  normalizePath(dirname(master_path), mustWork = FALSE)
}

get_theta_summary_files <- function(tables_dir) {
  list.files(
    tables_dir,
    pattern = "^per_chr_theta_summary_.*\\.tsv$",
    full.names = TRUE
  )
}

natural_chr_order <- function(chroms) {
  chroms <- unique(chroms)
  nums <- suppressWarnings(as.numeric(gsub("[^0-9]", "", chroms)))
  chroms[order(ifelse(is.na(nums), 999999, nums), chroms)]
}

pairwise_wilcox_dt <- function(dt, value_col, group_col, variable_label, group_prefix = "ancestry") {
  out <- list()
  pw <- pairwise.wilcox.test(dt[[value_col]], dt[[group_col]], p.adjust.method = "BH")
  pw_mat <- pw$p.value
  if (is.null(pw_mat)) return(NULL)

  for (i in seq_len(nrow(pw_mat))) {
    for (j in seq_len(ncol(pw_mat))) {
      if (!is.na(pw_mat[i, j])) {
        out[[length(out) + 1]] <- data.table(
          variable = variable_label,
          grouping = paste0(group_prefix, ":", rownames(pw_mat)[i], "_vs_", colnames(pw_mat)[j]),
          test = "pairwise_Wilcoxon_BH",
          n = nrow(dt),
          n_groups = 2L,
          statistic = NA_real_,
          p_value = NA_real_,
          adj_p_value = pw_mat[i, j],
          direction = ""
        )
      }
    }
  }

  if (length(out) == 0) return(NULL)
  rbindlist(out, fill = TRUE)
}

# ─────────────────────────────────────────────────────────────────────────────
# Load sample-level master table
# ─────────────────────────────────────────────────────────────────────────────
d <- fread(master_file)
d[, het_genomewide := safe_as_num(het_genomewide)]
d[, froh := safe_as_num(froh)]
d[, roh_total_bp := safe_as_num(roh_total_bp)]
d[, longest_roh := safe_as_num(longest_roh)]
if ("n_tracts" %in% names(d)) d[, n_tracts := safe_as_num(n_tracts)]

# ── Merge ancestry if available ──────────────────────────────────────────
if (!is.null(anc_file)) {
  anc <- fread(anc_file)
  if (ncol(anc) >= 2) {
    setnames(anc, 1:2, c("sample", "ancestry"))
    d <- merge(d, anc[, .(sample, ancestry)], by = "sample", all.x = TRUE)
  }
}
if (!"ancestry" %in% names(d)) d[, ancestry := "all"]

# ─────────────────────────────────────────────────────────────────────────────
# 1. Original sample-level stats
# ─────────────────────────────────────────────────────────────────────────────

# ── 1A. Spearman correlations (sample-level) ─────────────────────────────
cor_results <- list()
cor_pairs <- list(
  c("het_genomewide", "froh"),
  c("het_genomewide", "roh_total_bp"),
  c("het_genomewide", "longest_roh")
)

for (pair in cor_pairs) {
  v1 <- pair[1]
  v2 <- pair[2]
  if (!(v1 %in% names(d)) || !(v2 %in% names(d))) next

  x <- d[[v1]]
  y <- d[[v2]]
  ok <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
  if (sum(ok) < 5) next

  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  cor_results[[length(cor_results) + 1]] <- data.table(
    level = "sample",
    variable = paste0(v1, "_vs_", v2),
    context = "all_samples",
    n = sum(ok),
    rho = unname(ct$estimate),
    statistic = unname(ct$statistic),
    p_value = ct$p.value,
    direction = ifelse(unname(ct$estimate) > 0, "positive", "negative")
  )
}
if (length(cor_results) > 0) {
  cor_dt <- rbindlist(cor_results)
  safe_fwrite(cor_dt, file.path(stats_dir, "spearman_correlations.tsv"))
  cat("Wrote spearman_correlations.tsv\n")
}

# ── 1B. Group comparisons (sample-level) ─────────────────────────────────
groups <- unique(d$ancestry[!is.na(d$ancestry)])
if (length(groups) > 1) {
  group_results <- list()
  test_vars <- c("het_genomewide", "froh", "roh_total_bp", "longest_roh", "n_tracts")

  for (v in test_vars) {
    if (!v %in% names(d)) next
    vals <- d[[v]]
    ok <- !is.na(vals) & is.finite(vals) & !is.na(d$ancestry)
    if (sum(ok) < 5) next

    dsub <- d[ok]
    ng <- length(unique(dsub$ancestry))

    if (ng == 2) {
      wt <- wilcox.test(as.formula(paste(v, "~ ancestry")), data = dsub)
      group_results[[length(group_results) + 1]] <- data.table(
        level = "sample",
        variable = v,
        grouping = "ancestry",
        test = "Wilcoxon_rank_sum",
        n = nrow(dsub),
        n_groups = ng,
        statistic = unname(wt$statistic),
        p_value = wt$p.value,
        adj_p_value = NA_real_,
        direction = ""
      )
    } else if (ng > 2) {
      kt <- kruskal.test(as.formula(paste(v, "~ ancestry")), data = dsub)
      group_results[[length(group_results) + 1]] <- data.table(
        level = "sample",
        variable = v,
        grouping = "ancestry",
        test = "Kruskal_Wallis",
        n = nrow(dsub),
        n_groups = ng,
        statistic = unname(kt$statistic),
        p_value = kt$p.value,
        adj_p_value = NA_real_,
        direction = ""
      )

      pw_dt <- pairwise_wilcox_dt(dsub, v, "ancestry", v, "ancestry")
      if (!is.null(pw_dt)) {
        pw_dt[, level := "sample"]
        group_results[[length(group_results) + 1]] <- pw_dt
      }
    }
  }

  if (length(group_results) > 0) {
    group_dt <- rbindlist(group_results, fill = TRUE)
    safe_fwrite(group_dt, file.path(stats_dir, "group_comparisons.tsv"))
    cat("Wrote group_comparisons.tsv\n")
  }
}

# ── 1C. Descriptive summary (sample-level) ───────────────────────────────
desc_vars <- c("het_genomewide", "froh", "roh_total_bp", "longest_roh", "n_tracts")
desc_list <- list()
for (v in desc_vars) {
  if (!v %in% names(d)) next
  x <- safe_as_num(d[[v]])
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) == 0) next
  desc_list[[length(desc_list) + 1]] <- data.table(
    level = "sample",
    variable = v,
    n = length(x),
    mean = mean(x),
    median = median(x),
    sd = sd(x),
    min = min(x),
    max = max(x),
    q25 = as.numeric(quantile(x, 0.25)),
    q75 = as.numeric(quantile(x, 0.75))
  )
}
if (length(desc_list) > 0) {
  desc_dt <- rbindlist(desc_list)
  safe_fwrite(desc_dt, file.path(stats_dir, "descriptive_summary.tsv"))
  cat("Wrote descriptive_summary.tsv\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 2. New chromosome-level ROH stats
# ─────────────────────────────────────────────────────────────────────────────
tables_dir <- detect_tables_dir(master_file)
per_chr_file <- file.path(tables_dir, "per_chr_roh_summary.tsv")

if (file.exists(per_chr_file)) {
  chr_d <- fread(per_chr_file)
  chr_d[, roh_bp := safe_as_num(roh_bp)]
  chr_d[, n_tracts := safe_as_num(n_tracts)]
  chr_d[, longest_tract := safe_as_num(longest_tract)]
  chr_d[, callable_bp_chr := safe_as_num(callable_bp_chr)]
  chr_d[, froh_chr := safe_as_num(froh_chr)]

  # ancestry merge
  if (!is.null(anc_file)) {
    anc <- fread(anc_file)
    if (ncol(anc) >= 2) {
      setnames(anc, 1:2, c("sample", "ancestry"))
      chr_d <- merge(chr_d, anc[, .(sample, ancestry)], by = "sample", all.x = TRUE)
    }
  }
  if (!"ancestry" %in% names(chr_d)) chr_d[, ancestry := "all"]
  chr_d[, chrom := factor(chrom, levels = natural_chr_order(chrom))]

  # ── 2A. Chromosome-level descriptive summary ───────────────────────────
  chr_desc <- chr_d[
    ,
    .(
      n_samples = .N,
      n_with_roh = sum(roh_bp > 0, na.rm = TRUE),
      frac_with_roh = sum(roh_bp > 0, na.rm = TRUE) / .N,
      mean_froh_chr = mean(froh_chr, na.rm = TRUE),
      median_froh_chr = median(froh_chr, na.rm = TRUE),
      sd_froh_chr = sd(froh_chr, na.rm = TRUE),
      min_froh_chr = min(froh_chr, na.rm = TRUE),
      max_froh_chr = max(froh_chr, na.rm = TRUE),
      mean_roh_bp_chr = mean(roh_bp, na.rm = TRUE),
      median_roh_bp_chr = median(roh_bp, na.rm = TRUE),
      mean_n_tracts_chr = mean(n_tracts, na.rm = TRUE),
      mean_longest_tract_chr = mean(longest_tract, na.rm = TRUE)
    ),
    by = chrom
  ]
  safe_fwrite(chr_desc, file.path(stats_dir, "per_chr_descriptive_summary.tsv"))
  cat("Wrote per_chr_descriptive_summary.tsv\n")

  # ── 2B. Chromosome-level correlations ──────────────────────────────────
  chr_cor_list <- list()
  chr_cor_pairs <- list(
    c("froh_chr", "roh_bp"),
    c("froh_chr", "n_tracts"),
    c("froh_chr", "longest_tract"),
    c("roh_bp", "n_tracts")
  )

  for (chr in unique(chr_d$chrom)) {
    sub <- chr_d[chrom == chr]
    for (pair in chr_cor_pairs) {
      v1 <- pair[1]
      v2 <- pair[2]
      if (!(v1 %in% names(sub)) || !(v2 %in% names(sub))) next
      x <- sub[[v1]]
      y <- sub[[v2]]
      ok <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
      if (sum(ok) < 5) next

      ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
      chr_cor_list[[length(chr_cor_list) + 1]] <- data.table(
        level = "chromosome",
        chrom = as.character(chr),
        variable = paste0(v1, "_vs_", v2),
        n = sum(ok),
        rho = unname(ct$estimate),
        statistic = unname(ct$statistic),
        p_value = ct$p.value,
        direction = ifelse(unname(ct$estimate) > 0, "positive", "negative")
      )
    }
  }

  if (length(chr_cor_list) > 0) {
    chr_cor_dt <- rbindlist(chr_cor_list, fill = TRUE)
    safe_fwrite(chr_cor_dt, file.path(stats_dir, "per_chr_correlations.tsv"))
    cat("Wrote per_chr_correlations.tsv\n")
  }

  # ── 2C. Chromosome-level group comparisons ─────────────────────────────
  chr_groups <- unique(chr_d$ancestry[!is.na(chr_d$ancestry)])
  if (length(chr_groups) > 1) {
    chr_group_list <- list()
    chr_test_vars <- c("froh_chr", "roh_bp", "n_tracts", "longest_tract")

    for (chr in unique(chr_d$chrom)) {
      sub <- chr_d[chrom == chr]
      ng <- length(unique(sub$ancestry[!is.na(sub$ancestry)]))
      if (ng < 2) next

      for (v in chr_test_vars) {
        if (!v %in% names(sub)) next
        vals <- sub[[v]]
        ok <- !is.na(vals) & is.finite(vals) & !is.na(sub$ancestry)
        if (sum(ok) < 5) next
        dsub <- sub[ok]
        ng2 <- length(unique(dsub$ancestry))

        if (ng2 == 2) {
          wt <- wilcox.test(as.formula(paste(v, "~ ancestry")), data = dsub)
          chr_group_list[[length(chr_group_list) + 1]] <- data.table(
            level = "chromosome",
            chrom = as.character(chr),
            variable = v,
            grouping = "ancestry",
            test = "Wilcoxon_rank_sum",
            n = nrow(dsub),
            n_groups = ng2,
            statistic = unname(wt$statistic),
            p_value = wt$p.value,
            adj_p_value = NA_real_,
            direction = ""
          )
        } else if (ng2 > 2) {
          kt <- kruskal.test(as.formula(paste(v, "~ ancestry")), data = dsub)
          chr_group_list[[length(chr_group_list) + 1]] <- data.table(
            level = "chromosome",
            chrom = as.character(chr),
            variable = v,
            grouping = "ancestry",
            test = "Kruskal_Wallis",
            n = nrow(dsub),
            n_groups = ng2,
            statistic = unname(kt$statistic),
            p_value = kt$p.value,
            adj_p_value = NA_real_,
            direction = ""
          )

          pw_dt <- pairwise_wilcox_dt(dsub, v, "ancestry", v, "ancestry")
          if (!is.null(pw_dt)) {
            pw_dt[, `:=`(level = "chromosome", chrom = as.character(chr))]
            chr_group_list[[length(chr_group_list) + 1]] <- pw_dt
          }
        }
      }
    }

    if (length(chr_group_list) > 0) {
      chr_group_dt <- rbindlist(chr_group_list, fill = TRUE)
      safe_fwrite(chr_group_dt, file.path(stats_dir, "per_chr_group_comparisons.tsv"))
      cat("Wrote per_chr_group_comparisons.tsv\n")
    }
  }
} else {
  cat("Skipping chromosome-level ROH stats (per_chr_roh_summary.tsv not found)\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. Optional per-chromosome theta stats for every scale
# ─────────────────────────────────────────────────────────────────────────────
theta_files <- get_theta_summary_files(tables_dir)

if (length(theta_files) > 0) {
  for (theta_fp in theta_files) {
    theta_base <- basename(theta_fp)
    scale_label <- sub("^per_chr_theta_summary_", "", theta_base)
    scale_label <- sub("\\.tsv$", "", scale_label)

    theta_d <- fread(theta_fp)
    if ("sample" %in% names(theta_d) && !is.null(anc_file)) {
      anc <- fread(anc_file)
      if (ncol(anc) >= 2) {
        setnames(anc, 1:2, c("sample", "ancestry"))
        theta_d <- merge(theta_d, anc[, .(sample, ancestry)], by = "sample", all.x = TRUE)
      }
    }
    if (!"ancestry" %in% names(theta_d)) theta_d[, ancestry := "all"]

    num_cols <- intersect(
      c(
        "n_windows", "mean_theta", "median_theta", "sd_theta", "min_theta", "max_theta",
        "mean_theta_in_roh", "mean_theta_out_roh", "ratio_theta_in_out"
      ),
      names(theta_d)
    )
    for (cc in num_cols) theta_d[, (cc) := safe_as_num(get(cc))]
    theta_d[, chrom := factor(chrom, levels = natural_chr_order(chrom))]

    # ── 3A. descriptive summary per chromosome ───────────────────────────
    theta_desc <- theta_d[
      ,
      .(
        n_samples = .N,
        mean_mean_theta = mean(mean_theta, na.rm = TRUE),
        median_mean_theta = median(mean_theta, na.rm = TRUE),
        sd_mean_theta = sd(mean_theta, na.rm = TRUE),
        min_mean_theta = min(mean_theta, na.rm = TRUE),
        max_mean_theta = max(mean_theta, na.rm = TRUE),
        mean_ratio_theta_in_out = mean(ratio_theta_in_out, na.rm = TRUE)
      ),
      by = chrom
    ]
    safe_fwrite(
      theta_desc,
      file.path(stats_dir, paste0("per_chr_theta_descriptive_summary_", scale_label, ".tsv"))
    )
    cat("Wrote", paste0("per_chr_theta_descriptive_summary_", scale_label, ".tsv"), "\n")

    # ── 3B. group comparisons per chromosome ─────────────────────────────
    theta_groups <- unique(theta_d$ancestry[!is.na(theta_d$ancestry)])
    if (length(theta_groups) > 1) {
      theta_group_list <- list()
      theta_test_vars <- intersect(
        c("mean_theta", "median_theta", "sd_theta", "mean_theta_in_roh", "mean_theta_out_roh", "ratio_theta_in_out"),
        names(theta_d)
      )

      for (chr in unique(theta_d$chrom)) {
        sub <- theta_d[chrom == chr]
        ng <- length(unique(sub$ancestry[!is.na(sub$ancestry)]))
        if (ng < 2) next

        for (v in theta_test_vars) {
          vals <- sub[[v]]
          ok <- !is.na(vals) & is.finite(vals) & !is.na(sub$ancestry)
          if (sum(ok) < 5) next
          dsub <- sub[ok]
          ng2 <- length(unique(dsub$ancestry))

          if (ng2 == 2) {
            wt <- wilcox.test(as.formula(paste(v, "~ ancestry")), data = dsub)
            theta_group_list[[length(theta_group_list) + 1]] <- data.table(
              level = "chromosome_theta",
              scale = scale_label,
              chrom = as.character(chr),
              variable = v,
              grouping = "ancestry",
              test = "Wilcoxon_rank_sum",
              n = nrow(dsub),
              n_groups = ng2,
              statistic = unname(wt$statistic),
              p_value = wt$p.value,
              adj_p_value = NA_real_,
              direction = ""
            )
          } else if (ng2 > 2) {
            kt <- kruskal.test(as.formula(paste(v, "~ ancestry")), data = dsub)
            theta_group_list[[length(theta_group_list) + 1]] <- data.table(
              level = "chromosome_theta",
              scale = scale_label,
              chrom = as.character(chr),
              variable = v,
              grouping = "ancestry",
              test = "Kruskal_Wallis",
              n = nrow(dsub),
              n_groups = ng2,
              statistic = unname(kt$statistic),
              p_value = kt$p.value,
              adj_p_value = NA_real_,
              direction = ""
            )

            pw_dt <- pairwise_wilcox_dt(dsub, v, "ancestry", v, "ancestry")
            if (!is.null(pw_dt)) {
              pw_dt[, `:=`(
                level = "chromosome_theta",
                scale = scale_label,
                chrom = as.character(chr)
              )]
              theta_group_list[[length(theta_group_list) + 1]] <- pw_dt
            }
          }
        }
      }

      if (length(theta_group_list) > 0) {
        theta_group_dt <- rbindlist(theta_group_list, fill = TRUE)
        safe_fwrite(
          theta_group_dt,
          file.path(stats_dir, paste0("per_chr_theta_group_comparisons_", scale_label, ".tsv"))
        )
        cat("Wrote", paste0("per_chr_theta_group_comparisons_", scale_label, ".tsv"), "\n")
      }
    }

    # ── 3C. theta vs theta/ROH correlations per chromosome ───────────────
    theta_cor_list <- list()
    theta_cor_pairs <- list(
      c("mean_theta", "ratio_theta_in_out"),
      c("mean_theta_in_roh", "mean_theta_out_roh")
    )

    for (chr in unique(theta_d$chrom)) {
      sub <- theta_d[chrom == chr]
      for (pair in theta_cor_pairs) {
        v1 <- pair[1]
        v2 <- pair[2]
        if (!(v1 %in% names(sub)) || !(v2 %in% names(sub))) next
        x <- sub[[v1]]
        y <- sub[[v2]]
        ok <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
        if (sum(ok) < 5) next

        ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
        theta_cor_list[[length(theta_cor_list) + 1]] <- data.table(
          level = "chromosome_theta",
          scale = scale_label,
          chrom = as.character(chr),
          variable = paste0(v1, "_vs_", v2),
          n = sum(ok),
          rho = unname(ct$estimate),
          statistic = unname(ct$statistic),
          p_value = ct$p.value,
          direction = ifelse(unname(ct$estimate) > 0, "positive", "negative")
        )
      }
    }

    if (length(theta_cor_list) > 0) {
      theta_cor_dt <- rbindlist(theta_cor_list, fill = TRUE)
      safe_fwrite(
        theta_cor_dt,
        file.path(stats_dir, paste0("per_chr_theta_correlations_", scale_label, ".tsv"))
      )
      cat("Wrote", paste0("per_chr_theta_correlations_", scale_label, ".tsv"), "\n")
    }
  }
} else {
  cat("Skipping per-chromosome theta stats (no per_chr_theta_summary_*.tsv files found)\n")
}

cat("Stats written to:", stats_dir, "\n")
