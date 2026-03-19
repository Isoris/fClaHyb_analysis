#!/usr/bin/env Rscript
# =============================================================================
# run_stats.R — Simple robust statistics for het/ROH summaries
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

d <- fread(master_file)
d[, het_genomewide := as.numeric(het_genomewide)]
d[, froh := as.numeric(froh)]
d[, roh_total_bp := as.numeric(roh_total_bp)]
d[, longest_roh := as.numeric(longest_roh)]

# ── Merge ancestry if available ──────────────────────────────────────────
if (!is.null(anc_file)) {
  anc <- fread(anc_file)
  if (ncol(anc) >= 2) {
    setnames(anc, 1:2, c("sample", "ancestry"))
    d <- merge(d, anc[, .(sample, ancestry)], by = "sample", all.x = TRUE)
  }
}
if (!"ancestry" %in% names(d)) d[, ancestry := "all"]

# ── 1. Spearman correlations ────────────────────────────────────────────
cor_results <- list()
cor_pairs <- list(
  c("het_genomewide", "froh"),
  c("het_genomewide", "roh_total_bp"),
  c("het_genomewide", "longest_roh")
)

for (pair in cor_pairs) {
  v1 <- pair[1]; v2 <- pair[2]
  x <- d[[v1]]; y <- d[[v2]]
  ok <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
  if (sum(ok) < 5) next
  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  cor_results[[length(cor_results) + 1]] <- data.table(
    var1 = v1, var2 = v2, n = sum(ok),
    rho = ct$estimate, statistic = ct$statistic,
    p_value = ct$p.value, direction = ifelse(ct$estimate > 0, "positive", "negative")
  )
}
if (length(cor_results) > 0) {
  cor_dt <- rbindlist(cor_results)
  fwrite(cor_dt, file.path(stats_dir, "spearman_correlations.tsv"), sep = "\t")
  cat("Wrote spearman_correlations.tsv\n")
}

# ── 2. Group comparisons (if ancestry groups exist) ─────────────────────
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
      # Wilcoxon rank-sum
      wt <- wilcox.test(as.formula(paste(v, "~ ancestry")), data = dsub)
      group_results[[length(group_results) + 1]] <- data.table(
        variable = v, grouping = "ancestry", test = "Wilcoxon_rank_sum",
        n = nrow(dsub), n_groups = ng,
        statistic = wt$statistic, p_value = wt$p.value,
        adj_p_value = NA_real_, direction = ""
      )
    } else if (ng > 2) {
      # Kruskal-Wallis
      kt <- kruskal.test(as.formula(paste(v, "~ ancestry")), data = dsub)
      group_results[[length(group_results) + 1]] <- data.table(
        variable = v, grouping = "ancestry", test = "Kruskal_Wallis",
        n = nrow(dsub), n_groups = ng,
        statistic = kt$statistic, p_value = kt$p.value,
        adj_p_value = NA_real_, direction = ""
      )

      # Post-hoc pairwise Wilcoxon with BH correction
      pw <- pairwise.wilcox.test(dsub[[v]], dsub$ancestry, p.adjust.method = "BH")
      pw_mat <- pw$p.value
      for (i in seq_len(nrow(pw_mat))) {
        for (j in seq_len(ncol(pw_mat))) {
          if (!is.na(pw_mat[i, j])) {
            group_results[[length(group_results) + 1]] <- data.table(
              variable = v,
              grouping = paste0(rownames(pw_mat)[i], "_vs_", colnames(pw_mat)[j]),
              test = "pairwise_Wilcoxon_BH",
              n = nrow(dsub), n_groups = 2,
              statistic = NA_real_, p_value = NA_real_,
              adj_p_value = pw_mat[i, j], direction = ""
            )
          }
        }
      }
    }
  }

  if (length(group_results) > 0) {
    group_dt <- rbindlist(group_results, fill = TRUE)
    fwrite(group_dt, file.path(stats_dir, "group_comparisons.tsv"), sep = "\t")
    cat("Wrote group_comparisons.tsv\n")
  }
}

# ── 3. Descriptive summary ──────────────────────────────────────────────
desc_vars <- c("het_genomewide", "froh", "roh_total_bp", "longest_roh", "n_tracts")
desc_list <- list()
for (v in desc_vars) {
  if (!v %in% names(d)) next
  x <- as.numeric(d[[v]])
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) == 0) next
  desc_list[[length(desc_list) + 1]] <- data.table(
    variable = v, n = length(x),
    mean = mean(x), median = median(x),
    sd = sd(x), min = min(x), max = max(x),
    q25 = quantile(x, 0.25), q75 = quantile(x, 0.75)
  )
}
if (length(desc_list) > 0) {
  desc_dt <- rbindlist(desc_list)
  fwrite(desc_dt, file.path(stats_dir, "descriptive_summary.tsv"), sep = "\t")
  cat("Wrote descriptive_summary.tsv\n")
}

cat("Stats written to:", stats_dir, "\n")
