#!/usr/bin/env Rscript
# v2 note: added per-chr descriptive, correlations, group comparisons;
#          added theta-scale-aware chr stats from per_chr_theta_summary_*.tsv
# =============================================================================
# run_stats.R — Robust statistics for het/ROH summaries (sample + chromosome)
# =============================================================================
# Usage:
#   Rscript run_stats.R <master_summary.tsv> <stats_dir> [ancestry_labels.tsv] \
#     [per_chr_roh_summary.tsv] [breeding_tables_dir]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: run_stats.R <master.tsv> <stats_dir> [ancestry] [per_chr] [breeding_dir]")

master_file <- args[1]
stats_dir   <- args[2]
anc_file    <- if (length(args) >= 3 && file.exists(args[3])) args[3] else NULL
chr_file    <- if (length(args) >= 4 && file.exists(args[4])) args[4] else NULL
breed_dir   <- if (length(args) >= 5 && dir.exists(args[5])) args[5] else NULL

dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

d <- fread(master_file)
for (col in c("het_genomewide","froh","roh_total_bp","longest_roh")) {
  if (col %in% names(d)) d[, (col) := as.numeric(get(col))]
}

# ── Merge ancestry ───────────────────────────────────────────────────────
if (!is.null(anc_file)) {
  anc <- fread(anc_file)
  if (ncol(anc) >= 2) {
    setnames(anc, 1:2, c("sample", "ancestry"))
    d <- merge(d, anc[, .(sample, ancestry)], by = "sample", all.x = TRUE)
  }
}
if (!"ancestry" %in% names(d)) d[, ancestry := "all"]

# ═══════════════════════════════════════════════════════════════════════════
# ORIGINAL SAMPLE-LEVEL STATS (unchanged)
# ═══════════════════════════════════════════════════════════════════════════

# ── 1. Spearman correlations ────────────────────────────────────────────
cor_results <- list()
cor_pairs <- list(c("het_genomewide","froh"), c("het_genomewide","roh_total_bp"),
                  c("het_genomewide","longest_roh"))
for (pair in cor_pairs) {
  v1 <- pair[1]; v2 <- pair[2]
  if (!all(c(v1,v2) %in% names(d))) next
  x <- d[[v1]]; y <- d[[v2]]
  ok <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
  if (sum(ok) < 5) next
  ct <- cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
  cor_results[[length(cor_results)+1]] <- data.table(
    var1=v1, var2=v2, n=sum(ok), rho=ct$estimate, statistic=ct$statistic,
    p_value=ct$p.value, direction=ifelse(ct$estimate>0,"positive","negative"))
}
if (length(cor_results) > 0) {
  fwrite(rbindlist(cor_results), file.path(stats_dir, "spearman_correlations.tsv"), sep="\t")
  cat("Wrote spearman_correlations.tsv\n")
}

# ── 2. Group comparisons ────────────────────────────────────────────────
groups <- unique(d$ancestry[!is.na(d$ancestry)])
if (length(groups) > 1) {
  group_results <- list()
  test_vars <- c("het_genomewide","froh","roh_total_bp","longest_roh","n_tracts")
  for (v in test_vars) {
    if (!v %in% names(d)) next
    vals <- as.numeric(d[[v]])
    ok <- !is.na(vals) & is.finite(vals) & !is.na(d$ancestry)
    if (sum(ok) < 5) next
    dsub <- d[ok]
    ng <- length(unique(dsub$ancestry))
    if (ng == 2) {
      wt <- wilcox.test(as.formula(paste(v, "~ ancestry")), data=dsub)
      group_results[[length(group_results)+1]] <- data.table(
        variable=v, grouping="ancestry", test="Wilcoxon_rank_sum",
        n=nrow(dsub), n_groups=ng, statistic=wt$statistic, p_value=wt$p.value,
        adj_p_value=NA_real_, direction="")
    } else if (ng > 2) {
      kt <- kruskal.test(as.formula(paste(v, "~ ancestry")), data=dsub)
      group_results[[length(group_results)+1]] <- data.table(
        variable=v, grouping="ancestry", test="Kruskal_Wallis",
        n=nrow(dsub), n_groups=ng, statistic=kt$statistic, p_value=kt$p.value,
        adj_p_value=NA_real_, direction="")
      pw <- pairwise.wilcox.test(dsub[[v]], dsub$ancestry, p.adjust.method="BH")
      pm <- pw$p.value
      for (i in seq_len(nrow(pm))) for (j in seq_len(ncol(pm))) {
        if (!is.na(pm[i,j])) {
          group_results[[length(group_results)+1]] <- data.table(
            variable=v, grouping=paste0(rownames(pm)[i],"_vs_",colnames(pm)[j]),
            test="pairwise_Wilcoxon_BH", n=nrow(dsub), n_groups=2,
            statistic=NA_real_, p_value=NA_real_, adj_p_value=pm[i,j], direction="")
        }
      }
    }
  }
  if (length(group_results) > 0) {
    fwrite(rbindlist(group_results, fill=TRUE), file.path(stats_dir, "group_comparisons.tsv"), sep="\t")
    cat("Wrote group_comparisons.tsv\n")
  }
}

# ── 3. Descriptive summary ──────────────────────────────────────────────
desc_vars <- c("het_genomewide","froh","roh_total_bp","longest_roh","n_tracts")
desc_list <- list()
for (v in desc_vars) {
  if (!v %in% names(d)) next
  x <- as.numeric(d[[v]]); x <- x[!is.na(x) & is.finite(x)]
  if (length(x) == 0) next
  desc_list[[length(desc_list)+1]] <- data.table(
    variable=v, n=length(x), mean=mean(x), median=median(x),
    sd=sd(x), min=min(x), max=max(x), q25=quantile(x,0.25), q75=quantile(x,0.75))
}
if (length(desc_list) > 0) {
  fwrite(rbindlist(desc_list), file.path(stats_dir, "descriptive_summary.tsv"), sep="\t")
  cat("Wrote descriptive_summary.tsv\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# v2: CHROMOSOME-LEVEL STATS
# ═══════════════════════════════════════════════════════════════════════════
if (!is.null(chr_file)) {
  chr_d <- fread(chr_file)
  for (col in c("roh_bp","froh_chr","n_tracts","longest_tract","callable_bp_chr")) {
    if (col %in% names(chr_d)) chr_d[, (col) := as.numeric(get(col))]
  }
  # Merge ancestry
  if (!is.null(anc_file)) {
    chr_d <- merge(chr_d, unique(d[, .(sample, ancestry)]), by="sample", all.x=TRUE)
  }
  if (!"ancestry" %in% names(chr_d)) chr_d[, ancestry := "all"]

  # Per-chr descriptive summary
  chr_desc <- list()
  for (v in c("froh_chr","roh_bp","n_tracts","longest_tract")) {
    if (!v %in% names(chr_d)) next
    tmp <- chr_d[!is.na(get(v)), .(
      n=.N, mean=mean(get(v),na.rm=TRUE), median=median(get(v),na.rm=TRUE),
      sd=sd(get(v),na.rm=TRUE), min=min(get(v),na.rm=TRUE), max=max(get(v),na.rm=TRUE)
    ), by=chrom]
    tmp[, variable := v]
    chr_desc[[length(chr_desc)+1]] <- tmp
  }
  if (length(chr_desc) > 0) {
    fwrite(rbindlist(chr_desc, fill=TRUE), file.path(stats_dir, "per_chr_descriptive_summary.tsv"), sep="\t")
    cat("Wrote per_chr_descriptive_summary.tsv\n")
  }

  # Per-chr group comparisons
  if (length(groups) > 1 && "froh_chr" %in% names(chr_d)) {
    chr_grp <- list()
    for (chrom_val in unique(chr_d$chrom)) {
      dsub <- chr_d[chrom == chrom_val & !is.na(froh_chr) & !is.na(ancestry)]
      ng <- length(unique(dsub$ancestry))
      if (nrow(dsub) < 5 || ng < 2) next
      if (ng == 2) {
        wt <- tryCatch(wilcox.test(froh_chr ~ ancestry, data=dsub), error=function(e) NULL)
        if (!is.null(wt)) {
          chr_grp[[length(chr_grp)+1]] <- data.table(
            chrom=chrom_val, variable="froh_chr", test="Wilcoxon", n=nrow(dsub),
            n_groups=ng, statistic=wt$statistic, p_value=wt$p.value)
        }
      } else {
        kt <- tryCatch(kruskal.test(froh_chr ~ ancestry, data=dsub), error=function(e) NULL)
        if (!is.null(kt)) {
          chr_grp[[length(chr_grp)+1]] <- data.table(
            chrom=chrom_val, variable="froh_chr", test="Kruskal_Wallis", n=nrow(dsub),
            n_groups=ng, statistic=kt$statistic, p_value=kt$p.value)
        }
      }
    }
    if (length(chr_grp) > 0) {
      fwrite(rbindlist(chr_grp, fill=TRUE), file.path(stats_dir, "per_chr_group_comparisons.tsv"), sep="\t")
      cat("Wrote per_chr_group_comparisons.tsv\n")
    }
  }
}

# ═══════════════════════════════════════════════════════════════════════════
# v2: THETA-SCALE CHROMOSOME STATS (scan for per_chr_theta_summary_*.tsv)
# ═══════════════════════════════════════════════════════════════════════════
if (!is.null(breed_dir)) {
  theta_files <- list.files(breed_dir, pattern="^per_chr_theta_summary_.*\\.tsv$", full.names=TRUE)
  for (tf in theta_files) {
    scale_label <- sub(".*per_chr_theta_summary_(.*)\\.tsv$", "\\1", basename(tf))
    td <- fread(tf)
    for (col in c("mean_theta","median_theta","sd_theta","frac_low_theta_windows")) {
      if (col %in% names(td)) td[, (col) := as.numeric(get(col))]
    }

    # Descriptive per chr
    td_desc <- td[!is.na(mean_theta), .(
      n_samples=.N, mean_mean_theta=mean(mean_theta,na.rm=TRUE),
      median_mean_theta=median(mean_theta,na.rm=TRUE),
      sd_mean_theta=sd(mean_theta,na.rm=TRUE),
      mean_frac_low=mean(as.numeric(frac_low_theta_windows),na.rm=TRUE)
    ), by=chrom]
    fwrite(td_desc, file.path(stats_dir, paste0("per_chr_theta_descriptive_summary_",scale_label,".tsv")), sep="\t")
    cat("Wrote per_chr_theta_descriptive_summary_", scale_label, ".tsv\n", sep="")

    # Theta correlations per chr (mean_theta vs frac_low)
    if (all(c("mean_theta","frac_low_theta_windows") %in% names(td))) {
      td[, frac_low_theta_windows := as.numeric(frac_low_theta_windows)]
      tcor <- list()
      for (chrom_val in unique(td$chrom)) {
        dsub <- td[chrom == chrom_val & !is.na(mean_theta) & !is.na(frac_low_theta_windows)]
        if (nrow(dsub) < 5) next
        ct <- tryCatch(cor.test(dsub$mean_theta, dsub$frac_low_theta_windows,
                                method="spearman", exact=FALSE), error=function(e) NULL)
        if (!is.null(ct)) {
          tcor[[length(tcor)+1]] <- data.table(
            chrom=chrom_val, var1="mean_theta", var2="frac_low_theta_windows",
            n=nrow(dsub), rho=ct$estimate, p_value=ct$p.value)
        }
      }
      if (length(tcor) > 0) {
        fwrite(rbindlist(tcor), file.path(stats_dir, paste0("per_chr_theta_correlations_",scale_label,".tsv")), sep="\t")
        cat("Wrote per_chr_theta_correlations_", scale_label, ".tsv\n", sep="")
      }
    }
  }
}

cat("Stats complete:", stats_dir, "\n")
