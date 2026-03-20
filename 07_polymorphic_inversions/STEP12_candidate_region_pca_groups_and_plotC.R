#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop(
    paste(
      "Usage:",
      "Rscript STEP12_candidate_region_pca_groups_and_plotC.R",
      "<candidate_regions.tsv.gz>",
      "<candidate_id>",
      "<sites.tsv.gz>",
      "<dosage.tsv.gz>",
      "<sample_het_windows.tsv.gz>",
      "<outprefix>",
      "<het_sample_col>",
      "<het_value_col>",
      "[het_chrom_col=chrom]",
      "[het_start_col=start_bp]",
      "[het_end_col=end_bp]"
    )
  )
}

candidate_file <- args[1]
candidate_id_value <- args[2]
sites_file <- args[3]
dosage_file <- args[4]
het_file <- args[5]
outprefix <- args[6]
het_sample_col <- args[7]
het_value_col <- args[8]

het_chrom_col <- if (length(args) >= 9) args[9] else "chrom"
het_start_col <- if (length(args) >= 10) args[10] else "start_bp"
het_end_col   <- if (length(args) >= 11) args[11] else "end_bp"

message("[INFO] Reading candidate regions")
cand <- fread(candidate_file)

if (!("candidate_id" %in% names(cand))) {
  stop("candidate_regions file must contain candidate_id")
}

cand_row <- cand[as.character(candidate_id) == as.character(candidate_id_value)]
if (nrow(cand_row) != 1) {
  stop("Expected exactly one candidate row for candidate_id=", candidate_id_value, "; found ", nrow(cand_row))
}

region_chrom <- cand_row$chrom[1]
region_start <- as.numeric(cand_row$start_bp[1])
region_end   <- as.numeric(cand_row$end_bp[1])

message("[INFO] Candidate region: ", region_chrom, ":", region_start, "-", region_end)

message("[INFO] Reading sites")
sites <- fread(sites_file)
req_sites <- c("marker", "chrom", "pos", "allele1", "allele2")
miss_sites <- setdiff(req_sites, names(sites))
if (length(miss_sites) > 0) {
  stop("Sites file missing required columns: ", paste(miss_sites, collapse = ", "))
}

message("[INFO] Reading dosage")
dos <- fread(dosage_file)
if (!("marker" %in% names(dos))) stop("Dosage file must contain marker column")

if (!identical(sites$marker, dos$marker)) {
  setkeyv(sites, "marker")
  setkeyv(dos, "marker")
  dos <- dos[sites$marker]
  if (!identical(sites$marker, dos$marker)) stop("Sites and dosage markers do not match")
}

sample_names <- setdiff(names(dos), "marker")

# keep only SNPs in candidate region
keep <- which(sites$chrom == region_chrom & sites$pos >= region_start & sites$pos <= region_end)
if (length(keep) < 10) {
  stop("Too few SNPs in candidate region after filtering: ", length(keep))
}

sites_reg <- sites[keep]
dos_reg <- dos[keep]

X <- as.matrix(dos_reg[, ..sample_names])
storage.mode(X) <- "double"

# PCA on individuals: transpose so rows = individuals, cols = SNPs
# center columns (SNPs)
pca <- prcomp(t(X), center = TRUE, scale. = FALSE)

pcs <- data.table(
  sample = sample_names,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

# kmeans into 3 groups on PC1 only
set.seed(1)
km <- kmeans(matrix(pcs$PC1, ncol = 1), centers = 3, nstart = 50)

pcs[, raw_group := km$cluster]

# order groups by PC1 mean so group1=left, group2=middle, group3=right
grp_order <- pcs[, .(pc1_mean = mean(PC1)), by = raw_group][order(pc1_mean)]
grp_order[, ordered_group := 1:.N]
pcs <- merge(pcs, grp_order[, .(raw_group, ordered_group)], by = "raw_group", all.x = TRUE)
setorder(pcs, match(sample, sample_names))

# regional heterozygosity from dosage:
# heterozygosity per SNP per sample ≈ dosage * (2 - dosage)
# this is 0 at 0/2 and max at 1
H <- X * (2 - X)
regional_het <- colMeans(H, na.rm = TRUE)

pcs[, regional_het := regional_het[match(sample, sample_names)]]

# Optional SNP orientation to PC1 axis
# For each SNP, if correlation(dosage, PC1) < 0, flip dosage to 2-dosage
pc1_vec <- pcs$PC1[match(sample_names, pcs$sample)]

flip_flag <- rep(FALSE, nrow(X))
for (i in seq_len(nrow(X))) {
  xi <- X[i, ]
  ok <- is.finite(xi) & is.finite(pc1_vec)
  if (sum(ok) >= 3) {
    rr <- suppressWarnings(cor(xi[ok], pc1_vec[ok]))
    if (is.finite(rr) && rr < 0) flip_flag[i] <- TRUE
  }
}

X_oriented <- X
X_oriented[flip_flag, ] <- 2 - X_oriented[flip_flag, ]

regional_hap_score <- colMeans(X_oriented, na.rm = TRUE)
pcs[, regional_hap_score := regional_hap_score[match(sample, sample_names)]]

# Make nice group labels
pcs[, group_label := fifelse(
  ordered_group == 1, "G1_left",
  fifelse(ordered_group == 2, "G2_middle", "G3_right")
)]

# Highlight representative individuals nearest group centers in PCA space
centers <- pcs[, .(PC1c = mean(PC1), PC2c = mean(PC2)), by = group_label]
pcs2 <- merge(pcs, centers, by = "group_label", all.x = TRUE)
pcs2[, dist_center := sqrt((PC1 - PC1c)^2 + (PC2 - PC2c)^2)]
reps <- pcs2[order(dist_center), .SD[1], by = group_label]

# Read chromosome-wide heterozygosity windows
message("[INFO] Reading sample heterozygosity windows")
het <- fread(het_file)

req_het <- c(het_sample_col, het_value_col, het_chrom_col, het_start_col, het_end_col)
miss_het <- setdiff(req_het, names(het))
if (length(miss_het) > 0) {
  stop("Heterozygosity window file missing required columns: ", paste(miss_het, collapse = ", "))
}

setnames(het, het_sample_col, "sample")
setnames(het, het_value_col, "het_value")
setnames(het, het_chrom_col, "chrom")
setnames(het, het_start_col, "start_bp")
setnames(het, het_end_col, "end_bp")

het[, start_bp := as.numeric(start_bp)]
het[, end_bp := as.numeric(end_bp)]
het[, mid_bp := (start_bp + end_bp) / 2]

# merge group labels onto heterozygosity windows
het2 <- merge(
  het,
  pcs[, .(sample, group_label, ordered_group, regional_het, regional_hap_score)],
  by = "sample",
  all.x = FALSE
)

# keep only chromosome of interest for plot C
het_chr <- het2[chrom == region_chrom]

# summarize mean heterozygosity by group and window midpoint
plotc_dt <- het_chr[
  ,
  .(
    het_mean = mean(het_value, na.rm = TRUE),
    het_sd = sd(het_value, na.rm = TRUE),
    n = .N
  ),
  by = .(group_label, ordered_group, mid_bp)
][order(ordered_group, mid_bp)]

# ---- outputs tables ----
pcs_out <- paste0(outprefix, ".candidate_", candidate_id_value, ".regional_pca_samples.tsv.gz")
sites_out <- paste0(outprefix, ".candidate_", candidate_id_value, ".regional_sites.tsv.gz")
plotc_out <- paste0(outprefix, ".candidate_", candidate_id_value, ".plotC_summary.tsv.gz")
rep_out <- paste0(outprefix, ".candidate_", candidate_id_value, ".representatives.tsv.gz")

fwrite(pcs, pcs_out, sep = "\t")
fwrite(sites_reg, sites_out, sep = "\t")
fwrite(plotc_dt, plotc_out, sep = "\t")
fwrite(reps, rep_out, sep = "\t")

# ---- plots ----

# Panel A: plain PCA
p_a <- ggplot(pcs, aes(x = PC1, y = PC2)) +
  geom_point(shape = 1, alpha = 0.6, size = 2) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("Candidate ", candidate_id_value, " regional PCA"),
    subtitle = paste0(region_chrom, ":", format(region_start, big.mark=","), "-", format(region_end, big.mark=",")),
    x = "PC1",
    y = "PC2"
  )

# Panel B: PCA colored by regional heterozygosity, labeled by group number
p_b <- ggplot(pcs, aes(x = PC1, y = PC2, color = regional_het, label = ordered_group)) +
  geom_text(size = 3) +
  geom_point(
    data = reps,
    aes(x = PC1, y = PC2),
    inherit.aes = FALSE,
    shape = 1,
    size = 6,
    stroke = 1.2,
    color = "black"
  ) +
  theme_bw(base_size = 12) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = paste0("Candidate ", candidate_id_value, " PCA colored by regional heterozygosity"),
    subtitle = "Circled points are representative individuals nearest each group center",
    x = "PC1",
    y = "PC2",
    color = "Regional\nheterozygosity"
  )

# Panel C: chromosome-wide heterozygosity curves by group
p_c <- ggplot(plotc_dt, aes(x = mid_bp / 1e6, y = het_mean, color = group_label)) +
  geom_line(linewidth = 1) +
  annotate(
    "rect",
    xmin = region_start / 1e6,
    xmax = region_end / 1e6,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.15,
    fill = "grey60"
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0(region_chrom, " heterozygosity by inferred regional PCA group"),
    subtitle = paste0("Shaded region = candidate ", candidate_id_value),
    x = "Chromosome position (Mb)",
    y = "Mean heterozygosity"
  )

# save plots
pdf_a <- paste0(outprefix, ".candidate_", candidate_id_value, ".panelA_regional_PCA.pdf")
pdf_b <- paste0(outprefix, ".candidate_", candidate_id_value, ".panelB_regional_PCA_het.pdf")
pdf_c <- paste0(outprefix, ".candidate_", candidate_id_value, ".panelC_chrom_het_by_group.pdf")

ggsave(pdf_a, p_a, width = 6, height = 5)
ggsave(pdf_b, p_b, width = 7, height = 5.5)
ggsave(pdf_c, p_c, width = 7, height = 4.8)

message("[DONE] Wrote:")
message("  ", pcs_out)
message("  ", sites_out)
message("  ", plotc_out)
message("  ", rep_out)
message("  ", pdf_a)
message("  ", pdf_b)
message("  ", pdf_c)
