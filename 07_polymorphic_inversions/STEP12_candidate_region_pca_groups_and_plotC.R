#!/usr/bin/env Rscript

# =============================================================================
# STEP12_candidate_region_pca_groups_and_plotC.R
#
# For one chosen candidate region:
#   - extract SNPs from the regional dosage matrix
#   - run a regional PCA on individuals
#   - assign 3 genotype-like groups using k-means on PC1
#   - compute regional heterozygosity per individual from dosage
#   - optionally orient SNP dosages to the regional PC1 axis
#   - merge group labels with chromosome-wide tP diversity windows
#   - generate panels A, B, C
#
# Usage:
#   Rscript STEP12_candidate_region_pca_groups_and_plotC.R \
#     <candidate_regions.tsv.gz> \
#     <candidate_id> \
#     <sites.tsv.gz> \
#     <dosage.tsv.gz> \
#     <sample_tP_windows.tsv.gz> \
#     <outprefix> \
#     [n_groups=3]
#
# The sample_tP_windows file is the output of STEP10b, with columns:
#   sample, chrom, WinStart, WinStop, WinCenter, tP, nSites
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop(paste(
    "Usage: Rscript STEP12_candidate_region_pca_groups_and_plotC.R",
    "<candidate_regions.tsv.gz> <candidate_id> <sites.tsv.gz>",
    "<dosage.tsv.gz> <sample_tP_windows.tsv.gz> <outprefix>",
    "[n_groups=3]"
  ))
}

candidate_file     <- args[1]
candidate_id_value <- args[2]
sites_file         <- args[3]
dosage_file        <- args[4]
tP_windows_file    <- args[5]
outprefix          <- args[6]
n_groups           <- if (length(args) >= 7) as.integer(args[7]) else 3L

# ── Read candidate region ──────────────────────────────────────────────────
message("[INFO] Reading candidate regions")
cand <- fread(candidate_file)

if (!("candidate_id" %in% names(cand))) {
  stop("candidate_regions file must contain 'candidate_id' column")
}

cand_row <- cand[as.character(candidate_id) == as.character(candidate_id_value)]
if (nrow(cand_row) != 1) {
  stop("Expected exactly one candidate row for candidate_id=", candidate_id_value,
       "; found ", nrow(cand_row))
}

region_chrom <- cand_row$chrom[1]
region_start <- as.numeric(cand_row$start_bp[1])
region_end   <- as.numeric(cand_row$end_bp[1])

message("[INFO] Candidate region: ", region_chrom, ":",
        format(region_start, big.mark = ","), "-",
        format(region_end, big.mark = ","))

# ── Read sites + dosage ───────────────────────────────────────────────────
message("[INFO] Reading sites")
sites <- fread(sites_file)
req_sites <- c("marker", "chrom", "pos", "allele1", "allele2")
miss_sites <- setdiff(req_sites, names(sites))
if (length(miss_sites) > 0) {
  stop("Sites file missing required columns: ", paste(miss_sites, collapse = ", "))
}

message("[INFO] Reading dosage")
dos <- fread(dosage_file)
if (!("marker" %in% names(dos))) stop("Dosage file must contain 'marker' column")

# Align
if (!identical(sites$marker, dos$marker)) {
  setkeyv(sites, "marker")
  setkeyv(dos, "marker")
  dos <- dos[sites$marker]
  if (!identical(sites$marker, dos$marker)) stop("Sites and dosage markers do not match")
}

sample_names <- setdiff(names(dos), "marker")
message("[INFO] Samples: ", length(sample_names))

# ── Filter to candidate region ─────────────────────────────────────────────
keep <- which(sites$chrom == region_chrom &
              sites$pos >= region_start &
              sites$pos <= region_end)
if (length(keep) < 10) {
  stop("Too few SNPs in candidate region after filtering: ", length(keep))
}
message("[INFO] SNPs in candidate region: ", length(keep))

sites_reg <- sites[keep]
dos_reg   <- dos[keep]

X <- as.matrix(dos_reg[, ..sample_names])
storage.mode(X) <- "double"

# ── Regional PCA on individuals ────────────────────────────────────────────
# t(X): rows = individuals, cols = SNPs
pca <- prcomp(t(X), center = TRUE, scale. = FALSE)

pcs <- data.table(
  sample = sample_names,
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2]
)

# Variance explained
var_expl <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1_pct <- round(var_expl[1], 1)
pc2_pct <- round(var_expl[2], 1)

# ── k-means clustering on PC1 ─────────────────────────────────────────────
set.seed(1)
km <- kmeans(matrix(pcs$PC1, ncol = 1), centers = n_groups, nstart = 50)
pcs[, raw_group := km$cluster]

# Order groups by PC1 mean: group 1 = leftmost, group N = rightmost
grp_order <- pcs[, .(pc1_mean = mean(PC1)), by = raw_group][order(pc1_mean)]
grp_order[, ordered_group := seq_len(.N)]
pcs <- merge(pcs, grp_order[, .(raw_group, ordered_group)], by = "raw_group", all.x = TRUE)

# Restore original sample order
pcs <- pcs[match(sample_names, sample)]

# Group labels
group_labels <- c("Homo_1", "Het", "Homo_2")
if (n_groups != 3) {
  group_labels <- paste0("G", seq_len(n_groups))
}
pcs[, group_label := group_labels[ordered_group]]
pcs[, group_label := factor(group_label, levels = group_labels)]

message("[INFO] Group sizes: ", paste(pcs[, .N, by = group_label][order(group_label), paste0(group_label, "=", N)], collapse = ", "))

# ── Regional heterozygosity from dosage ────────────────────────────────────
# H per SNP per sample ≈ dosage × (2 - dosage)  (0 at homo, max at het)
H <- X * (2 - X)
regional_het <- colMeans(H, na.rm = TRUE)
pcs[, regional_het := regional_het[match(sample, sample_names)]]

# ── Optional SNP orientation to PC1 axis ──────────────────────────────────
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

# ── Representative individuals (nearest to group centroid) ────────────────
centers <- pcs[, .(PC1c = mean(PC1), PC2c = mean(PC2)), by = group_label]
pcs2 <- merge(pcs, centers, by = "group_label", all.x = TRUE)
pcs2[, dist_center := sqrt((PC1 - PC1c)^2 + (PC2 - PC2c)^2)]
reps <- pcs2[order(dist_center), .SD[1], by = group_label]

# Carry region coordinates into reps for STEP13
reps[, candidate_id       := as.integer(candidate_id_value)]
reps[, candidate_chrom    := region_chrom]
reps[, candidate_start_bp := region_start]
reps[, candidate_end_bp   := region_end]

# ── Read per-sample tP windows ────────────────────────────────────────────
message("[INFO] Reading per-sample tP windows")
tP_all <- fread(tP_windows_file)

req_tP <- c("sample", "chrom", "WinStart", "WinStop", "WinCenter", "tP")
miss_tP <- setdiff(req_tP, names(tP_all))
if (length(miss_tP) > 0) {
  stop("tP windows file missing required columns: ", paste(miss_tP, collapse = ", "))
}

# Keep only samples in our analysis and the chromosome of interest
tP_chr <- tP_all[sample %in% sample_names & chrom == region_chrom]

if (nrow(tP_chr) == 0) {
  stop("No tP data found for chromosome ", region_chrom, " and the analysis samples")
}

message("[INFO] tP windows for this chromosome: ", nrow(tP_chr),
        " (", uniqueN(tP_chr$sample), " samples)")

# Merge group labels onto tP windows
tP_chr <- merge(
  tP_chr,
  pcs[, .(sample, group_label, ordered_group, regional_het, regional_hap_score)],
  by = "sample",
  all.x = FALSE
)

# ── Summarize tP by group and window for plot C ──────────────────────────
plotc_dt <- tP_chr[
  is.finite(tP),
  .(
    tP_mean = mean(tP, na.rm = TRUE),
    tP_sd   = sd(tP, na.rm = TRUE),
    n       = .N
  ),
  by = .(group_label, ordered_group, WinCenter)
][order(ordered_group, WinCenter)]

# ── Write output tables ───────────────────────────────────────────────────
cid <- candidate_id_value
pcs_out   <- paste0(outprefix, ".candidate_", cid, ".regional_pca_samples.tsv.gz")
sites_out <- paste0(outprefix, ".candidate_", cid, ".regional_sites.tsv.gz")
plotc_out <- paste0(outprefix, ".candidate_", cid, ".plotC_summary.tsv.gz")
rep_out   <- paste0(outprefix, ".candidate_", cid, ".representatives.tsv.gz")

fwrite(pcs, pcs_out, sep = "\t")
fwrite(sites_reg, sites_out, sep = "\t")
fwrite(plotc_dt, plotc_out, sep = "\t")
fwrite(reps, rep_out, sep = "\t")

# ── Panel A: plain regional PCA ──────────────────────────────────────────
p_a <- ggplot(pcs, aes(x = PC1, y = PC2)) +
  geom_point(shape = 1, alpha = 0.6, size = 2) +
  theme_bw(base_size = 12) +
  labs(
    title    = paste0("Candidate ", cid, " — regional PCA"),
    subtitle = paste0(region_chrom, ": ",
                      format(region_start, big.mark = ","), " – ",
                      format(region_end, big.mark = ","),
                      "  (", length(keep), " SNPs)"),
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)")
  )

# ── Panel B: PCA colored by regional heterozygosity, labeled by group ────
p_b <- ggplot(pcs, aes(x = PC1, y = PC2, color = regional_het)) +
  geom_point(aes(shape = group_label), size = 2.5, alpha = 0.8) +
  geom_point(
    data = reps,
    aes(x = PC1, y = PC2),
    inherit.aes = FALSE,
    shape = 1, size = 6, stroke = 1.2, color = "black"
  ) +
  geom_text(
    data = reps,
    aes(x = PC1, y = PC2, label = group_label),
    inherit.aes = FALSE,
    nudge_y = diff(range(pcs$PC2)) * 0.06,
    color = "black", size = 3.5, fontface = "bold"
  ) +
  theme_bw(base_size = 12) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title    = paste0("Candidate ", cid, " — PCA colored by dosage heterozygosity"),
    subtitle = "Circled = representative individual nearest each group center",
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)"),
    color = "Regional\nhet (dosage)",
    shape = "Group"
  )

# ── Panel C: chromosome-wide tP by group ─────────────────────────────────
p_c <- ggplot(plotc_dt, aes(x = WinCenter / 1e6, y = tP_mean, color = group_label)) +
  annotate(
    "rect",
    xmin = region_start / 1e6, xmax = region_end / 1e6,
    ymin = -Inf, ymax = Inf,
    alpha = 0.15, fill = "grey60"
  ) +
  geom_line(linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title    = paste0(region_chrom, " — pairwise θ (tP) by inferred PCA group"),
    subtitle = paste0("Shaded region = candidate ", cid),
    x = "Chromosome position (Mb)",
    y = "Mean tP (pairwise θ)",
    color = "Group"
  )

# ── Save plots ────────────────────────────────────────────────────────────
pdf_a <- paste0(outprefix, ".candidate_", cid, ".panelA_regional_PCA.pdf")
pdf_b <- paste0(outprefix, ".candidate_", cid, ".panelB_regional_PCA_het.pdf")
pdf_c <- paste0(outprefix, ".candidate_", cid, ".panelC_chrom_tP_by_group.pdf")

ggsave(pdf_a, p_a, width = 6, height = 5)
ggsave(pdf_b, p_b, width = 7.5, height = 5.5)
ggsave(pdf_c, p_c, width = 8, height = 4.8)

message("[DONE] Wrote:")
message("  ", pcs_out)
message("  ", sites_out)
message("  ", plotc_out)
message("  ", rep_out)
message("  ", pdf_a)
message("  ", pdf_b)
message("  ", pdf_c)
