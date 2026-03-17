#!/usr/bin/env Rscript
###############################################################################
# plot_relatedness_3panel.R
#
# Produces a 3-panel figure from ngsRelate output:
#   A) Kinship heatmap (continuous theta values)
#   B) Relationship category heatmap (discrete classes)
#   C) Bar plot of pairwise relationship counts
#
# Input:  ngsRelate output file (catfish_226_relatedness.res)
# Output: PDF and PNG of the 3-panel figure
#
# Usage:
#   Rscript plot_relatedness_3panel.R \
#     --input /path/to/catfish_226_relatedness.res \
#     --outdir /path/to/output_directory \
#     --prefix catfish_relatedness
#
# Or source in R and call:
#   plot_relatedness(input_file, outdir, prefix)
###############################################################################

# ── Load libraries ──────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(patchwork)
  library(optparse)
})

# ── Command-line arguments ──────────────────────────────────────────────────
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "ngsRelate output file (.res)"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Output directory [default: .]"),
  make_option(c("-p", "--prefix"), type = "character", default = "relatedness",
              help = "Output file prefix [default: relatedness]"),
  make_option(c("--width"), type = "numeric", default = 20,
              help = "Figure width in inches [default: 20]"),
  make_option(c("--height"), type = "numeric", default = 8,
              help = "Figure height in inches [default: 8]"),
  make_option(c("--max_label"), type = "integer", default = 60,
              help = "Max samples to show axis labels [default: 60]. Set 0 to always show."),
  make_option(c("--theta_col"), type = "integer", default = 5,
              help = "1-indexed column number for theta in ngsRelate output [default: 5]"),
  make_option(c("--ida_col"), type = "integer", default = 3,
              help = "1-indexed column for sample ID a [default: 3]"),
  make_option(c("--idb_col"), type = "integer", default = 4,
              help = "1-indexed column for sample ID b [default: 4]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("--input is required", call. = FALSE)
}

# ── Read ngsRelate output ───────────────────────────────────────────────────
cat("[INFO] Reading ngsRelate output:", opt$input, "\n")

res <- read.table(opt$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("[INFO] Columns detected:\n")
cat(paste0("  ", seq_along(colnames(res)), ": ", colnames(res), "\n"))
cat("[INFO] Number of pairwise comparisons:", nrow(res), "\n")

# ── Extract relevant columns ───────────────────────────────────────────────
# Adjust column indices (user provides 1-indexed)
ida_col <- opt$ida_col
idb_col <- opt$idb_col
theta_col <- opt$theta_col

df <- data.frame(
  ind1  = as.character(res[, ida_col]),
  ind2  = as.character(res[, idb_col]),
  theta = as.numeric(res[, theta_col]),
  stringsAsFactors = FALSE
)

cat("[INFO] Sample of extracted data:\n")
print(head(df))

# ── Classify relationships ──────────────────────────────────────────────────
# Standard KING kinship thresholds applied to theta
classify_relationship <- function(theta) {
  case_when(
    theta >= 0.354   ~ "Duplicate/MZ",
    theta >= 0.177   ~ "First degree",
    theta >= 0.0884  ~ "Second degree",
    theta >= 0.0442  ~ "Third degree",
    TRUE             ~ "Unrelated"
  )
}

df$relationship <- classify_relationship(df$theta)

# Order factor for consistent plotting
rel_levels <- c("Duplicate/MZ", "First degree", "Second degree",
                "Third degree", "Unrelated")
df$relationship <- factor(df$relationship, levels = rel_levels)

cat("[INFO] Relationship counts:\n")
print(table(df$relationship))

# ── Get unique sample list ──────────────────────────────────────────────────
samples <- sort(unique(c(df$ind1, df$ind2)))
n_samples <- length(samples)
cat("[INFO] Number of unique samples:", n_samples, "\n")

# ── Build symmetric matrices ────────────────────────────────────────────────
cat("[INFO] Building kinship matrix...\n")

# Initialize matrices
theta_mat <- matrix(NA, nrow = n_samples, ncol = n_samples,
                    dimnames = list(samples, samples))
rel_mat   <- matrix(NA, nrow = n_samples, ncol = n_samples,
                    dimnames = list(samples, samples))

# Fill upper + lower triangle
for (i in seq_len(nrow(df))) {
  a <- df$ind1[i]
  b <- df$ind2[i]
  theta_mat[a, b] <- df$theta[i]
  theta_mat[b, a] <- df$theta[i]
  rel_mat[a, b]   <- as.character(df$relationship[i])
  rel_mat[b, a]   <- as.character(df$relationship[i])
}

# Diagonal = self (theta ~ 0.5 for kinship, or NA)
diag(theta_mat) <- 0.5
diag(rel_mat)   <- "Self"

# ── Optional: cluster samples by kinship for better visualization ───────────
cat("[INFO] Clustering samples by kinship...\n")

# Replace NA with 0 for clustering
theta_clust <- theta_mat
theta_clust[is.na(theta_clust)] <- 0

# Distance = 1 - theta (more related = closer)
d <- as.dist(1 - theta_clust)
hc <- hclust(d, method = "average")
sample_order <- hc$labels[hc$order]

# ── Melt matrices for ggplot ────────────────────────────────────────────────
cat("[INFO] Preparing plot data...\n")

# Theta heatmap data
theta_melt <- melt(theta_mat, varnames = c("Ind1", "Ind2"), value.name = "Theta")
theta_melt$Ind1 <- factor(theta_melt$Ind1, levels = sample_order)
theta_melt$Ind2 <- factor(theta_melt$Ind2, levels = rev(sample_order))

# Relationship heatmap data
rel_melt <- melt(rel_mat, varnames = c("Ind1", "Ind2"), value.name = "Relationship")
rel_melt$Ind1 <- factor(rel_melt$Ind1, levels = sample_order)
rel_melt$Ind2 <- factor(rel_melt$Ind2, levels = rev(sample_order))

# Factor levels for relationship including Self
all_rel_levels <- c("Self", "Duplicate/MZ", "First degree", "Second degree",
                     "Third degree", "Unrelated")
rel_melt$Relationship <- factor(rel_melt$Relationship, levels = all_rel_levels)

# ── Define color palettes matching reference figure ─────────────────────────

# Panel A: Yellow-green-purple continuous scale (like viridis but warmer)
# The reference uses a yellow(high) -> green(mid) -> dark purple(low) scheme

# Panel B: Discrete categorical colors matching reference
rel_colors <- c(
  "Self"          = "#FFFFFF",
  "Duplicate/MZ"  = "#3B0F70",    # dark purple
  "First degree"  = "#8C6BB1",    # medium purple
  "Second degree" = "#DED5C4",    # beige/cream
  "Third degree"  = "#B5CF6B",    # yellow-green
  "Unrelated"     = "#2D3436"     # dark slate/charcoal
)

# Panel C: Bar colors matching reference
bar_colors <- c(
  "Duplicate/MZ"  = "#3B0F70",    # dark purple
  "First degree"  = "#8C6BB1",    # medium purple
  "Second degree" = "#DED5C4",    # beige/cream
  "Third degree"  = "#B5CF6B",    # yellow-green
  "Unrelated"     = "#2D3436"     # dark olive/charcoal
)

# ── Decide whether to show axis labels ──────────────────────────────────────
show_labels <- (opt$max_label == 0) | (n_samples <= opt$max_label)

if (!show_labels) {
  cat("[INFO] Too many samples (", n_samples, ") for axis labels.",
      "Set --max_label 0 to force labels.\n")
}

# ── Panel A: Continuous kinship heatmap ─────────────────────────────────────
cat("[INFO] Building Panel A (kinship heatmap)...\n")

pA <- ggplot(theta_melt, aes(x = Ind1, y = Ind2, fill = Theta)) +
  geom_tile(color = NA) +
  scale_fill_gradientn(
    colours = c("#0D0887", "#3B0F70", "#7E03A8", "#B63679",
                "#F0F921", "#FCFFA4"),
    # Use a custom viridis-like scale going dark purple -> green -> yellow
    values  = scales::rescale(c(-1, -0.5, -0.1, 0.0, 0.3, 0.5)),
    name    = "Kinship\n(theta)",
    na.value = "grey90",
    limits  = c(min(theta_melt$Theta, na.rm = TRUE),
                max(theta_melt$Theta, na.rm = TRUE))
  ) +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = if (show_labels) element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4)
                  else element_blank(),
    axis.text.y = if (show_labels) element_text(size = 4) else element_blank(),
    axis.ticks  = element_blank(),
    panel.grid  = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    plot.tag = element_text(face = "bold", size = 16),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  coord_fixed()

# ── Panel B: Categorical relationship heatmap ───────────────────────────────
cat("[INFO] Building Panel B (relationship heatmap)...\n")

pB <- ggplot(rel_melt, aes(x = Ind1, y = Ind2, fill = Relationship)) +
  geom_tile(color = NA) +
  scale_fill_manual(
    values = rel_colors,
    name   = "Relationships",
    na.value = "grey90",
    breaks = c("First degree", "Second degree", "Third degree", "Unrelated"),
    drop = FALSE
  ) +
  labs(x = NULL, y = NULL, tag = "B") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = if (show_labels) element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4)
                  else element_blank(),
    axis.text.y = if (show_labels) element_text(size = 4) else element_blank(),
    axis.ticks  = element_blank(),
    panel.grid  = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.5, "cm"),
    legend.key.width  = unit(0.5, "cm"),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    plot.tag = element_text(face = "bold", size = 16),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  coord_fixed()

# ── Panel C: Bar plot of relationship counts ────────────────────────────────
cat("[INFO] Building Panel C (relationship bar plot)...\n")

# Count pairs (exclude self comparisons)
pair_counts <- df %>%
  filter(relationship != "Self") %>%
  count(relationship, .drop = FALSE) %>%
  # Ensure all levels present
  complete(relationship = factor(rel_levels, levels = rel_levels), fill = list(n = 0))

# Remove Duplicate/MZ if 0 pairs (common in typical datasets)
if (pair_counts$n[pair_counts$relationship == "Duplicate/MZ"] == 0) {
  pair_counts <- pair_counts %>% filter(relationship != "Duplicate/MZ")
  bar_colors_use <- bar_colors[names(bar_colors) != "Duplicate/MZ"]
} else {
  bar_colors_use <- bar_colors
}

pC <- ggplot(pair_counts, aes(x = relationship, y = n, fill = relationship)) +
  geom_col(width = 0.7, color = NA) +
  geom_text(aes(label = n), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = bar_colors_use) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Pairs", tag = "C") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y  = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position = "none",
    plot.tag = element_text(face = "bold", size = 16),
    plot.margin = margin(5, 15, 5, 5)
  )

# ── Combine panels ──────────────────────────────────────────────────────────
cat("[INFO] Combining panels...\n")

# Layout: A and B same width, C narrower
combined <- pA + pB + pC +
  plot_layout(widths = c(1, 1, 0.6))

# ── Save outputs ────────────────────────────────────────────────────────────
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

pdf_out <- file.path(opt$outdir, paste0(opt$prefix, "_3panel.pdf"))
png_out <- file.path(opt$outdir, paste0(opt$prefix, "_3panel.png"))

cat("[INFO] Saving PDF:", pdf_out, "\n")
ggsave(pdf_out, combined, width = opt$width, height = opt$height, dpi = 300)

cat("[INFO] Saving PNG:", png_out, "\n")
ggsave(png_out, combined, width = opt$width, height = opt$height, dpi = 300)

# ── Also save summary statistics ────────────────────────────────────────────
stats_out <- file.path(opt$outdir, paste0(opt$prefix, "_summary_stats.txt"))
sink(stats_out)
cat("=== Relatedness Summary Statistics ===\n\n")
cat("Total samples:", n_samples, "\n")
cat("Total pairwise comparisons:", nrow(df), "\n\n")

cat("--- Theta (kinship coefficient) distribution ---\n")
cat("Min:    ", min(df$theta, na.rm = TRUE), "\n")
cat("Q1:     ", quantile(df$theta, 0.25, na.rm = TRUE), "\n")
cat("Median: ", median(df$theta, na.rm = TRUE), "\n")
cat("Mean:   ", mean(df$theta, na.rm = TRUE), "\n")
cat("Q3:     ", quantile(df$theta, 0.75, na.rm = TRUE), "\n")
cat("Max:    ", max(df$theta, na.rm = TRUE), "\n\n")

cat("--- Relationship classification ---\n")
rel_tab <- table(df$relationship)
for (r in names(rel_tab)) {
  cat(sprintf("  %-15s %6d pairs (%5.1f%%)\n", r, rel_tab[r],
              100 * rel_tab[r] / sum(rel_tab)))
}
cat("\n")

cat("--- Samples involved in close relationships (theta >= 0.177) ---\n")
close_rels <- df %>% filter(theta >= 0.177)
if (nrow(close_rels) > 0) {
  close_samples <- sort(unique(c(close_rels$ind1, close_rels$ind2)))
  cat("Number of pairs:", nrow(close_rels), "\n")
  cat("Number of unique samples involved:", length(close_samples), "\n")
  cat("Samples:", paste(close_samples, collapse = ", "), "\n\n")
  cat("--- Close relative pairs ---\n")
  close_rels_sorted <- close_rels %>% arrange(desc(theta))
  for (i in seq_len(nrow(close_rels_sorted))) {
    cat(sprintf("  %s -- %s : theta = %.4f (%s)\n",
                close_rels_sorted$ind1[i], close_rels_sorted$ind2[i],
                close_rels_sorted$theta[i], close_rels_sorted$relationship[i]))
  }
} else {
  cat("No pairs with theta >= 0.177 found.\n")
}
sink()

cat("[INFO] Summary stats saved:", stats_out, "\n")

# ── Additional plot: theta distribution histogram ───────────────────────────
hist_out <- file.path(opt$outdir, paste0(opt$prefix, "_theta_histogram.pdf"))

p_hist <- ggplot(df, aes(x = theta)) +
  geom_histogram(bins = 100, fill = "#3B0F70", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = c(0.0442, 0.0884, 0.177, 0.354),
             linetype = "dashed", color = c("#B5CF6B", "#DED5C4", "#8C6BB1", "#3B0F70"),
             linewidth = 0.7) +
  annotate("text", x = 0.354, y = Inf, label = "Dup/MZ", vjust = 1.5, hjust = -0.1, size = 3) +
  annotate("text", x = 0.177, y = Inf, label = "1st°", vjust = 1.5, hjust = -0.1, size = 3) +
  annotate("text", x = 0.0884, y = Inf, label = "2nd°", vjust = 1.5, hjust = -0.1, size = 3) +
  annotate("text", x = 0.0442, y = Inf, label = "3rd°", vjust = 1.5, hjust = -0.1, size = 3) +
  labs(x = expression(paste("Kinship coefficient (", theta, ")")),
       y = "Number of pairs",
       title = "Distribution of pairwise kinship coefficients") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(hist_out, p_hist, width = 8, height = 5, dpi = 300)
cat("[INFO] Theta histogram saved:", hist_out, "\n")

cat("[DONE] All outputs saved to:", opt$outdir, "\n")
