#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Usage: plot_theta_ideogram_simple.R theta_window.tsv chrom_sizes.tsv out.pdf")
}

theta_file <- args[1]
sizes_file <- args[2]
out_pdf <- args[3]

# thetaStat do_stat output is whitespace-delimited; keep flexible
t <- read.table(theta_file, header=TRUE, stringsAsFactors=FALSE)
# Common columns in thetaStat outputs include: chr, winCenter, tW, tP, ... (varies by version)
# We'll try to infer columns:
# Expect at least: chr, start, end, and a theta-like column.
# If your file differs, edit these 3 lines.

# Heuristic: try to find columns
cn <- colnames(t)
chr_col <- which(cn %in% c("chr","Chromo","chrom","contig"))[1]
start_col <- which(cn %in% c("start","Start","beg","winStart"))[1]
end_col <- which(cn %in% c("end","End","winEnd"))[1]

# theta-like column preference order
theta_col <- which(cn %in% c("tP","thetaPi","ThetaPi","pi","tW","thetaW","ThetaW"))[1]

if (any(is.na(c(chr_col,start_col,end_col,theta_col)))) {
  stop("Could not auto-detect columns. Please open the TSV and map chr/start/end/theta manually.")
}

t$chr <- t[[chr_col]]
t$start <- t[[start_col]]
t$end <- t[[end_col]]
t$theta <- t[[theta_col]]
t$mid <- (t$start + t$end) / 2

# chromosome sizes
s <- read.table(sizes_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(s) <- c("chr","len")

# order chromosomes as in sizes file
t$chr <- factor(t$chr, levels=s$chr)

pdf(out_pdf, width=10, height=4)
par(mar=c(4,4,2,1))
plot(t$mid, t$theta, pch=16, cex=0.4,
     xlab="Genomic position (within chromosome)", ylab="Theta (windowed)",
     main="Theta ideogram-like track", xaxt="n")
# This minimal plot is “within chromosome” only; for true genome-wide ideogram,
# you would cum-sum chromosome offsets. Keep it simple for now:
dev.off()
