#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: plot_heterozygosity_boxplot.R heterozygosity_per_sample.tsv out.pdf [group.tsv(optional)]")
}

het_file <- args[1]
out_pdf  <- args[2]
group_file <- ifelse(length(args) >= 3, args[3], NA)

d <- read.table(het_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Optional grouping file: two columns -> sample \t group
if (!is.na(group_file)) {
  g <- read.table(group_file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  colnames(g) <- c("sample","group")
  d <- merge(d, g, by="sample", all.x=TRUE)
} else {
  d$group <- "all"
}

pdf(out_pdf, width=7, height=4)
boxplot(het_saf1 ~ group, data=d, las=2, ylab="Heterozygosity (ANGSD single-sample SFS)", xlab="")
stripchart(het_saf1 ~ group, data=d, vertical=TRUE, method="jitter", pch=16, add=TRUE)
dev.off()
