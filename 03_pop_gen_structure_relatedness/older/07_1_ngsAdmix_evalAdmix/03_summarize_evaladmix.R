###
# 03_summarize_evaladmix.R does this:
# Reads evalAdmix residual-like output files, computes a single residual-strength
# summary per run, writes a summary table, and generates plots of residual
# strength across K and seeds.
###
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- EDIT THESE ----
setwd("C:/Users/quent/Desktop/MS_Inversions")  # local where you copied eval_* folder
TAG <- "thin500"
EVALDIR <- file.path("05_pcangsd_byLG_ALL_2026-02-10", "05_pcangsd_global", paste0("eval_", TAG))
# If you instead copied from scratch: set EVALDIR to that folder path.
# -------------------

# evalAdmix outputs vary by version; common files include:
#   *.residuals.gz or *.residuals
#   *.corr or *.cov-like
# We'll search for anything that looks like residual matrices or correlations.
files <- list.files(EVALDIR, recursive=TRUE, full.names=TRUE)

# Helper: parse K and seed from basename like thin500_K04_seed2.*
parse_K_seed <- function(x){
  b <- basename(x)
  K <- as.integer(sub(".*_K([0-9]{2}).*", "\\1", b))
  seed <- as.integer(sub(".*_seed([0-9]+).*", "\\1", b))
  if (is.na(K)) K <- NA_integer_
  if (is.na(seed)) seed <- NA_integer_
  list(K=K, seed=seed)
}

# Try to find a "residual correlation" output file per run.
# Many evalAdmix versions write: <out>.res or <out>.residuals.gz
cand <- files[grepl("resid|residual", basename(files), ignore.case=TRUE)]
if (length(cand) == 0) stop("No residual-like files found in: ", EVALDIR)

# Function to compute "residual strength" from an NxN matrix:
# Use mean absolute off-diagonal correlation (robust single-number summary)
res_strength <- function(mat){
  mat[is.na(mat)] <- 0
  diag(mat) <- NA_real_
  mean(abs(mat), na.rm=TRUE)
}

read_any_matrix <- function(f){
  # supports plain, .gz
  if (grepl("\\.gz$", f)) {
    con <- gzfile(f, "rt")
    on.exit(close(con))
    as.matrix(fread(con, header=FALSE))
  } else {
    as.matrix(fread(f, header=FALSE))
  }
}

rows <- list()
for (f in cand) {
  meta <- parse_K_seed(f)
  if (is.na(meta$K) || is.na(meta$seed)) next
  # Try reading as matrix; if fails, skip
  mat <- try(read_any_matrix(f), silent=TRUE)
  if (inherits(mat, "try-error")) next
  if (nrow(mat) < 10 || ncol(mat) < 10) next
  rows[[length(rows)+1]] <- data.table(
    file=f, K=meta$K, seed=meta$seed,
    strength=res_strength(mat),
    n=nrow(mat)
  )
}

dt <- rbindlist(rows, fill=TRUE)
if (nrow(dt)==0) stop("Could not read any residual matrices. Check evalAdmix outputs; list files in EVALDIR.")

# Save summary
out_tsv <- file.path(EVALDIR, "summary_residuals.tsv")
fwrite(dt[order(K, seed)], out_tsv, sep="\t")

# Plot: mean ± sd across seeds
sumK <- dt[, .(mean_strength=mean(strength), sd_strength=sd(strength)), by=.(K)]
p1 <- ggplot(sumK, aes(x=K, y=mean_strength)) +
  geom_line() +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_strength-sd_strength, ymax=mean_strength+sd_strength), width=0.2) +
  theme_bw(base_size=12) +
  labs(
    title=paste0("evalAdmix residual strength vs K (", TAG, ")"),
    subtitle="Mean(|off-diagonal residual correlation|) ± SD across seeds",
    x="K", y="Residual strength"
  )
ggsave(file.path(EVALDIR, "plot_residual_strength_byK.png"), p1, width=7, height=4, dpi=300)

# Plot: each seed
p2 <- ggplot(dt, aes(x=K, y=strength, color=factor(seed))) +
  geom_line() + geom_point(size=1.8) +
  theme_bw(base_size=12) +
  labs(
    title=paste0("evalAdmix residual strength vs K by seed (", TAG, ")"),
    x="K", y="Residual strength", color="seed"
  )
ggsave(file.path(EVALDIR, "plot_residual_strength_byK_seed.png"), p2, width=7, height=4, dpi=300)

message("Wrote: ", out_tsv)
message("Done.")
