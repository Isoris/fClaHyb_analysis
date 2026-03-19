#!/usr/bin/env bash
# =============================================================================
# 00_delly_config.sh — Central configuration for DELLY DEL pipeline
# =============================================================================
# Source this from any pipeline script:
#   source "$(dirname "$0")/00_delly_config.sh"
# =============================================================================

# ── Module loads (required for DELLY) ───────────────────────────────────────
module load HTSlib/1.17-cpeGNU-23.03
module load Boost/1.81.0-cpeGNU-23.03

# ── DELLY binary (compiled from source, v1.7.3) ────────────────────────────
DELLY_BIN="/project/lt200308-agbsci/01-catfish_assembly/delly/src/delly"

# ── Project root ────────────────────────────────────────────────────────────
DELLY_PROJECT="${DELLY_PROJECT:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

# ── Reference genome ────────────────────────────────────────────────────────
REF="${DELLY_PROJECT}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${REF}.fai"

# ── Exclusion BED (built from callable data by 01_prep) ─────────────────────
EXCL_BED="${DELLY_PROJECT}/delly_sv/exclude.minimal.bed"

# ── BAM directory ───────────────────────────────────────────────────────────
# >>> EDIT to match your actual directory and suffix <<<
BAMDIR="${DELLY_PROJECT}/02-bwa2"
BAM_SUFFIX=".markdup.bam"

# ── Sample lists ────────────────────────────────────────────────────────────
SAMPLES_ALL="${DELLY_PROJECT}/delly_sv/samples_all_226.txt"
SAMPLES_UNRELATED="${DELLY_PROJECT}/delly_sv/samples_unrelated_81.txt"
NATORA_KEEP="${DELLY_PROJECT}/popstruct_thin/05_ngsrelate/catfish_first_degree_pairwise_toKeep.txt"

# ── PA-Roary callable data (for building exclude BED) ──────────────────────
PA_CALLABLE_BP="${DELLY_PROJECT}/pa_roary_results/02_coverage_matrix_callable/callable_bp_per_bin.tsv"
PA_WINDOW_ABSENT="${DELLY_PROJECT}/pa_roary_results/02_coverage_matrix_callable/window_absent_counts.tsv"

# ── Annotation files ───────────────────────────────────────────────────────
GFF3="${DELLY_PROJECT}/00-samples/fClaHyb_Gar_LG.from_CGAR.gff3_polished.sorted.gff3"
ANNOT_DIR="${DELLY_PROJECT}/pa_roary_results/annotation_flag_beds"
GENE_BED="${ANNOT_DIR}/features/GENE.bed"
EXON_BED="${ANNOT_DIR}/features/EXON.bed"
CDS_BED="${ANNOT_DIR}/features/CDS.bed"
REPEAT_BED="${DELLY_PROJECT}/fClaHyb_Gar_LG.mask_regions.softacgt.bed"
SIMPLE_REPEAT_BED=""
SEGDUP_BED=""

# ── Output directories ─────────────────────────────────────────────────────
OUTDIR="${DELLY_PROJECT}/delly_sv"
DIR_DISC="${OUTDIR}/01_discovery"
DIR_SITES="${OUTDIR}/02_merged_sites"
DIR_GENO="${OUTDIR}/03_genotyped"
DIR_MERGED="${OUTDIR}/04_merged_cohort"
DIR_SUBSET81="${OUTDIR}/05_subset_81"
DIR_FILTERED="${OUTDIR}/06_germline_filtered"
DIR_FINAL="${OUTDIR}/07_final_catalogs"
DIR_ANNOT="${OUTDIR}/08_annotation"
DIR_DEPTH="${OUTDIR}/09_depth_support"
DIR_MATDIST="${OUTDIR}/10_mate_distance_qc"
DIR_SUMMARY="${OUTDIR}/11_summary"
DIR_LOGS="${OUTDIR}/logs"

# ── Performance ─────────────────────────────────────────────────────────────
THREADS="${SLURM_CPUS_PER_TASK:-127}"

# DELLY supports -h <threads> per call (default 4)
DELLY_THREADS_PER_CALL=4
# Parallel DELLY calls: 127 cores / 4 threads = ~31, use 30 for safety
DELLY_PARALLEL=30

# ── mosdepth (depth support layer) ─────────────────────────────────────────
DEPTH_WINDOW=500
DEPTH_MAPQ=30
DEPTH_THREADS=4

# ── Mate distance thresholds ───────────────────────────────────────────────
MATE_WARN_KB=20
MATE_SUSPICIOUS_KB=50
MATE_EXTREME_KB=100

# ── Exclude BED thresholds (from callable_bp_per_bin.tsv) ──────────────────
# 50-kb bins with callable_bp below this → uncallable
EXCL_MIN_CALLABLE_BP=500
# Merged uncallable blocks must be at least this big for exclude BED
EXCL_MIN_BLOCK_BP=50000

# =============================================================================
# Helper functions
# =============================================================================
dv_timestamp() { date '+%F %T'; }
dv_log() { echo "[$(dv_timestamp)] [DELLY-DEL] $*"; }
dv_err() { echo "[$(dv_timestamp)] [DELLY-DEL] [ERROR] $*" >&2; }
dv_die() { dv_err "$@"; exit 1; }

dv_check_file() {
  local f="$1" label="${2:-file}"
  [[ -e "$f" ]] || dv_die "Missing ${label}: ${f}"
}

dv_check_cmd() {
  command -v "$1" &>/dev/null || dv_die "Command not found: $1"
}

dv_init_dirs() {
  mkdir -p \
    "$OUTDIR" \
    "$DIR_DISC" "$DIR_SITES" "$DIR_GENO" \
    "$DIR_MERGED" "$DIR_SUBSET81" "$DIR_FILTERED" "$DIR_FINAL" \
    "$DIR_ANNOT" "$DIR_DEPTH" "$DIR_MATDIST" "$DIR_SUMMARY" "$DIR_LOGS"
}

dv_sample_from_bam() {
  basename "$1" "${BAM_SUFFIX}"
}
