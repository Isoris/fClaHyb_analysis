#!/usr/bin/env bash
# =============================================================================
# 00_delly_config.sh — Central configuration for DELLY DEL pipeline (v3)
# =============================================================================
#
# CHANGES FROM v2:
#   - Manifest-driven BAM sourcing (no flat BAMDIR glob)
#   - Duplicate-marked BAM output directory
#   - Unconditional 50-kb chromosome-end masking variable
#   - All variables exported (use set -a / set +a when sourcing)
#   - Removed old BAMDIR/BAM_SUFFIX glob-based logic
#
# HOW TO SOURCE THIS (in SLURM scripts):
#   SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
#   CONFIG="${SCRIPT_DIR}/00_delly_config.sh"
#   [[ -f "${CONFIG}" ]] || { echo "Missing config" >&2; exit 1; }
#   set -a
#   source "${CONFIG}"
#   set +a
# =============================================================================

# ── Module loads (required for compiled DELLY) ──────────────────────────────
module load HTSlib/1.17-cpeGNU-23.03
module load Boost/1.81.0-cpeGNU-23.03

# ── DELLY binary (compiled from source, v1.7.3) ────────────────────────────
DELLY_BIN="/project/lt200308-agbsci/01-catfish_assembly/delly/src/delly"

# ── Project root ────────────────────────────────────────────────────────────
DELLY_PROJECT="${DELLY_PROJECT:-/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04}"

# ── Reference genome ────────────────────────────────────────────────────────
REF="${DELLY_PROJECT}/00-samples/fClaHyb_Gar_LG.fa"
REF_FAI="${REF}.fai"

# ── Exclusion BED (built by 01_prep from callable + chromosome-end masks) ──
EXCL_BED="${DELLY_PROJECT}/delly_sv/exclude.minimal.bed"

# ── BAM manifest (manifest-driven, NOT flat directory glob) ─────────────────
# TSV format: sample_id <tab> unfiltered_bam <tab> filtered_bam
# Column 1 = sample ID
# Column 2 = unfiltered/raw BAM path  → USE THIS for DELLY (after markdup)
# Column 3 = P99/TLEN/MAPQ30 filtered BAM → IGNORE for DELLY
BAM_MANIFEST="${DELLY_PROJECT}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"

# ── Duplicate-marked BAM output directory ───────────────────────────────────
# DELLY requires duplicate-marked BAMs. Since manifest column 2 BAMs are raw
# minimap2 alignments (no markdup), we create markdup copies here.
# The 01_prep step runs samtools markdup on each column-2 BAM and writes to:
DIR_MARKDUP="${DELLY_PROJECT}/delly_sv/00_markdup"

# ── Sample lists ────────────────────────────────────────────────────────────
SAMPLES_ALL="${DELLY_PROJECT}/delly_sv/samples_all_226.txt"
SAMPLES_UNRELATED="${DELLY_PROJECT}/delly_sv/samples_unrelated_81.txt"
NATORA_KEEP="${DELLY_PROJECT}/popstruct_thin/05_ngsrelate/catfish_first_degree_pairwise_toKeep.txt"

# ── PA-Roary callable data (for building exclude BED) ──────────────────────
PA_CALLABLE_BP="${DELLY_PROJECT}/pa_roary_results/02_coverage_matrix_callable/callable_bp_per_bin.tsv"

# ── Chromosome-end masking (unconditional) ──────────────────────────────────
# Always exclude the first and last CHR_END_MASK_BP of every chromosome.
# This catches telomeric/centromeric junk regardless of callable data.
CHR_END_MASK_BP=50000

# ── Callable-based exclusion thresholds ─────────────────────────────────────
# 50-kb bins with callable_bp below this → uncallable
EXCL_MIN_CALLABLE_BP=500
# Merged uncallable blocks must be at least this big
EXCL_MIN_BLOCK_BP=50000

# ── Annotation files ───────────────────────────────────────────────────────
GFF3="${DELLY_PROJECT}/00-samples/fClaHyb_Gar_LG.from_CGAR.gff3_polished.sorted.gff3"
ANNOT_DIR="${DELLY_PROJECT}/pa_roary_results/annotation_flag_beds"
GENE_BED="${ANNOT_DIR}/features/GENE.bed"
EXON_BED="${ANNOT_DIR}/features/EXON.bed"
CDS_BED="${ANNOT_DIR}/features/CDS.bed"
REPEAT_BED="${DELLY_PROJECT}/fClaHyb_Gar_LG.mask_regions.softacgt.renamed.bed"
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
DIR_PLOTS="${OUTDIR}/12_plots"
DIR_LOGS="${OUTDIR}/logs"

# ── Performance ─────────────────────────────────────────────────────────────
THREADS="${SLURM_CPUS_PER_TASK:-127}"
DELLY_THREADS_PER_CALL=4
DELLY_PARALLEL=30

# ── mosdepth (depth support layer) ─────────────────────────────────────────
DEPTH_WINDOW=500
DEPTH_MAPQ=30
DEPTH_THREADS=4

# ── Mate distance thresholds ───────────────────────────────────────────────
MATE_WARN_KB=20
MATE_SUSPICIOUS_KB=50
MATE_EXTREME_KB=100

# ── Strict final catalog filter for marker-grade DEL calls ─────────────────
STRICT_DEL_REQUIRE_PASS=1
STRICT_DEL_REQUIRE_PRECISE=1
STRICT_DEL_MIN_QUAL=500
STRICT_DEL_MIN_PE=5   # PE > 4

# ── Ancestry / admixture (for plotting) ────────────────────────────────────
PA_BEST_SEED_TABLE="${DELLY_PROJECT}/popstruct_thin/05_ngsadmix_global/best_seed_by_K.tsv"
PA_NGSADMIX_DIR="${DELLY_PROJECT}/popstruct_thin/05_ngsadmix_global/runs_thin500"

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
    "$OUTDIR" "$DIR_MARKDUP" \
    "$DIR_DISC" "$DIR_SITES" "$DIR_GENO" \
    "$DIR_MERGED" "$DIR_SUBSET81" "$DIR_FILTERED" "$DIR_FINAL" \
    "$DIR_ANNOT" "$DIR_DEPTH" "$DIR_MATDIST" "$DIR_SUMMARY" \
    "$DIR_PLOTS" "$DIR_LOGS"
}

# ── Manifest reader ────────────────────────────────────────────────────────
# Reads the BAM manifest TSV and populates parallel arrays.
# After calling this, the following arrays are available:
#   MANIFEST_SIDS[]    — sample IDs
#   MANIFEST_RAW_BAMS[]  — unfiltered BAM paths (column 2)
# And the count:
#   MANIFEST_N_SAMPLES

dv_read_manifest() {
  [[ -s "$BAM_MANIFEST" ]] || dv_die "BAM manifest not found or empty: $BAM_MANIFEST"

  MANIFEST_SIDS=()
  MANIFEST_RAW_BAMS=()

  {
    read -r _ _ _
    while IFS=$'\t' read -r sid raw_bam _ignore; do
      [[ -z "$sid" ]] && continue
      MANIFEST_SIDS+=("$sid")
      MANIFEST_RAW_BAMS+=("$raw_bam")
    done
  } < "$BAM_MANIFEST"

  MANIFEST_N_SAMPLES=${#MANIFEST_SIDS[@]}
  [[ $MANIFEST_N_SAMPLES -gt 0 ]] || dv_die "Zero samples loaded from manifest: $BAM_MANIFEST"

  dv_log "Manifest loaded: ${MANIFEST_N_SAMPLES} samples from ${BAM_MANIFEST}"
}

# Get the markdup BAM path for a given sample ID
dv_markdup_bam() {
  local sid="$1"
  echo "${DIR_MARKDUP}/${sid}.markdup.bam"
}
