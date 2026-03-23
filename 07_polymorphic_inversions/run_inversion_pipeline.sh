#!/usr/bin/env bash
# =============================================================================
# run_inversion_pipeline.sh — Master runner for the inversion local-PCA workflow
#
# This is NOT meant to be run as a single script. It is a reference for the
# exact commands, in order, with real paths from your project.
#
# Run each section manually or adapt into individual sbatch submissions.
# =============================================================================

set -euo pipefail

# ── Project paths ───────────────────────────────────────────────────────────
BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR="${BASE}/inversion_localpca"
HETDIR="${BASE}/het_roh"
SCRIPTS="${INVDIR}/scripts"   # put the STEP*.{py,sh,R,slurm} files here

REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
BAMLIST="${HETDIR}/01_inputs_check/bamlist_qcpass.txt"
SAMPLES_IND="${HETDIR}/01_inputs_check/samples.ind"

mkdir -p "${INVDIR}"/{chunks,global_sfs,02_snps_beagle,03_sites,04_dosage_by_chr,05_local_pca,06_mds_candidates,07_het_overlap,08_regional_pca,09_combined_plots,logs}

# =============================================================================
# PHASE 1: Build callable mask (one-time, if not already done)
# =============================================================================

# Check if the inversion callable mask already exists from a previous run
MASK_PREFIX="${BASE}/fClaHyb_Gar_LG.mask_regions"
if [[ ! -f "${MASK_PREFIX}.inversion_acgt_allcase.1-based.angsd.idx" ]]; then
    echo "=== STEP01 + STEP02: Build inversion callable mask ==="
    bash "${SCRIPTS}/STEP02_make_inversion_callable_angsd_mask.sh" \
        --fasta "$REF" \
        --mask-script "${SCRIPTS}/STEP01_mask_regions_from_fasta.py" \
        --prefix "$MASK_PREFIX"
else
    echo "[SKIP] Inversion callable mask already exists"
fi

# =============================================================================
# PHASE 2: SFS prior estimation (STEP04 → STEP05 → STEP06)
# =============================================================================

# You need a chunk_rf.list file — one path per line, each pointing to an RF file
# that defines a chunk of chromosomes/contigs.
# Example: create it from your chromosome names:
#   awk '{print $1}' ${REF}.fai | while read chr; do
#     echo "$chr" > ${INVDIR}/chunks/rf_files/${chr}.rf.txt
#     echo "${INVDIR}/chunks/rf_files/${chr}.rf.txt"
#   done > ${INVDIR}/chunk_rf.list

CHUNK_RF_LIST="${INVDIR}/chunk_rf.list"
N_CHUNKS=$(wc -l < "$CHUNK_RF_LIST")

echo "=== STEP04: Run ANGSD SAF per chunk (array job) ==="
echo "sbatch --array=0-$((N_CHUNKS - 1))%8 ${SCRIPTS}/STEP04_run_angsd_saf_chunks_inversion.slurm ${CHUNK_RF_LIST}"

echo "=== STEP05: Merge chunk SAFs ==="
echo "sbatch ${SCRIPTS}/STEP05_merge_chunk_saf_to_global_inversion.slurm"

echo "=== STEP06: Make global folded SFS + mean pest ==="
echo "sbatch ${SCRIPTS}/STEP06_make_global_folded_sfs_and_mean_pest_inversion.slurm"

# =============================================================================
# PHASE 3: SNP calling and dosage (STEP07 → STEP08)
# =============================================================================

echo "=== STEP07: Call SNPs + BEAGLE (array job) ==="
echo "sbatch --array=0-$((N_CHUNKS - 1))%8 ${SCRIPTS}/STEP07_angsd_call_snps_and_beagle_inversion.slurm ${CHUNK_RF_LIST}"

# After STEP07, you may want to concatenate per-chunk beagle files if you
# ran per-chromosome chunks. Otherwise, if each chunk is one chromosome,
# you can feed each .beagle.gz directly to STEP08.

echo "=== STEP08: Convert BEAGLE to dosage per chromosome ==="
DOSAGE_DIR="${INVDIR}/04_dosage_by_chr"
# Run for each beagle file:
for beagle in "${INVDIR}"/02_snps_beagle/*.beagle.gz; do
    echo "python3 ${SCRIPTS}/STEP08_beagle_to_dosage_by_chr.py --beagle ${beagle} --outdir ${DOSAGE_DIR}"
done

# =============================================================================
# PHASE 4: Local PCA + MDS + candidate discovery (STEP09 → STEP10)
# =============================================================================

echo "=== STEP09: Local PCA per chromosome ==="
LPCA_DIR="${INVDIR}/05_local_pca"
WINSIZE=100
NPC=2
# Run per chromosome:
for sites in "${DOSAGE_DIR}"/*.sites.tsv.gz; do
    chrom=$(basename "$sites" .sites.tsv.gz)
    dosage="${DOSAGE_DIR}/${chrom}.dosage.tsv.gz"
    echo "Rscript ${SCRIPTS}/STEP09_local_pca_windows_by_chr.R ${sites} ${dosage} ${LPCA_DIR}/STEP09_${chrom} ${NPC} ${WINSIZE}"
done

echo "=== STEP10: Window MDS + outlier detection ==="
MDS_DIR="${INVDIR}/06_mds_candidates"
echo "Rscript ${SCRIPTS}/STEP10_window_mds_outliers.R ${LPCA_DIR} ${MDS_DIR}/inversion_localpca ${NPC} 20 3 1000000 3"

# =============================================================================
# PHASE 5: Bridge het_roh theta to inversion pipeline (STEP10b)
# =============================================================================

echo "=== STEP10b: Parse .pestPG to per-sample tP windows ==="
PESTPG_DIR="${HETDIR}/02_heterozygosity/03_theta/multiscale"
BRIDGE_DIR="${INVDIR}/07_het_overlap"
echo "Rscript ${SCRIPTS}/STEP10b_parse_pestPG_to_sample_theta_windows.R \
    ${PESTPG_DIR} ${SAMPLES_IND} ${BRIDGE_DIR}/het_bridge \
    win50000.step10000"

# Also make a population-level mean tP summary for STEP11:
# (you can do this from the STEP10b output)
cat << 'RSCRIPT'
# Quick R to make population-mean tP per window:
library(data.table)
dt <- fread("${BRIDGE_DIR}/het_bridge.sample_tP_windows.tsv.gz")
pop_mean <- dt[, .(tP = mean(tP, na.rm=TRUE), nSamples = .N),
               by = .(chrom, WinStart, WinStop, WinCenter)]
fwrite(pop_mean, "${BRIDGE_DIR}/population_mean_tP_windows.tsv.gz", sep="\t")
RSCRIPT

# =============================================================================
# PHASE 6: Overlap + validation (STEP11 → STEP12 → STEP13)
# =============================================================================

echo "=== STEP11: Overlap candidate regions with population tP ==="
CAND_FILE="${MDS_DIR}/inversion_localpca.candidate_regions.tsv.gz"
POP_TP="${BRIDGE_DIR}/population_mean_tP_windows.tsv.gz"
echo "Rscript ${SCRIPTS}/STEP11_overlap_candidate_regions_with_theta_and_het.R \
    ${CAND_FILE} ${POP_TP} ${BRIDGE_DIR}/inversion_localpca \
    chrom WinStart WinStop tP"

echo "=== STEP12: Regional PCA + plot C per candidate ==="
# For each candidate_id (inspect the candidate_regions table first):
CANDIDATE_ID=1  # adjust based on STEP10 output
CHR="C_gar_LG12"  # adjust to the chromosome of the candidate
REGIONAL_DIR="${INVDIR}/08_regional_pca"
echo "Rscript ${SCRIPTS}/STEP12_candidate_region_pca_groups_and_plotC.R \
    ${CAND_FILE} \
    ${CANDIDATE_ID} \
    ${DOSAGE_DIR}/${CHR}.sites.tsv.gz \
    ${DOSAGE_DIR}/${CHR}.dosage.tsv.gz \
    ${BRIDGE_DIR}/het_bridge.sample_tP_windows.tsv.gz \
    ${REGIONAL_DIR}/inversion_localpca"

echo "=== STEP13: Combined panels A/B/C ==="
COMBINED_DIR="${INVDIR}/09_combined_plots"
echo "Rscript ${SCRIPTS}/STEP13_make_combined_panels_ABC.R \
    ${REGIONAL_DIR}/inversion_localpca \
    ${CANDIDATE_ID} \
    ${COMBINED_DIR}/inversion_localpca"

echo ""
echo "Pipeline reference complete."
echo "Run each step sequentially after confirming the previous step succeeded."
