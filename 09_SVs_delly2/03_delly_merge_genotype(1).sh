#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 127
#SBATCH --mem=237G
#SBATCH -t 2-00:00:00
#SBATCH -J delly_geno
#SBATCH -o logs/03_merge_geno.%j.out
#SBATCH -e logs/03_merge_geno.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_delly_config.sh"

dv_init_dirs
dv_log "=== STEPS 3–8: Merge → Regenotype → Cohort merge → Germline filter ==="

dv_check_file "$REF" "Reference"
dv_check_file "$EXCL_BED" "Exclusion BED"
dv_check_file "$SAMPLES_ALL" "Sample list (all)"
dv_check_file "$SAMPLES_UNRELATED" "Sample list (unrelated 81)"
[[ -x "${DELLY_BIN}" ]] || dv_die "DELLY binary not executable: ${DELLY_BIN}"

N_ALL=$(wc -l < "$SAMPLES_ALL")
N_UNREL=$(wc -l < "$SAMPLES_UNRELATED")
dv_log "Full cohort: ${N_ALL} | Unrelated subset: ${N_UNREL}"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3. delly merge — shared DEL site list
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- STEP 3: delly merge (shared DEL site list) ---"

DISC_LIST="${DIR_SITES}/disc_bcf_list.txt"
ls "${DIR_DISC}"/*.disc.DEL.bcf > "${DISC_LIST}"
N_DISC=$(wc -l < "${DISC_LIST}")
dv_log "  Input: ${N_DISC} per-sample DEL discovery BCFs"

# Ensure all indexed
while IFS= read -r bcf; do
  [[ -f "${bcf}.csi" ]] || bcftools index "$bcf"
done < "${DISC_LIST}"

SITES_BCF="${DIR_SITES}/sites.DEL.bcf"

"${DELLY_BIN}" merge \
  -o "${SITES_BCF}" \
  $(cat "${DISC_LIST}")

bcftools index "${SITES_BCF}"

N_SITES=$(bcftools view -H "${SITES_BCF}" | wc -l)
dv_log "  Merged DEL site list: ${N_SITES} candidate sites"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4. Regenotype merged sites in all 226 samples
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- STEP 4: Regenotype ${N_ALL} samples ---"

run_genotype() {
  local sid="$1"
  local bam="${BAMDIR}/${sid}${BAM_SUFFIX}"
  local out="${DIR_GENO}/${sid}.geno.DEL.bcf"
  local logf="${DIR_LOGS}/geno_${sid}.log"

  if [[ -f "${out}" && -f "${out}.csi" ]]; then
    echo "[$(date '+%T')] ${sid}: already genotyped, skipping"
    return 0
  fi

  echo "[$(date '+%T')] ${sid}: regenotyping..."

  # -v sites.bcf for regenotyping, -t DEL, -h threads
  "${DELLY_BIN}" call \
    -t DEL \
    -g "${REF}" \
    -x "${EXCL_BED}" \
    -v "${SITES_BCF}" \
    -h "${DELLY_THREADS_PER_CALL}" \
    -o "${out}" \
    "${bam}" \
    2>>"${logf}"

  bcftools index "${out}" 2>>"${logf}"
  echo "[$(date '+%T')] ${sid}: done"
}

export -f run_genotype
export BAMDIR BAM_SUFFIX REF EXCL_BED DELLY_BIN DELLY_THREADS_PER_CALL SITES_BCF DIR_GENO DIR_LOGS

parallel -j "${DELLY_PARALLEL}" --joblog "${DIR_LOGS}/genotype_parallel.log" \
  run_genotype {} \
  :::: "${SAMPLES_ALL}"

# Verify
FAIL=0
while IFS= read -r sid; do
  if [[ ! -f "${DIR_GENO}/${sid}.geno.DEL.bcf" ]]; then
    dv_err "Missing genotyped BCF: ${sid}"
    ((FAIL++)) || true
  fi
done < "$SAMPLES_ALL"
[[ ${FAIL} -eq 0 ]] || dv_die "${FAIL} failed regenotyping."
dv_log "  Regenotyping complete (${N_ALL} samples)"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5. bcftools merge — 226-sample cohort BCF
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- STEP 5: bcftools merge (226 cohort) ---"

GENO_LIST="${DIR_MERGED}/geno_bcf_list.txt"
ls "${DIR_GENO}"/*.geno.DEL.bcf > "${GENO_LIST}"

COHORT_226="${DIR_MERGED}/cohort_226.DEL.merged.bcf"

bcftools merge \
  -m id \
  -O b \
  -o "${COHORT_226}" \
  --threads "${THREADS}" \
  $(cat "${GENO_LIST}")

bcftools index "${COHORT_226}"

N_226=$(bcftools view -H "${COHORT_226}" | wc -l)
dv_log "  226 cohort: ${N_226} DEL sites"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6. Subset to 81 unrelated
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- STEP 6: Subset to ${N_UNREL} unrelated ---"

COHORT_81="${DIR_SUBSET81}/cohort_81.DEL.merged.bcf"

bcftools view \
  -S "${SAMPLES_UNRELATED}" \
  --force-samples \
  -O b \
  -o "${COHORT_81}" \
  "${COHORT_226}"

bcftools index "${COHORT_81}"

# Remove monomorphic sites in the 81 subset
COHORT_81_TRIM="${DIR_SUBSET81}/cohort_81.DEL.merged.trimmed.bcf"
bcftools view \
  -i 'COUNT(GT="alt")>0' \
  -O b \
  -o "${COHORT_81_TRIM}" \
  "${COHORT_81}"

bcftools index "${COHORT_81_TRIM}"

N_81=$(bcftools view -H "${COHORT_81_TRIM}" | wc -l)
dv_log "  81 subset: ${N_81} DEL sites (polymorphic in 81)"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 7. delly filter -f germline on 81 unrelated
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- STEP 7: delly filter -f germline (81 unrelated) ---"

GERMLINE_81="${DIR_FILTERED}/cohort_81.DEL.germline.bcf"

if "${DELLY_BIN}" filter \
  -f germline \
  -o "${GERMLINE_81}" \
  "${COHORT_81_TRIM}" 2>"${DIR_LOGS}/delly_filter_germline.log"; then

  bcftools index "${GERMLINE_81}"
  N_GERM=$(bcftools view -H "${GERMLINE_81}" | wc -l)
  dv_log "  Germline filter PASS: ${N_GERM} DEL sites"
else
  dv_log "  WARNING: delly filter -f germline returned non-zero (common with family structure)"
  dv_log "  Falling back to merged+trimmed BCF."
fi

# ─────────────────────────────────────────────────────────────────────────────
# STEP 8. Final catalogs + GT matrices
# ─────────────────────────────────────────────────────────────────────────────
dv_log "--- STEP 8: Final DEL catalogs + genotype matrices ---"

# --- 226 catalog ---
FINAL_226_VCF="${DIR_FINAL}/catalog_226.DEL.vcf.gz"
bcftools view \
  -i 'INFO/SVTYPE="DEL"' \
  -O z -o "${FINAL_226_VCF}" \
  "${COHORT_226}"
tabix -p vcf "${FINAL_226_VCF}"

# 226 GT matrix
MATRIX_226="${DIR_FINAL}/catalog_226.DEL.GT_matrix.tsv"
{
  echo -e "CHROM\tPOS\tEND\tID\tSVLEN\t$(bcftools query -l "${FINAL_226_VCF}" | tr '\n' '\t' | sed 's/\t$//')"
  bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN[\t%GT]\n' "${FINAL_226_VCF}"
} > "${MATRIX_226}"

N_226_FINAL=$(tail -n +2 "${MATRIX_226}" | wc -l)
dv_log "  226 catalog: ${N_226_FINAL} DELs"

# --- 81 catalog ---
if [[ -f "${GERMLINE_81}" ]]; then
  INPUT_81="${GERMLINE_81}"
  LABEL_81="germline"
else
  INPUT_81="${COHORT_81_TRIM}"
  LABEL_81="merged_trimmed"
fi

FINAL_81_VCF="${DIR_FINAL}/catalog_81.DEL.${LABEL_81}.vcf.gz"
bcftools view -i 'INFO/SVTYPE="DEL"' -O z -o "${FINAL_81_VCF}" "${INPUT_81}"
tabix -p vcf "${FINAL_81_VCF}"

# 81 PASS
FINAL_81_PASS="${DIR_FINAL}/catalog_81.DEL.${LABEL_81}.PASS.vcf.gz"
bcftools view -f PASS -i 'INFO/SVTYPE="DEL"' -O z -o "${FINAL_81_PASS}" "${INPUT_81}"
tabix -p vcf "${FINAL_81_PASS}"

# 81 GT matrix (PASS)
MATRIX_81="${DIR_FINAL}/catalog_81.DEL.${LABEL_81}.PASS.GT_matrix.tsv"
{
  echo -e "CHROM\tPOS\tEND\tID\tSVLEN\t$(bcftools query -l "${FINAL_81_PASS}" | tr '\n' '\t' | sed 's/\t$//')"
  bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN[\t%GT]\n' "${FINAL_81_PASS}"
} > "${MATRIX_81}"

N_81_PASS=$(bcftools view -H "${FINAL_81_PASS}" 2>/dev/null | wc -l)
dv_log "  81 catalog (${LABEL_81}): ${N_81} total, ${N_81_PASS} PASS"

# --- DEL BEDs for annotation ---
for vcf_path in "${FINAL_226_VCF}" "${FINAL_81_PASS}"; do
  bed_out="${vcf_path%.vcf.gz}.bed"
  bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/SVLEN\n' "${vcf_path}" \
    | awk 'BEGIN{OFS="\t"} { s=$2; e=$3; if(s>e){t=s;s=e;e=t} print $1,s,e,$4,$5 }' \
    > "${bed_out}"
  dv_log "  BED: $(basename "${bed_out}") ($(wc -l < "${bed_out}") sites)"
done

dv_log "=== STEPS 3–8 COMPLETE ==="
dv_log "Next: submit 04_annotation_layers.sh"
