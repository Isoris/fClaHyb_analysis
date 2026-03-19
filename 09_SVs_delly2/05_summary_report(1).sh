#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH -t 0-02:00:00
#SBATCH -J delly_sum
#SBATCH -o logs/05_summary.%j.out
#SBATCH -e logs/05_summary.%j.err
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_delly_config.sh"

dv_init_dirs
dv_log "=== FINAL SUMMARY REPORT ==="

REPORT="${DIR_SUMMARY}/delly_DEL_summary_report.txt"
COUNTS="${DIR_SUMMARY}/delly_DEL_counts.tsv"

# ── Gather counts ───────────────────────────────────────────────────────────

# 226 catalog
FINAL_226_VCF="${DIR_FINAL}/catalog_226.DEL.vcf.gz"
N_226=0
[[ -f "${FINAL_226_VCF}" ]] && N_226=$(bcftools view -H "${FINAL_226_VCF}" | wc -l)

# 81 catalog
FINAL_81_VCF=$(ls "${DIR_FINAL}"/catalog_81.DEL.*.vcf.gz 2>/dev/null | grep -v PASS | head -1)
FINAL_81_PASS=$(ls "${DIR_FINAL}"/catalog_81.DEL.*.PASS.vcf.gz 2>/dev/null | head -1)
N_81=0; N_81_PASS=0
[[ -n "${FINAL_81_VCF}" && -f "${FINAL_81_VCF}" ]] && N_81=$(bcftools view -H "${FINAL_81_VCF}" | wc -l)
[[ -n "${FINAL_81_PASS}" && -f "${FINAL_81_PASS}" ]] && N_81_PASS=$(bcftools view -H "${FINAL_81_PASS}" | wc -l)

# Functional
FUNC="${DIR_ANNOT}/catalog_226.functional_class.tsv"
N_CDS=0; N_EXON=0; N_INTRONIC=0; N_INTERGENIC=0; N_GENES=0
if [[ -f "${FUNC}" ]]; then
  N_CDS=$(awk -F'\t' 'NR>1 && $10=="CDS_overlap"' "${FUNC}" | wc -l)
  N_EXON=$(awk -F'\t' 'NR>1 && $10=="exon_overlap"' "${FUNC}" | wc -l)
  N_INTRONIC=$(awk -F'\t' 'NR>1 && $10=="intronic"' "${FUNC}" | wc -l)
  N_INTERGENIC=$(awk -F'\t' 'NR>1 && $10=="intergenic"' "${FUNC}" | wc -l)
  N_GENES=$(awk -F'\t' 'NR>1 && $10!="intergenic" {print $11}' "${FUNC}" \
    | tr ',' '\n' | grep -v '^\.$' | sort -u | wc -l)
fi

# Repeats
N_REP=0; N_NOREP=0
[[ -f "${DIR_ANNOT}/catalog_226.DELs_in_repeats.bed" ]] && N_REP=$(wc -l < "${DIR_ANNOT}/catalog_226.DELs_in_repeats.bed")
[[ -f "${DIR_ANNOT}/catalog_226.DELs_not_in_repeats.bed" ]] && N_NOREP=$(wc -l < "${DIR_ANNOT}/catalog_226.DELs_not_in_repeats.bed")

# Depth
N_STRONG=0; N_MODERATE=0; N_WEAK=0; N_NODEPTH=0; N_NODATA=0
DEPTH_FILE="${DIR_DEPTH}/depth_support_226.tsv"
if [[ -f "${DEPTH_FILE}" ]]; then
  N_STRONG=$(awk -F'\t' 'NR>1 && $5=="strong_DEL_depth_support"' "${DEPTH_FILE}" | wc -l)
  N_MODERATE=$(awk -F'\t' 'NR>1 && $5=="moderate_DEL_depth_support"' "${DEPTH_FILE}" | wc -l)
  N_WEAK=$(awk -F'\t' 'NR>1 && $5=="weak_DEL_depth_support"' "${DEPTH_FILE}" | wc -l)
  N_NODEPTH=$(awk -F'\t' 'NR>1 && $5=="no_depth_support"' "${DEPTH_FILE}" | wc -l)
  N_NODATA=$(awk -F'\t' 'NR>1 && $5=="depth_no_data"' "${DEPTH_FILE}" | wc -l)
fi

# Mate distance
N_NORMAL=0; N_WARN=0; N_SUSP=0; N_EXTREME=0
MATE_FILE="${DIR_MATDIST}/mate_distance_qc_226.tsv"
if [[ -f "${MATE_FILE}" ]]; then
  N_NORMAL=$(awk -F'\t' 'NR>1 && $7=="normal"' "${MATE_FILE}" | wc -l)
  N_WARN=$(awk -F'\t' 'NR>1 && $7=="warning_large"' "${MATE_FILE}" | wc -l)
  N_SUSP=$(awk -F'\t' 'NR>1 && $7=="very_suspicious"' "${MATE_FILE}" | wc -l)
  N_EXTREME=$(awk -F'\t' 'NR>1 && $7=="extreme_artifact_candidate"' "${MATE_FILE}" | wc -l)
fi

# Per-sample + private/shared
MATRIX_226="${DIR_FINAL}/catalog_226.DEL.GT_matrix.tsv"
PER_SAMPLE="${DIR_SUMMARY}/per_sample_DEL_counts.tsv"
PRIV_SHARED="${DIR_SUMMARY}/private_vs_shared_DEL.tsv"

if [[ -f "${MATRIX_226}" ]]; then
  python3 << PYEOF
matrix = "${MATRIX_226}"
per_out = "${PER_SAMPLE}"
ps_out  = "${PRIV_SHARED}"

with open(matrix) as f:
    header = f.readline().strip().split('\t')
    samples = header[5:]
    ns = len(samples)
    sc = [0]*ns
    sharing = []
    for line in f:
        p = line.strip().split('\t')
        gts = p[5:]
        nr = 0
        for i, gt in enumerate(gts):
            if gt not in ('./.','0/0','./0','0/.','.'):
                sc[i] += 1; nr += 1
        sharing.append(nr)

with open(per_out, 'w') as o:
    o.write("sample\tn_DELs\n")
    for n, c in zip(samples, sc):
        o.write(f"{n}\t{c}\n")

priv = sum(1 for s in sharing if s == 1)
s2_5 = sum(1 for s in sharing if 2 <= s <= 5)
s6_20 = sum(1 for s in sharing if 6 <= s <= 20)
s21p = sum(1 for s in sharing if s > 20)
fixed = sum(1 for s in sharing if s == ns)

with open(ps_out, 'w') as o:
    o.write("category\tcount\n")
    o.write(f"private (1 sample)\t{priv}\n")
    o.write(f"shared (2-5)\t{s2_5}\n")
    o.write(f"shared (6-20)\t{s6_20}\n")
    o.write(f"shared (21+)\t{s21p}\n")
    o.write(f"fixed (all {ns})\t{fixed}\n")
PYEOF
fi

# ── Write report ────────────────────────────────────────────────────────────
cat > "${REPORT}" << REOF
================================================================================
DELLY DEL CATALOG — SUMMARY REPORT
Generated: $(date '+%F %T')
DELLY version: $("${DELLY_BIN}" 2>&1 | grep -oP 'Version: \K[0-9.]+' || echo "1.7.3")
Reference: fClaHyb_Gar_LG.fa (C. gariepinus haplotype)
Cohort: 226 catfish (81 approximately unrelated, theta < 0.177)
Median depth: ~9X (lcWGS)
================================================================================

1. CATALOG SIZES
   226-sample merged DEL catalog:       ${N_226} sites
   81-sample catalog:                    ${N_81} sites
   81-sample PASS DELs:                  ${N_81_PASS} sites

2. FUNCTIONAL OVERLAP (226 catalog)
   CDS-overlapping:                      ${N_CDS}
   Exon-overlapping (non-CDS):           ${N_EXON}
   Intronic:                             ${N_INTRONIC}
   Intergenic:                           ${N_INTERGENIC}
   Unique genes overlapped:              ${N_GENES}

   NOTE: Gene overlap ≠ functional impact. Intronic DELs = uncertain.
   CDS overlap = prioritize cautiously. No phenotype inference.

3. REPEAT OVERLAP (226 catalog)
   ≥50% repeat overlap:                  ${N_REP}
   No major repeat overlap:              ${N_NOREP}

4. DEPTH SUPPORT (226 catalog)
   Strong (ratio < 0.3):                 ${N_STRONG}
   Moderate (0.3–0.6):                   ${N_MODERATE}
   Weak (0.6–0.85):                      ${N_WEAK}
   None (≥ 0.85):                        ${N_NODEPTH}
   No data:                              ${N_NODATA}

5. SVLEN / MATE DISTANCE QC (226 catalog)
   Normal (<20 kb):                      ${N_NORMAL}
   Warning (20–50 kb):                   ${N_WARN}
   Very suspicious (50–100 kb):          ${N_SUSP}
   Extreme artifact (>100 kb):           ${N_EXTREME}

6. SHARING → See: ${PRIV_SHARED}
   Per-sample → See: ${PER_SAMPLE}

================================================================================
CAVEATS: Catalog of putative deletions. ~9X limits het sensitivity.
Gene names ≠ phenotype. CDS overlap ≠ pathogenicity.
Repeat-overlapping ≠ false. No breeding decisions from this alone.
================================================================================
REOF

# Counts TSV
cat > "${COUNTS}" << CEOF
metric	value
total_DEL_226	${N_226}
total_DEL_81	${N_81}
PASS_DEL_81	${N_81_PASS}
CDS_overlap	${N_CDS}
exon_overlap	${N_EXON}
intronic	${N_INTRONIC}
intergenic	${N_INTERGENIC}
unique_genes	${N_GENES}
repeat_50pct	${N_REP}
no_repeat	${N_NOREP}
depth_strong	${N_STRONG}
depth_moderate	${N_MODERATE}
depth_weak	${N_WEAK}
depth_none	${N_NODEPTH}
depth_nodata	${N_NODATA}
mate_normal	${N_NORMAL}
mate_warning	${N_WARN}
mate_suspicious	${N_SUSP}
mate_extreme	${N_EXTREME}
CEOF

dv_log "=== SUMMARY COMPLETE ==="
cat "${REPORT}"
