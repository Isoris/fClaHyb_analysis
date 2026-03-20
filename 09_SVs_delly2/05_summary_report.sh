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

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/00_delly_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

dv_init_dirs
dv_log "=== FINAL SUMMARY ==="

REPORT="${DIR_SUMMARY}/delly_DEL_summary_report.txt"
COUNTS="${DIR_SUMMARY}/delly_DEL_counts.tsv"

# Counts
FINAL_226_VCF="${DIR_FINAL}/catalog_226.DEL.vcf.gz"
N_226=0; [[ -f "${FINAL_226_VCF}" ]] && N_226=$(bcftools view -H "${FINAL_226_VCF}" | wc -l)
FINAL_81_PASS=$(ls "${DIR_FINAL}"/catalog_81.DEL.*.PASS.vcf.gz 2>/dev/null | head -1)
N_81_PASS=0; [[ -n "${FINAL_81_PASS}" && -f "${FINAL_81_PASS}" ]] && N_81_PASS=$(bcftools view -H "${FINAL_81_PASS}" | wc -l)

FUNC="${DIR_ANNOT}/catalog_226.functional_class.tsv"
N_CDS=0;N_EXON=0;N_INTRONIC=0;N_INTERGENIC=0;N_GENES=0
if [[ -f "${FUNC}" ]]; then
  N_CDS=$(awk -F'\t' 'NR>1&&$10=="CDS_overlap"' "${FUNC}"|wc -l)
  N_EXON=$(awk -F'\t' 'NR>1&&$10=="exon_overlap"' "${FUNC}"|wc -l)
  N_INTRONIC=$(awk -F'\t' 'NR>1&&$10=="intronic"' "${FUNC}"|wc -l)
  N_INTERGENIC=$(awk -F'\t' 'NR>1&&$10=="intergenic"' "${FUNC}"|wc -l)
  N_GENES=$(awk -F'\t' 'NR>1&&$10!="intergenic"{print $11}' "${FUNC}"|tr ',' '\n'|grep -v '^\.$'|sort -u|wc -l)
fi

N_REP=0;N_NOREP=0
[[ -f "${DIR_ANNOT}/catalog_226.DELs_in_repeats.bed" ]] && N_REP=$(wc -l < "${DIR_ANNOT}/catalog_226.DELs_in_repeats.bed")
[[ -f "${DIR_ANNOT}/catalog_226.DELs_not_in_repeats.bed" ]] && N_NOREP=$(wc -l < "${DIR_ANNOT}/catalog_226.DELs_not_in_repeats.bed")

# Per-sample + sharing (also generates plotting intermediate tables)
MATRIX_226="${DIR_FINAL}/catalog_226.DEL.GT_matrix.tsv"
PER_SAMPLE="${DIR_SUMMARY}/per_sample_DEL_counts.tsv"
PRIV_SHARED="${DIR_SUMMARY}/private_vs_shared_DEL.tsv"
PER_CHR="${DIR_SUMMARY}/per_chromosome_DEL_counts.tsv"
SVLEN_DIST="${DIR_SUMMARY}/DEL_svlen_distribution.tsv"
WINDOW_COUNTS="${DIR_SUMMARY}/DEL_window_counts_1Mb.tsv"
PAIRWISE_SHARED="${DIR_SUMMARY}/pairwise_shared_DEL.tsv"

if [[ -f "${MATRIX_226}" ]]; then
  python3 << PYEOF
import sys
from collections import defaultdict

matrix = "${MATRIX_226}"
per_out = "${PER_SAMPLE}"
ps_out  = "${PRIV_SHARED}"
chr_out = "${PER_CHR}"
svlen_out = "${SVLEN_DIST}"
window_out = "${WINDOW_COUNTS}"
pairwise_out = "${PAIRWISE_SHARED}"

with open(matrix) as f:
    header = f.readline().strip().split('\t')
    samples = header[5:]
    ns = len(samples)
    sc = [0]*ns        # per-sample DEL count
    sharing = []       # n_samples carrying each DEL
    chr_counts = defaultdict(int)
    svlens = []
    window_data = defaultdict(int)  # (chrom, window_start) -> count
    # For pairwise sharing: build binary genotype matrix
    gt_matrix = []     # list of binary vectors (1 = has DEL)

    for line in f:
        p = line.strip().split('\t')
        chrom, pos, end, vid, svlen = p[0], int(p[1]), int(p[2]), p[3], p[4]
        gts = p[5:]
        binary = []
        nr = 0
        for i, gt in enumerate(gts):
            has_alt = gt not in ('./.','0/0','./0','0/.','.','0|0')
            if has_alt:
                sc[i] += 1; nr += 1
            binary.append(1 if has_alt else 0)
        sharing.append(nr)
        chr_counts[chrom] += 1
        gt_matrix.append(binary)

        # SVLEN
         svlens.append(abs(end - pos) + 1)
#        try:
#            svlens.append(abs(int(svlen)))
#        except:
#            pass

        # 1-Mb window
        w_start = (pos // 1000000) * 1000000
        window_data[(chrom, w_start)] += 1

# Per-sample
with open(per_out, 'w') as o:
    o.write("sample\tn_DELs\n")
    for n,c in zip(samples, sc): o.write(f"{n}\t{c}\n")

# Private vs shared
priv = sum(1 for s in sharing if s==1)
s2_5 = sum(1 for s in sharing if 2<=s<=5)
s6_20 = sum(1 for s in sharing if 6<=s<=20)
s21p = sum(1 for s in sharing if s>20)
fixed = sum(1 for s in sharing if s==ns)
with open(ps_out, 'w') as o:
    o.write("category\tcount\n")
    o.write(f"private\t{priv}\nshared_2_5\t{s2_5}\nshared_6_20\t{s6_20}\nshared_21plus\t{s21p}\nfixed\t{fixed}\n")

# Per chromosome
with open(chr_out, 'w') as o:
    o.write("chrom\tn_DELs\n")
    for c in sorted(chr_counts, key=lambda x: (x.replace('C_gar_LG',''), x)):
        o.write(f"{c}\t{chr_counts[c]}\n")

# SVLEN distribution
with open(svlen_out, 'w') as o:
    o.write("svlen_bp\n")
    for s in svlens: o.write(f"{s}\n")

# 1-Mb window counts (for heatmap panel A)
with open(window_out, 'w') as o:
    o.write("chrom\twindow_start\twindow_end\tn_DELs\n")
    for (c,ws), cnt in sorted(window_data.items()):
        o.write(f"{c}\t{ws}\t{ws+1000000}\t{cnt}\n")

# Pairwise shared DEL counts (for heatmap panel D)
# This can be large for 226 samples, so compute efficiently
print("Computing pairwise shared DELs...")
import numpy as np
gt_arr = np.array(gt_matrix, dtype=np.uint8)  # (n_dels, n_samples)
# Shared = dot product of transposed binary matrix
shared = gt_arr.T @ gt_arr  # (n_samples, n_samples)

with open(pairwise_out, 'w') as o:
    o.write("sample_i\tsample_j\tn_shared\n")
    for i in range(ns):
        for j in range(i, ns):
            o.write(f"{samples[i]}\t{samples[j]}\t{shared[i,j]}\n")

print(f"Per-sample: {per_out}")
print(f"Sharing: {ps_out}")
print(f"Per-chrom: {chr_out}")
print(f"SVLEN: {svlen_out}")
print(f"Windows: {window_out}")
print(f"Pairwise: {pairwise_out}")
PYEOF
fi

# Per-sample genotype matrix for PCA (binary 0/1 presence)
PCA_INPUT="${DIR_SUMMARY}/DEL_binary_genotype_matrix.tsv"
if [[ -f "${MATRIX_226}" ]]; then
  python3 << PYEOF
matrix = "${MATRIX_226}"
out = "${PCA_INPUT}"
with open(matrix) as f, open(out, 'w') as o:
    header = f.readline().strip().split('\t')
    samples = header[5:]
    o.write("DEL_ID\t" + "\t".join(samples) + "\n")
    for line in f:
        p = line.strip().split('\t')
        did = f"{p[0]}:{p[1]}-{p[2]}"
        gts = p[5:]
        binary = []
        for gt in gts:
            binary.append("1" if gt not in ('./.','0/0','./0','0/.','.','0|0') else "0")
        o.write(did + "\t" + "\t".join(binary) + "\n")
PYEOF
fi

# Write report
cat > "${REPORT}" << REOF
================================================================================
DELLY DEL CATALOG — SUMMARY REPORT
Generated: $(date '+%F %T')
Reference: fClaHyb_Gar_LG.fa | Cohort: 226 (81 unrelated) | ~9X lcWGS
================================================================================
226 catalog:     ${N_226} DELs
81 PASS catalog: ${N_81_PASS} DELs

Functional: CDS=${N_CDS} exon=${N_EXON} intronic=${N_INTRONIC} intergenic=${N_INTERGENIC}
Unique genes: ${N_GENES}
Repeat >=50%: ${N_REP} | No repeat: ${N_NOREP}

Caveats: Catalog only. ~9X limits het sensitivity. Gene overlap ≠ phenotype.
================================================================================
REOF

cat > "${COUNTS}" << CEOF
metric	value
total_DEL_226	${N_226}
PASS_DEL_81	${N_81_PASS}
CDS_overlap	${N_CDS}
exon_overlap	${N_EXON}
intronic	${N_INTRONIC}
intergenic	${N_INTERGENIC}
unique_genes	${N_GENES}
repeat_50pct	${N_REP}
no_repeat	${N_NOREP}
CEOF

dv_log "=== SUMMARY COMPLETE ==="
cat "${REPORT}"
