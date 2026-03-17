# ============================================================
# ANGSD per-sample heterozygosity + local theta tracks
# Produces:
#   1) boxplot-ready per-sample heterozygosity table
#   2) per-sample theta windows (for ideograms)
#
# Assumes:
#   - You have one BAM per sample
#   - You have ref.fa (no true ancestral) -> use folded mode
#   - You have a callable/nonrepeat BED mask (recommended)
#
# Requirements:
#   angsd, realSFS, thetaStat (from ANGSD), bedtools, bgzip+tabix
# ============================================================

# ---------------------------
# USER EDIT
# ---------------------------
REF="ref.fa"
BAMLIST="bam_list.txt"                 # one BAM path per line
MASK_BED="callable_nonrepeat.bed"      # recommended: callable minus repeats
OUTDIR="ANGSD_HET_THETA"
THREADS=8

# theta windows for ideogram
WIN=500000          # 500 kb ideogram bins (change to 1000000 for 1 Mb)
STEP=500000

# QC filters (reasonable defaults)
MINQ=20
MINMAPQ=30
CLIP=50

mkdir -p "${OUTDIR}"/{01_saf,02_sfs,03_theta,04_summary,logs}

echo -e "sample\thet_saf1" > "${OUTDIR}/04_summary/heterozygosity_per_sample.tsv"

# ---------------------------
# LOOP: per-sample SAF -> SFS -> HET -> local theta
# ---------------------------
while read -r BAM; do
  [[ -z "${BAM}" ]] && continue
  [[ "${BAM}" =~ ^# ]] && continue

  SAMPLE=$(basename "${BAM}")
  SAMPLE=${SAMPLE%.bam}
  echo "==> ${SAMPLE}"

  # 1) SAF (folded, because no ancestral)
  # Note: -fold 1 requires -anc; using ref as proxy is standard.
  # Masking: using -sites to restrict to callable/nonrepeat sites is strongly recommended.
  angsd \
    -i "${BAM}" \
    -ref "${REF}" \
    -anc "${REF}" \
    -GL 1 \
    -doSaf 1 \
    -fold 1 \
    -minQ "${MINQ}" \
    -minMapQ "${MINMAPQ}" \
    -C "${CLIP}" \
    -sites "${MASK_BED}" \
    -P "${THREADS}" \
    -out "${OUTDIR}/01_saf/${SAMPLE}" \
    2> "${OUTDIR}/logs/${SAMPLE}.saf.log"

  # 2) Single-sample SFS
  realSFS "${OUTDIR}/01_saf/${SAMPLE}.saf.idx" > "${OUTDIR}/02_sfs/${SAMPLE}.est.ml"

  # 3) Heterozygosity = second bin / sum bins (diploid single sample)
  HET=$(awk '{
      sum=0; for(i=1;i<=NF;i++) sum+=$i;
      if(sum>0) printf "%.12g", $2/sum; else printf "NA";
    }' "${OUTDIR}/02_sfs/${SAMPLE}.est.ml")

  echo -e "${SAMPLE}\t${HET}" >> "${OUTDIR}/04_summary/heterozygosity_per_sample.tsv"

  # 4) Local theta (per-site thetas + window summary)
  #
  # Use -pest sample.est.ml so ANGSD outputs posterior-based thetas for that sample.
  # This is the "more work" step needed for ideograms.
  angsd \
    -i "${BAM}" \
    -ref "${REF}" \
    -anc "${REF}" \
    -GL 1 \
    -doSaf 1 \
    -fold 1 \
    -pest "${OUTDIR}/02_sfs/${SAMPLE}.est.ml" \
    -doThetas 1 \
    -minQ "${MINQ}" \
    -minMapQ "${MINMAPQ}" \
    -C "${CLIP}" \
    -sites "${MASK_BED}" \
    -P "${THREADS}" \
    -out "${OUTDIR}/03_theta/${SAMPLE}" \
    2> "${OUTDIR}/logs/${SAMPLE}.theta.log"

  # 4a) Make a BED-style per-site theta file (thetaStat utility)
  # Produces files like: ${SAMPLE}.thetas.gz + ${SAMPLE}.thetas.idx (ANGSD) and bed outputs (thetaStat)
  thetaStat make_bed \
    -outnames "${OUTDIR}/03_theta/${SAMPLE}" \
    "${OUTDIR}/03_theta/${SAMPLE}.thetas.idx"

  # 4b) Windowed theta for ideograms (TSV)
  # This is the easiest thing to plot as an ideogram track.
  thetaStat do_stat \
    "${OUTDIR}/03_theta/${SAMPLE}.thetas.idx" \
    -win "${WIN}" -step "${STEP}" \
    > "${OUTDIR}/03_theta/${SAMPLE}.theta.win${WIN}.step${STEP}.tsv"

done < "${BAMLIST}"

echo "DONE:"
echo "  Boxplot table: ${OUTDIR}/04_summary/heterozygosity_per_sample.tsv"
echo "  Per-sample theta windows: ${OUTDIR}/03_theta/*.theta.win${WIN}.step${STEP}.tsv"
