VCF="07_final_catalogs/catalog_226.DEL.vcf.gz"
BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
IGV="${BASE}/pa_roary/IGV_2.19.7/igv.sh"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
GFF="${BASE}/00-samples/fClaHyb_Gar_LG.from_CGAR.gff3_polished.sorted.gff3.gz"
BAMBASE="${BASE}/02-merged_per_sample"

SAMPLE_HOM=$(bcftools query -i 'ID="DEL00002309"' -f '[%SAMPLE\t%GT\n]' "$VCF" | awk '$2=="1/1"{print $1; exit}')
SAMPLE_HET=$(bcftools query -i 'ID="DEL00002309"' -f '[%SAMPLE\t%GT\n]' "$VCF" | awk '$2=="0/1"{print $1; exit}')
SAMPLE_REF=$(bcftools query -i 'ID="DEL00002309"' -f '[%SAMPLE\t%GT\n]' "$VCF" | awk '$2=="0/0"{print $1; exit}')

echo "HOM=$SAMPLE_HOM"
echo "HET=$SAMPLE_HET"
echo "REF=$SAMPLE_REF"

IGV_JAVA_OPTS="-Dsun.java2d.uiScale=1.5" "$IGV" \
  -g "$REF" \
  "${BAMBASE}/${SAMPLE_HOM}/${SAMPLE_HOM}.merged.markdup.clip.pp.samechr.tlenP99.filtered.bam" \
  "${BAMBASE}/${SAMPLE_HET}/${SAMPLE_HET}.merged.markdup.clip.pp.samechr.tlenP99.filtered.bam" \
  "${BAMBASE}/${SAMPLE_REF}/${SAMPLE_REF}.merged.markdup.clip.pp.samechr.tlenP99.filtered.bam" \
  "${BASE}/delly_sv/07_final_catalogs/catalog_226.DEL.vcf.gz" \
  "${BASE}/delly_sv/07_final_catalogs/catalog_226.DEL.bed" \
  "$GFF"
