#!/usr/bin/env bash
set -euo pipefail

ROOT="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
LIST="${ROOT}/popstruct_thin/list_of_samples_one_per_line_same_bamfile_list.tsv"
OUT="${ROOT}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"

mkdir -p "$(dirname "$OUT")"

{
  echo -e "Sample\tBAMminimap2\tBAMP99TLENMAPQ30"

  while IFS= read -r sample; do
    [[ -z "$sample" ]] && continue

    raw_bam="$(find "${ROOT}/01-minimap2-bams" \
      -type f -name "${sample}.*.sr.bam" \
      ! -name "*.tmp.*.bam" \
      ! -path "*/qc_dedup/*" \
      | sort | head -n 1)"

    filt_bam="$(find "${ROOT}/02-merged_per_sample/${sample}" \
      -type f -name "*.filtered.bam" \
      | sort | head -n 1)"

    echo -e "${sample}\t${raw_bam:-NA}\t${filt_bam:-NA}"
  done < "$LIST"
} > "$OUT"

echo "Wrote: $OUT"
