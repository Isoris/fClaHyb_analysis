#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --mem=237GB
#SBATCH -t 2-00:00:00
#SBATCH -J qc_covmap_ctx
#SBATCH -o /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/02-merged_per_sample/results/popgen_qc/logs/qc_covmap_ctx.%j.out
#SBATCH -e /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/02-merged_per_sample/results/popgen_qc/logs/qc_covmap_ctx.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
PA_PROJECT="${BASE}"

MANIFEST="${BASE}/pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv"
REF="${BASE}/00-samples/fClaHyb_Gar_LG.fa"
FAI="${REF}.fai"

PA_REPEAT_BED="${PA_PROJECT}/fClaHyb_Gar_LG.mask_regions.softacgt.bed"
PA_CALLABLE_BED="${PA_PROJECT}/fClaHyb_Gar_LG.mask_regions.normalACGT.renamed.0based.bed"

OUTDIR="${BASE}/02-merged_per_sample/results/popgen_qc"
LOGDIR="${OUTDIR}/logs"
TMPROOT="${OUTDIR}/tmp_covmap_ctx"
WORKDIR="${TMPROOT}/work"
PER_SAMPLE_DIR="${TMPROOT}/per_sample"
PER_CHR_DIR="${TMPROOT}/per_chr"

mkdir -p "$OUTDIR" "$LOGDIR" "$TMPROOT" "$WORKDIR" "$PER_SAMPLE_DIR" "$PER_CHR_DIR"

TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-64}"
JOBS="${JOBS:-8}"
THREADS_PER_JOB="${THREADS_PER_JOB:-4}"

if (( JOBS * THREADS_PER_JOB > TOTAL_CPUS )); then
  echo "[WARN] JOBS*THREADS_PER_JOB exceeds allocated CPUs; resetting conservatively."
  THREADS_PER_JOB=4
  JOBS=$(( TOTAL_CPUS / THREADS_PER_JOB ))
  (( JOBS < 1 )) && JOBS=1
fi

echo "[INFO] TOTAL_CPUS=${TOTAL_CPUS}"
echo "[INFO] JOBS=${JOBS}"
echo "[INFO] THREADS_PER_JOB=${THREADS_PER_JOB}"

command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools not found"; exit 1; }
command -v mosdepth >/dev/null 2>&1 || { echo "[ERROR] mosdepth not found"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "[ERROR] bedtools not found"; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo "[ERROR] bgzip not found"; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "[ERROR] tabix not found"; exit 1; }
python3 - <<'PY'
import importlib.util, sys
mods = ["gzip", "csv", "statistics", "collections"]
for m in mods:
    if importlib.util.find_spec(m) is None:
        print(f"[ERROR] Missing python module: {m}")
        sys.exit(1)
print("[INFO] Python stdlib OK")
PY

if [[ ! -s "$FAI" ]]; then
  samtools faidx "$REF"
fi

# -----------------------------
# 1. Build reference-wide chromosome BED
# -----------------------------
CHR_BED="${WORKDIR}/ref.chromosomes.bed"
awk 'BEGIN{OFS="\t"} {print $1,0,$2,$1}' "$FAI" > "$CHR_BED"

# -----------------------------
# 2. Sort + complement region beds
# -----------------------------
REPEAT_SORTED="${WORKDIR}/repeat.sorted.bed"
CALLABLE_SORTED="${WORKDIR}/callable.sorted.bed"
NONREPEAT_BED="${WORKDIR}/nonrepeat.bed"
NONCALLABLE_BED="${WORKDIR}/noncallable.bed"

sort -k1,1 -k2,2n "$PA_REPEAT_BED" > "$REPEAT_SORTED"
sort -k1,1 -k2,2n "$PA_CALLABLE_BED" > "$CALLABLE_SORTED"

bedtools complement -i "$REPEAT_SORTED" -g "$FAI" > "$NONREPEAT_BED"
bedtools complement -i "$CALLABLE_SORTED" -g "$FAI" > "$NONCALLABLE_BED"

# -----------------------------
# 3. Combined mosdepth regions BED
#    One run per BAM, all region classes together
# -----------------------------
COMBINED_BED="${WORKDIR}/combined_regions.bed"

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"chr_all|" $1 "|" NR}' "$CHR_BED" > "$COMBINED_BED"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"repeat|" $1 "|" NR}' "$REPEAT_SORTED" >> "$COMBINED_BED"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"nonrepeat|" $1 "|" NR}' "$NONREPEAT_BED" >> "$COMBINED_BED"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"callable|" $1 "|" NR}' "$CALLABLE_SORTED" >> "$COMBINED_BED"
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"noncallable|" $1 "|" NR}' "$NONCALLABLE_BED" >> "$COMBINED_BED"

sort -k1,1 -k2,2n "$COMBINED_BED" -o "$COMBINED_BED"

# -----------------------------
# 4. Sample list from manifest
# -----------------------------
SAMPLE_BAMS="${WORKDIR}/sample_bams.tsv"
awk -F'\t' 'NR>1 && $1!="" && $3!="" {print $1 "\t" $3}' "$MANIFEST" > "$SAMPLE_BAMS"

NSAMPLES=$(wc -l < "$SAMPLE_BAMS")
echo "[INFO] Samples in manifest: ${NSAMPLES}"
if [[ "$NSAMPLES" -eq 0 ]]; then
  echo "[ERROR] No sample/BAM pairs found in manifest"
  exit 1
fi

# -----------------------------
# 5. Worker script
# -----------------------------
WORKER="${WORKDIR}/worker_qc_covmap.sh"
cat > "$WORKER" <<'BASH'
#!/usr/bin/env bash
set -euo pipefail

sample="$1"
bam="$2"
base="$3"
combined_bed="$4"
per_sample_dir="$5"
per_chr_dir="$6"
threads="$7"
workdir="$8"

prefix="${workdir}/${sample}"

if [[ ! -s "$bam" ]]; then
  echo "[WARN] Missing BAM for ${sample}: ${bam}" >&2
  exit 0
fi

# A. mosdepth over all region classes in one run
mosdepth -n -t "$threads" --by "$combined_bed" --thresholds 1,5,10 "$prefix" "$bam" >/dev/null 2>&1

# B. MAPQ distribution per chromosome from mapped primary alignments
#    Final filtered BAMs should already exclude much junk, but we still keep the filters explicit.
samtools view -@ "$threads" -F 0x904 "$bam" \
| awk 'BEGIN{OFS="\t"}
{
  chr=$3
  total++
  tot[chr]++
  if ($5==0) {mq0++; mq0c[chr]++}
  if ($5>=30) {mq30++; mq30c[chr]++}
}
END{
  print "TOTAL", total+0, mq0+0, mq30+0
  for (c in tot) print c, tot[c]+0, mq0c[c]+0, mq30c[c]+0
}' > "${prefix}.mapq_counts.tsv"

# C. Aggregate mosdepth outputs + mapq counts into one per-sample overall row + one per-chr table
python3 - <<'PY' "$sample" "$bam" "${prefix}.regions.bed.gz" "${prefix}.thresholds.bed.gz" "${prefix}.mapq_counts.tsv" "${per_sample_dir}/${sample}.overall.tsv" "${per_chr_dir}/${sample}.per_chrom.tsv"
import sys, gzip, csv
from collections import defaultdict

sample, bam, regions_gz, thresholds_gz, mapq_tsv, overall_out, perchr_out = sys.argv[1:]

# region stats accumulator
# key = (class_name, chrom)
acc = defaultdict(lambda: {"length":0, "depth_sum":0.0, "cov1":0, "cov5":0, "cov10":0})

with gzip.open(regions_gz, "rt") as f:
    for line in f:
        if not line.strip():
            continue
        chrom, start, end, name, mean_cov = line.rstrip("\n").split("\t")[:5]
        class_name, chrom_from_name, _ = name.split("|", 2)
        start, end = int(start), int(end)
        length = end - start
        acc[(class_name, chrom_from_name)]["length"] += length
        acc[(class_name, chrom_from_name)]["depth_sum"] += float(mean_cov) * length

with gzip.open(thresholds_gz, "rt") as f:
    for line in f:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        chrom, start, end, name = parts[:4]
        vals = parts[4:]
        class_name, chrom_from_name, _ = name.split("|", 2)
        # mosdepth thresholds are counts of bases >= threshold
        c1 = int(vals[0]) if len(vals) > 0 else 0
        c5 = int(vals[1]) if len(vals) > 1 else 0
        c10 = int(vals[2]) if len(vals) > 2 else 0
        acc[(class_name, chrom_from_name)]["cov1"] += c1
        acc[(class_name, chrom_from_name)]["cov5"] += c5
        acc[(class_name, chrom_from_name)]["cov10"] += c10

mapq = {}
with open(mapq_tsv) as f:
    for line in f:
        chrom, total, mq0, mq30 = line.rstrip("\n").split("\t")
        total = int(total)
        mq0 = int(mq0)
        mq30 = int(mq30)
        mapq[chrom] = {
            "mapped_primary_reads": total,
            "mapq0_reads": mq0,
            "mapq30plus_reads": mq30,
            "mapq0_frac": (mq0 / total) if total > 0 else 0.0,
            "mapq30plus_frac": (mq30 / total) if total > 0 else 0.0,
        }

classes = ["chr_all", "repeat", "nonrepeat", "callable", "noncallable"]

def metrics_for(class_name, chrom):
    d = acc.get((class_name, chrom), None)
    if not d or d["length"] == 0:
        return {
            f"{class_name}_region_bp": 0,
            f"{class_name}_mean_depth_x": "",
            f"{class_name}_breadth_1x": "",
            f"{class_name}_breadth_5x": "",
            f"{class_name}_breadth_10x": "",
        }
    L = d["length"]
    return {
        f"{class_name}_region_bp": L,
        f"{class_name}_mean_depth_x": d["depth_sum"] / L,
        f"{class_name}_breadth_1x": d["cov1"] / L,
        f"{class_name}_breadth_5x": d["cov5"] / L,
        f"{class_name}_breadth_10x": d["cov10"] / L,
    }

# per chromosome rows
chroms = sorted({chrom for (_, chrom) in acc.keys() if chrom != "TOTAL"})
with open(perchr_out, "w", newline="") as out:
    fieldnames = [
        "sample", "bam", "chrom",
        "mapped_primary_reads", "mapq0_reads", "mapq30plus_reads", "mapq0_frac", "mapq30plus_frac"
    ]
    for c in classes:
        fieldnames.extend([
            f"{c}_region_bp",
            f"{c}_mean_depth_x",
            f"{c}_breadth_1x",
            f"{c}_breadth_5x",
            f"{c}_breadth_10x",
        ])
    w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()

    for chrom in chroms:
        row = {
            "sample": sample,
            "bam": bam,
            "chrom": chrom,
            "mapped_primary_reads": mapq.get(chrom, {}).get("mapped_primary_reads", 0),
            "mapq0_reads": mapq.get(chrom, {}).get("mapq0_reads", 0),
            "mapq30plus_reads": mapq.get(chrom, {}).get("mapq30plus_reads", 0),
            "mapq0_frac": mapq.get(chrom, {}).get("mapq0_frac", 0.0),
            "mapq30plus_frac": mapq.get(chrom, {}).get("mapq30plus_frac", 0.0),
        }
        for c in classes:
            row.update(metrics_for(c, chrom))
        w.writerow(row)

# overall row
# aggregate across all chromosomes by class
overall_acc = defaultdict(lambda: {"length":0, "depth_sum":0.0, "cov1":0, "cov5":0, "cov10":0})
for (class_name, chrom), d in acc.items():
    overall_acc[class_name]["length"] += d["length"]
    overall_acc[class_name]["depth_sum"] += d["depth_sum"]
    overall_acc[class_name]["cov1"] += d["cov1"]
    overall_acc[class_name]["cov5"] += d["cov5"]
    overall_acc[class_name]["cov10"] += d["cov10"]

overall_row = {
    "sample": sample,
    "bam": bam,
    "mapped_primary_reads": mapq.get("TOTAL", {}).get("mapped_primary_reads", 0),
    "mapq0_reads": mapq.get("TOTAL", {}).get("mapq0_reads", 0),
    "mapq30plus_reads": mapq.get("TOTAL", {}).get("mapq30plus_reads", 0),
    "mapq0_frac": mapq.get("TOTAL", {}).get("mapq0_frac", 0.0),
    "mapq30plus_frac": mapq.get("TOTAL", {}).get("mapq30plus_frac", 0.0),
}

for c in classes:
    d = overall_acc[c]
    if d["length"] == 0:
        overall_row[f"{c}_region_bp"] = 0
        overall_row[f"{c}_mean_depth_x"] = ""
        overall_row[f"{c}_breadth_1x"] = ""
        overall_row[f"{c}_breadth_5x"] = ""
        overall_row[f"{c}_breadth_10x"] = ""
    else:
        L = d["length"]
        overall_row[f"{c}_region_bp"] = L
        overall_row[f"{c}_mean_depth_x"] = d["depth_sum"] / L
        overall_row[f"{c}_breadth_1x"] = d["cov1"] / L
        overall_row[f"{c}_breadth_5x"] = d["cov5"] / L
        overall_row[f"{c}_breadth_10x"] = d["cov10"] / L

fieldnames = list(overall_row.keys())
with open(overall_out, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerow(overall_row)
PY
BASH
chmod +x "$WORKER"

# -----------------------------
# 6. Run workers in parallel
# -----------------------------
running=0
while IFS=$'\t' read -r sample bam; do
  echo "[INFO] Launching ${sample}"
  bash "$WORKER" "$sample" "$bam" "$BASE" "$COMBINED_BED" "$PER_SAMPLE_DIR" "$PER_CHR_DIR" "$THREADS_PER_JOB" "$WORKDIR" &
  running=$((running+1))
  if (( running >= JOBS )); then
    wait -n
    running=$((running-1))
  fi
done < "$SAMPLE_BAMS"
wait

# -----------------------------
# 7. Merge per-sample files
# -----------------------------
ALL_OVERALL="${OUTDIR}/all_samples.coverage_context_qc.tsv"
ALL_PERCHR="${OUTDIR}/all_samples.per_chrom.coverage_context_qc.tsv"

first=1
for f in "$PER_SAMPLE_DIR"/*.overall.tsv; do
  [[ -e "$f" ]] || continue
  if (( first == 1 )); then
    cat "$f" > "$ALL_OVERALL"
    first=0
  else
    tail -n +2 "$f" >> "$ALL_OVERALL"
  fi
done

first=1
for f in "$PER_CHR_DIR"/*.per_chrom.tsv; do
  [[ -e "$f" ]] || continue
  if (( first == 1 )); then
    cat "$f" > "$ALL_PERCHR"
    first=0
  else
    tail -n +2 "$f" >> "$ALL_PERCHR"
  fi
done

# stable sorting
python3 - <<'PY' "$ALL_OVERALL" "$ALL_PERCHR"
import sys, csv

def sort_tsv(path, keys):
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
        fieldnames = r.fieldnames
    rows.sort(key=lambda x: tuple(x[k] for k in keys))
    with open(path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

sort_tsv(sys.argv[1], ["sample"])
sort_tsv(sys.argv[2], ["chrom", "sample"])
PY

# -----------------------------
# 8. Summary tables
# -----------------------------
SUMMARY_OVERALL="${OUTDIR}/summary.coverage_context_qc.tsv"
SUMMARY_PERCHR="${OUTDIR}/summary.per_chrom.coverage_context_qc.tsv"

python3 - <<'PY' "$ALL_OVERALL" "$ALL_PERCHR" "$SUMMARY_OVERALL" "$SUMMARY_PERCHR"
import sys, csv, statistics
from collections import defaultdict

overall_tsv, perchr_tsv, summary_overall, summary_perchr = sys.argv[1:]

def is_number(x):
    if x is None or x == "":
        return False
    try:
        float(x)
        return True
    except:
        return False

def qstats(vals):
    vals = sorted(vals)
    if not vals:
        return None
    def q(p):
        k = max(1, int(round(p * len(vals))))
        k = min(k, len(vals))
        return vals[k-1]
    return {
        "n": len(vals),
        "min": vals[0],
        "p10": q(0.10),
        "p25": q(0.25),
        "p50": q(0.50),
        "mean": sum(vals)/len(vals),
        "p75": q(0.75),
        "p90": q(0.90),
        "max": vals[-1],
    }

# overall summary
with open(overall_tsv) as f:
    r = csv.DictReader(f, delimiter="\t")
    rows = list(r)
    fields = r.fieldnames

numeric_fields = [x for x in fields if x not in ("sample", "bam")]
with open(summary_overall, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["metric","n","min","p10","p25","p50","mean","p75","p90","max"])
    for metric in numeric_fields:
        vals = [float(row[metric]) for row in rows if is_number(row[metric])]
        qs = qstats(vals)
        if qs is None:
            continue
        w.writerow([metric, qs["n"], qs["min"], qs["p10"], qs["p25"], qs["p50"], qs["mean"], qs["p75"], qs["p90"], qs["max"]])

# per chromosome summary
by_chr_metric = defaultdict(list)
with open(perchr_tsv) as f:
    r = csv.DictReader(f, delimiter="\t")
    fields = r.fieldnames
    for row in r:
        chrom = row["chrom"]
        for metric in fields:
            if metric in ("sample", "bam", "chrom"):
                continue
            if is_number(row[metric]):
                by_chr_metric[(chrom, metric)].append(float(row[metric]))

with open(summary_perchr, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["chrom","metric","n","min","p10","p25","p50","mean","p75","p90","max"])
    for (chrom, metric), vals in sorted(by_chr_metric.items()):
        qs = qstats(vals)
        if qs is None:
            continue
        w.writerow([chrom, metric, qs["n"], qs["min"], qs["p10"], qs["p25"], qs["p50"], qs["mean"], qs["p75"], qs["p90"], qs["max"]])
PY

echo "[DONE] Per-sample overall table:   $ALL_OVERALL"
echo "[DONE] Overall summary table:     $SUMMARY_OVERALL"
echo "[DONE] Per-sample per-chr table:  $ALL_PERCHR"
echo "[DONE] Per-chr summary table:     $SUMMARY_PERCHR"
echo "[DONE] Logs:                      $LOGDIR"
