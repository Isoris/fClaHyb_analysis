###
# summarize_bam_qc.slurm does this:
# Parses per-sample samtools stats and flagstat outputs from merged BAMs,
# writes a per-sample QC table, and generates a summary report with
# depth distributions and heuristic ANGSD depth cutoff suggestions.
###
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8GB
#SBATCH -t 01:00:00
#SBATCH -J bam_qc_sum
#SBATCH -o bam_qc_sum.%j.out
#SBATCH -e bam_qc_sum.%j.err

set -euo pipefail
source ~/.bashrc
mamba activate assembly >/dev/null 2>&1 || true

# Run this from: /scratch/lt200308-agbsci/Quentin_project
BASE="$(pwd)"
MERGED_DIR="${BASE}/02-merged_per_sample"

# Genome size you gave:
GENOME_SIZE=963905721

OUTDIR="${MERGED_DIR}/results/popgen_qc"
mkdir -p "$OUTDIR"

python3 - <<'PY' "$MERGED_DIR" "$OUTDIR" "$GENOME_SIZE"
import os, re, sys, glob, statistics
from math import floor

merged_dir, outdir, genome_size = sys.argv[1], sys.argv[2], int(sys.argv[3])

def parse_samtools_stats(path):
    """
    Parse SN lines from samtools stats output.
    Returns dict with keys: raw_total_sequences, reads_properly_paired, reads_paired,
    reads_duplicated, reads_MQ0, bases_mapped_cigar, error_rate
    """
    d = {}
    with open(path, "r", errors="replace") as f:
        for line in f:
            if not line.startswith("SN\t"):
                continue
            # SN <tab> key: <tab> value
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            key = parts[1].strip().rstrip(":")
            val = parts[2].strip()
            # map keys we care about
            if key == "raw total sequences":
                d["raw_total_sequences"] = int(val)
            elif key == "reads properly paired":
                d["reads_properly_paired"] = int(val)
            elif key == "reads paired":
                d["reads_paired"] = int(val)
            elif key == "reads duplicated":
                d["reads_duplicated"] = int(val)
            elif key == "reads MQ0":
                d["reads_MQ0"] = int(val)
            elif key == "bases mapped (cigar)":
                d["bases_mapped_cigar"] = int(val)
            elif key == "error rate":
                try:
                    d["error_rate"] = float(val)
                except:
                    pass
    return d

def parse_flagstat(path):
    """
    Parse key numbers + properly paired % from samtools flagstat.
    Returns dict with keys: total_reads, primary, supplementary, duplicates,
    mapped_pct, properly_paired_pct, mate_mapped_diff_chr
    """
    d = {}
    # Example lines:
    # 54073372 + 0 in total (QC-passed reads + QC-failed reads)
    # 53945835 + 0 primary
    # 127537 + 0 supplementary
    # 0 + 0 duplicates
    # 52520718 + 0 properly paired (97.36% : N/A)
    # 1234894 + 0 with mate mapped to a different chr
    rx_total = re.compile(r"^(\d+)\s+\+\s+\d+\s+in total")
    rx_primary = re.compile(r"^(\d+)\s+\+\s+\d+\s+primary$")
    rx_supp = re.compile(r"^(\d+)\s+\+\s+\d+\s+supplementary$")
    rx_dup = re.compile(r"^(\d+)\s+\+\s+\d+\s+duplicates$")
    rx_mapped = re.compile(r"^(\d+)\s+\+\s+\d+\s+mapped\s+\(([\d\.]+)%")
    rx_pp = re.compile(r"^(\d+)\s+\+\s+\d+\s+properly paired\s+\(([\d\.]+)%")
    rx_diffchr = re.compile(r"^(\d+)\s+\+\s+\d+\s+with mate mapped to a different chr$")

    with open(path, "r", errors="replace") as f:
        for line in f:
            line = line.strip()
            m = rx_total.match(line)
            if m: d["total_reads"] = int(m.group(1)); continue
            m = rx_primary.match(line)
            if m: d["primary"] = int(m.group(1)); continue
            m = rx_supp.match(line)
            if m: d["supplementary"] = int(m.group(1)); continue
            m = rx_dup.match(line)
            if m: d["duplicates"] = int(m.group(1)); continue
            m = rx_mapped.match(line)
            if m:
                d["mapped"] = int(m.group(1))
                d["mapped_pct"] = float(m.group(2))
                continue
            m = rx_pp.match(line)
            if m:
                d["properly_paired"] = int(m.group(1))
                d["properly_paired_pct"] = float(m.group(2))
                continue
            m = rx_diffchr.match(line)
            if m: d["mate_mapped_diff_chr"] = int(m.group(1)); continue
    return d

def quantiles(xs):
    xs = sorted(xs)
    if not xs:
        return {}
    def q(p):
        # simple nearest-rank quantile
        k = max(1, int(round(p*len(xs))))
        k = min(k, len(xs))
        return xs[k-1]
    return {
        "min": xs[0],
        "p10": q(0.10),
        "p25": q(0.25),
        "p50": q(0.50),
        "p75": q(0.75),
        "p90": q(0.90),
        "max": xs[-1],
        "mean": sum(xs)/len(xs),
        "median": statistics.median(xs),
        "n": len(xs),
    }

# collect sample dirs CGA*
sample_dirs = sorted([d for d in glob.glob(os.path.join(merged_dir, "CGA*")) if os.path.isdir(d)])
rows = []
missing = []

for sdir in sample_dirs:
    sample = os.path.basename(sdir)

    # prefer filtered.flagstat if present
    flagstat_path = os.path.join(sdir, f"{sample}.filtered.flagstat.txt")
    if not os.path.exists(flagstat_path):
        # fallback to non-filtered
        flagstat_path = os.path.join(sdir, f"{sample}.flagstat.txt")

    stats_path = os.path.join(sdir, f"{sample}.samtools.stats.txt")

    if not os.path.exists(stats_path) and not os.path.exists(flagstat_path):
        missing.append(sample)
        continue

    st = parse_samtools_stats(stats_path) if os.path.exists(stats_path) else {}
    fs = parse_flagstat(flagstat_path) if os.path.exists(flagstat_path) else {}

    # mean depth from bases_mapped_cigar / genome_size
    bases_cigar = st.get("bases_mapped_cigar", None)
    mean_depth = (bases_cigar / genome_size) if bases_cigar else None

    raw_total = st.get("raw_total_sequences", None)
    reads_dup = st.get("reads_duplicated", None)
    reads_mq0 = st.get("reads_MQ0", None)

    dup_frac = (reads_dup / raw_total) if (raw_total and reads_dup is not None) else None
    mq0_frac = (reads_mq0 / raw_total) if (raw_total and reads_mq0 is not None) else None

    rows.append({
        "sample": sample,
        "stats_file": os.path.basename(stats_path) if os.path.exists(stats_path) else "",
        "flagstat_file": os.path.basename(flagstat_path) if os.path.exists(flagstat_path) else "",
        "mean_depth_x": mean_depth,
        "bases_mapped_cigar": bases_cigar,
        "raw_total_sequences": raw_total,
        "reads_duplicated": reads_dup,
        "dup_frac": dup_frac,
        "reads_MQ0": reads_mq0,
        "mq0_frac": mq0_frac,
        "error_rate": st.get("error_rate", None),
        "filtered_total_reads": fs.get("total_reads", None),
        "filtered_primary": fs.get("primary", None),
        "filtered_supplementary": fs.get("supplementary", None),
        "filtered_duplicates": fs.get("duplicates", None),
        "filtered_mapped_pct": fs.get("mapped_pct", None),
        "filtered_properly_paired_pct": fs.get("properly_paired_pct", None),
        "filtered_mate_diff_chr": fs.get("mate_mapped_diff_chr", None),
    })

# write per-sample TSV
tsv_path = os.path.join(outdir, "all_samples.bam_qc.tsv")
with open(tsv_path, "w") as out:
    cols = [
        "sample",
        "mean_depth_x",
        "bases_mapped_cigar",
        "raw_total_sequences",
        "reads_duplicated",
        "dup_frac",
        "reads_MQ0",
        "mq0_frac",
        "error_rate",
        "filtered_total_reads",
        "filtered_primary",
        "filtered_supplementary",
        "filtered_duplicates",
        "filtered_mapped_pct",
        "filtered_properly_paired_pct",
        "filtered_mate_diff_chr",
        "stats_file",
        "flagstat_file",
    ]
    out.write("\t".join(cols) + "\n")
    for r in rows:
        line = []
        for c in cols:
            v = r.get(c, None)
            if v is None:
                line.append("")
            elif isinstance(v, float):
                line.append(f"{v:.6g}")
            else:
                line.append(str(v))
        out.write("\t".join(line) + "\n")

# summary + suggested ANGSD depth cutoffs
depths = [r["mean_depth_x"] for r in rows if r["mean_depth_x"] is not None]
ppcts  = [r["filtered_properly_paired_pct"] for r in rows if r["filtered_properly_paired_pct"] is not None]
dups   = [r["dup_frac"] for r in rows if r["dup_frac"] is not None]
mq0s   = [r["mq0_frac"] for r in rows if r["mq0_frac"] is not None]
errs   = [r["error_rate"] for r in rows if r["error_rate"] is not None]

qd = quantiles(depths)
qpp = quantiles(ppcts)
qdup = quantiles(dups)
qmq0 = quantiles(mq0s)
qerr = quantiles(errs)

# heuristic recommendations (you can override manually)
# - minDepthInd: 3 if p10 depth >= 6 else 2 (very common for ANGSD structure)
# - maxDepthInd: ~5x p90 depth, capped at 60
rec_minDepthInd = 3 if (qd.get("p10") is not None and qd["p10"] >= 6.0) else 2
rec_maxDepthInd = None
if qd.get("p90") is not None:
    rec_maxDepthInd = min(60, int(round(qd["p90"] * 5)))

summary_path = os.path.join(outdir, "summary.txt")
with open(summary_path, "w") as out:
    out.write(f"Genome size (bp): {genome_size}\n")
    out.write(f"Samples parsed   : {len(rows)}\n")
    if missing:
        out.write(f"Samples missing stats+flagstat: {len(missing)}\n")
        out.write("Missing (first 30): " + ", ".join(missing[:30]) + ("\n" if len(missing)<=30 else ", ...\n"))
    out.write("\n=== Mean depth (X) from bases_mapped(cigar)/genome_size ===\n")
    if qd:
        for k in ["n","min","p10","p25","p50","mean","p75","p90","max"]:
            out.write(f"{k}\t{qd[k]:.6g}\n" if isinstance(qd[k], float) else f"{k}\t{qd[k]}\n")
    else:
        out.write("No depth values found.\n")

    out.write("\n=== Properly paired % (from *.filtered.flagstat.txt if present) ===\n")
    if qpp:
        for k in ["n","min","p10","p25","p50","mean","p75","p90","max"]:
            out.write(f"{k}\t{qpp[k]:.6g}\n" if isinstance(qpp[k], float) else f"{k}\t{qpp[k]}\n")
    else:
        out.write("No properly paired % values found.\n")

    out.write("\n=== Duplicate fraction (reads_duplicated/raw_total_sequences; from samtools.stats) ===\n")
    if qdup:
        for k in ["n","min","p10","p25","p50","mean","p75","p90","max"]:
            out.write(f"{k}\t{qdup[k]:.6g}\n" if isinstance(qdup[k], float) else f"{k}\t{qdup[k]}\n")
    else:
        out.write("No duplicate fraction values found.\n")

    out.write("\n=== MQ0 fraction (reads_MQ0/raw_total_sequences; from samtools.stats) ===\n")
    if qmq0:
        for k in ["n","min","p10","p25","p50","mean","p75","p90","max"]:
            out.write(f"{k}\t{qmq0[k]:.6g}\n" if isinstance(qmq0[k], float) else f"{k}\t{qmq0[k]}\n")
    else:
        out.write("No MQ0 fraction values found.\n")

    out.write("\n=== Error rate (from samtools.stats) ===\n")
    if qerr:
        for k in ["n","min","p10","p25","p50","mean","p75","p90","max"]:
            out.write(f"{k}\t{qerr[k]:.6g}\n" if isinstance(qerr[k], float) else f"{k}\t{qerr[k]}\n")
    else:
        out.write("No error rate values found.\n")

    out.write("\n=== Suggested ANGSD per-individual depth cutoffs (heuristic) ===\n")
    out.write(f"setMinDepthInd\t{rec_minDepthInd}\n")
    out.write(f"setMaxDepthInd\t{rec_maxDepthInd if rec_maxDepthInd is not None else ''}\n")
    out.write("\nNotes:\n")
    out.write("- Depth here is computed from samtools.stats 'bases mapped (cigar)'.\n")
    out.write("- If your samtools.stats were generated on pre-filter BAMs, dup/MQ0 fractions reflect that BAM.\n")
    out.write("- Properly paired % comes from filtered.flagstat when available.\n")

print("[DONE] Wrote:")
print("  ", tsv_path)
print("  ", summary_path)
PY

echo "Outputs:"
echo "  $OUTDIR/all_samples.bam_qc.tsv"
echo "  $OUTDIR/summary.txt"
