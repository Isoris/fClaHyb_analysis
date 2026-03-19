#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_write_report.py — Generate manuscript-ready markdown Methods/Results/Limitations

Reads summary tables produced by the het/ROH workflow and auto-fills a
markdown template with actual numbers. Produces both:
  - het_roh_analysis.md        (template with placeholders)
  - het_roh_analysis_filled.md (filled with data)

Usage:
  python3 06_write_report.py \
    --tables-dir 09_final_tables \
    --stats-dir  08_stats \
    --out-dir    10_report
"""

import argparse
import csv
import os
import sys
from pathlib import Path


def read_tsv_dict(path, key_col=None):
    """Read TSV and return list of row dicts."""
    rows = []
    if not os.path.isfile(path):
        return rows
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def safe_float(val, fmt=".6f"):
    try:
        return format(float(val), fmt)
    except (ValueError, TypeError):
        return "NA"


def safe_int(val):
    try:
        return str(int(float(val)))
    except (ValueError, TypeError):
        return "NA"


def compute_stats_from_master(master_rows):
    """Compute descriptive stats from master summary."""
    stats = {}
    if not master_rows:
        return stats

    n = len(master_rows)
    stats["n_samples"] = str(n)

    # Het
    hets = [float(r["het_genomewide"]) for r in master_rows
            if r.get("het_genomewide") not in (None, "", "NA")]
    if hets:
        stats["mean_het"] = f"{sum(hets)/len(hets):.6e}"
        stats["median_het"] = f"{sorted(hets)[len(hets)//2]:.6e}"
        stats["min_het"] = f"{min(hets):.6e}"
        stats["max_het"] = f"{max(hets):.6e}"
        stats["range_het"] = f"{min(hets):.6e}–{max(hets):.6e}"

    # FROH
    frohs = [float(r["froh"]) for r in master_rows
             if r.get("froh") not in (None, "", "NA")]
    if frohs:
        stats["mean_froh"] = f"{sum(frohs)/len(frohs):.4f}"
        stats["median_froh"] = f"{sorted(frohs)[len(frohs)//2]:.4f}"
        stats["min_froh"] = f"{min(frohs):.4f}"
        stats["max_froh"] = f"{max(frohs):.4f}"
        stats["range_froh"] = f"{min(frohs):.4f}–{max(frohs):.4f}"

    # ROH tracts
    ntracts = [int(float(r["n_tracts"])) for r in master_rows
               if r.get("n_tracts") not in (None, "", "NA")]
    if ntracts:
        stats["total_roh_tracts"] = str(sum(ntracts))
        stats["mean_n_tracts"] = f"{sum(ntracts)/len(ntracts):.1f}"
        stats["max_n_tracts"] = str(max(ntracts))

    # Longest ROH
    longest = [int(float(r["longest_roh"])) for r in master_rows
               if r.get("longest_roh") not in (None, "", "NA")]
    if longest:
        stats["max_longest_roh_mb"] = f"{max(longest)/1e6:.2f}"
        stats["mean_longest_roh_mb"] = f"{sum(longest)/len(longest)/1e6:.2f}"

    # Callable
    cbp = master_rows[0].get("callable_bp", "NA") if master_rows else "NA"
    try:
        stats["callable_bp"] = f"{int(float(cbp)):,}"
        stats["callable_mb"] = f"{int(float(cbp))/1e6:.1f}"
    except (ValueError, TypeError):
        stats["callable_bp"] = cbp
        stats["callable_mb"] = "NA"

    return stats


# ═══════════════════════════════════════════════════════════════════════════
# TEMPLATE
# ═══════════════════════════════════════════════════════════════════════════
TEMPLATE = r"""# Heterozygosity and Runs of Homozygosity Analysis

## Methods

### Genome-wide heterozygosity

For genome-wide heterozygosity, we ran ANGSD per sample (n = {n_samples}) using the
Samtools genotype-likelihood model (`-GL 1`) and site allele frequency likelihood
estimation (`-doSaf 1`). We estimated a per-sample global SFS with `realSFS`
(convergence: maxIter = 2000, tole = 1e-16) and used the resulting `est.ml` file
as the prior for local heterozygosity inference.

Callable sites were restricted to uppercase A/C/G/T positions only, totaling
{callable_mb} Mb ({callable_bp} bp) of callable non-repetitive genome. Soft-masked
lowercase bases were excluded as potentially repetitive, low-confidence, or
paralogy-prone sequence. Hard-masked or unknown bases (N) were excluded as
non-informative sequence. The callable genome was therefore defined as the
complement of the union of soft-masked and hard-masked sequence, restricted to
uppercase A/C/G/T bases.

Per-sample heterozygosity estimation was performed independently and therefore
did not require prior culling of related individuals.

### Local diversity tracks

For local/per-site diversity estimation, we reran ANGSD with the matching
per-sample global prior (`-pest est.ml`) to obtain site-level posterior estimates
via `-doThetas 1`. Windowed summaries were computed with `thetaStat do_stat`
using 500 kb non-overlapping windows.

**Important**: Local theta tracks are best treated as local diversity /
heterozygosity proxies. They are not exactly literal per-site observed
heterozygosity.

### Runs of homozygosity

ROH tracts were inferred with ngsF-HMM using genotype-likelihood input (BEAGLE
format from ANGSD `-doGlf 2`) and site-position files matched in identical
genomic order. We ran {n_ngsf_reps} random-start replicates and retained the
solution with the highest log-likelihood.

We converted contiguous runs assigned to the HMM IBD state into per-sample ROH
interval tables using `convert_ibd.pl`.

ROH lengths were reported as genomic interval spans on the reference assembly,
not as direct homozygosity observations for every base within the interval.

### ROH-associated diversity

ROH-associated heterozygosity was calculated by intersecting local theta
(diversity proxy) estimates with ROH intervals and summarizing mean values
inside versus outside ROH.

## Results

### Genome-wide heterozygosity

We assessed genome-wide heterozygosity for {n_samples} QC-passing samples using
a genotype-likelihood framework in ANGSD rather than hard genotype calls.
Mean genome-wide heterozygosity was {mean_het} (range: {range_het};
median: {median_het}).

### Runs of homozygosity

We inferred runs of homozygosity as contiguous autozygous tracts using a
two-state HMM on informative non-repetitive callable sites. Across all
{n_samples} samples, a total of {total_roh_tracts} ROH tracts were identified.
Mean F_ROH was {mean_froh} (range: {range_froh}; median: {median_froh}).
The longest individual ROH tract was {max_longest_roh_mb} Mb; the mean per-sample
longest ROH was {mean_longest_roh_mb} Mb.

We partitioned heterozygosity into ROH and non-ROH compartments to distinguish
recent autozygosity from background genome-wide diversity. ROH tracts were
summarized as physical interval spans on the reference genome, while inference
was driven by informative callable sites within those intervals.

{chr_summary_text}

### Summary statistics

| Metric | Value |
|--------|-------|
| N samples | {n_samples} |
| Callable genome (Mb) | {callable_mb} |
| Mean heterozygosity | {mean_het} |
| Mean F_ROH | {mean_froh} |
| Median F_ROH | {median_froh} |
| Range F_ROH | {range_froh} |
| Total ROH tracts | {total_roh_tracts} |
| Max longest ROH (Mb) | {max_longest_roh_mb} |

## Limitations

- Genome-wide and local heterozygosity were estimated on a conservative callable
  genome and therefore may underestimate diversity in repetitive or
  difficult-to-map regions.

- Because soft-masked and hard-masked sequence was excluded, some true
  heterozygous sites and some true ROH boundaries in masked regions may have
  been missed.

- Reported ROH spans reflect interval-level inference from informative callable
  sites and do not imply direct homozygosity evidence for every base within
  each tract.

- Repetitive and low-confidence regions inside inferred ROH were treated as
  uninformative rather than as direct evidence of homozygosity.

- Although relatedness does not invalidate per-sample heterozygosity estimation,
  strong family structure may still affect downstream interpretation of
  sample-level summaries across the cohort.

- ROH inference from site-based HMM states depends on marker density and
  callable-site distribution, so tract boundary precision may be reduced in
  sparsely informative regions.

- Local theta tracks should be interpreted as local diversity / heterozygosity
  proxies, not exact per-site observed heterozygosity.

- Whole-population Hardy-Weinberg-derived Hobs was not used as the primary
  heterozygosity track because hatchery structure can confound HWE-based
  summaries.

## Quick Reference Checklist

- Genome-wide het: GL1 + doSaf1 + realSFS
- Local het: rerun with matching -pest est.ml
- Callable: uppercase ATGC only
- Exclude: lowercase masked bases + N
- ROH: ngsF-HMM on informative callable sites
- ROH table: contiguous state-1 runs
- Het in/out ROH: intersect local theta proxy with ROH intervals

## Figure and Table Inventory

{figure_table_inventory}
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--stats-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--ngsf-reps", type=int, default=10)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ── Read master summary ──────────────────────────────────────────────
    master_path = os.path.join(args.tables_dir, "master_summary.tsv")
    master_rows = read_tsv_dict(master_path)
    stats = compute_stats_from_master(master_rows)
    stats["n_ngsf_reps"] = str(args.ngsf_reps)

    # ── Per-chromosome summary text ──────────────────────────────────────
    chr_path = os.path.join(args.tables_dir, "per_chr_roh_summary.tsv")
    chr_rows = read_tsv_dict(chr_path)
    chr_text = ""
    if chr_rows:
        # Summarize: which chromosomes have highest mean FROH
        from collections import defaultdict
        chr_frohs = defaultdict(list)
        for r in chr_rows:
            try:
                chr_frohs[r["chrom"]].append(float(r["froh_chr"]))
            except (ValueError, KeyError):
                pass
        if chr_frohs:
            chr_means = [(c, sum(v)/len(v)) for c, v in chr_frohs.items()]
            chr_means.sort(key=lambda x: -x[1])
            top3 = chr_means[:3]
            top3_str = ", ".join(f"{c} (mean F_ROH = {v:.4f})" for c, v in top3)
            chr_text = f"The chromosomes with the highest mean F_ROH were {top3_str}."

    stats["chr_summary_text"] = chr_text

    # ── Figure/table inventory ───────────────────────────────────────────
    inventory_lines = []

    # Check for tables
    table_files = [
        "master_summary.tsv",
        "per_sample_roh.tsv",
        "per_chr_roh_summary.tsv",
        "per_sample_het_in_out_roh.tsv",
        "roh_tracts_all.bed",
        "catfish_roh.per_sample_roh.tsv",
        "catfish_roh.per_sample_roh_bins_long.tsv",
    ]
    inventory_lines.append("### Tables\n")
    for tf in table_files:
        path = os.path.join(args.tables_dir, tf)
        status = "present" if os.path.isfile(path) else "not found"
        inventory_lines.append(f"- `{tf}`: {status}")

    # Check for stats
    stat_files = ["spearman_correlations.tsv", "group_comparisons.tsv", "descriptive_summary.tsv"]
    inventory_lines.append("\n### Statistical tests\n")
    for sf in stat_files:
        path = os.path.join(args.stats_dir, sf)
        status = "present" if os.path.isfile(path) else "not found"
        inventory_lines.append(f"- `{sf}`: {status}")

    inventory_lines.append("\n### Plots\n")
    inventory_lines.append("- See `06_plots_core/` and `07_plots_metadata/` for full figure set")

    stats["figure_table_inventory"] = "\n".join(inventory_lines)

    # ── Fill template ────────────────────────────────────────────────────
    # Write template (with placeholders)
    template_path = os.path.join(args.out_dir, "het_roh_analysis.md")
    with open(template_path, "w") as f:
        f.write(TEMPLATE)
    print(f"Wrote template: {template_path}")

    # Fill placeholders
    filled = TEMPLATE
    for key, val in stats.items():
        filled = filled.replace("{" + key + "}", val)

    # Replace any unfilled placeholders with "N/A"
    import re
    filled = re.sub(r'\{[a-z_]+\}', 'N/A', filled)

    filled_path = os.path.join(args.out_dir, "het_roh_analysis_filled.md")
    with open(filled_path, "w") as f:
        f.write(filled)
    print(f"Wrote filled report: {filled_path}")

    # ── Figures/tables used inventory ────────────────────────────────────
    tables_used_path = os.path.join(args.out_dir, "tables_used.tsv")
    with open(tables_used_path, "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["file", "directory", "status"])
        for tf in table_files:
            path = os.path.join(args.tables_dir, tf)
            w.writerow([tf, args.tables_dir, "present" if os.path.isfile(path) else "missing"])
        for sf in stat_files:
            path = os.path.join(args.stats_dir, sf)
            w.writerow([sf, args.stats_dir, "present" if os.path.isfile(path) else "missing"])
    print(f"Wrote: {tables_used_path}")


if __name__ == "__main__":
    main()
