#!/usr/bin/env python3
"""
07_build_breeding_tables.py — Build multiscale breeding-ready derived tables.

Reads existing workflow outputs and produces:
  A. per_sample_roh_bins_wide_{fixedBins,adaptiveBins}.tsv
  B. per_sample_genome_roh_segments.tsv
  C. per_chr_roh_group_summary.tsv
  D. recurrence_roh_windows_{100kb..1Mb}.tsv
  E. per_sample_window_theta_{scale}.tsv
  F. per_chr_theta_summary_{scale}.tsv
  G. per_chr_segment_theta_summary_{seg}_{scale}.tsv
  H. per_sample_chr_breeding_features_{scale}.tsv
  I. per_sample_segment_breeding_features_{seg}_{scale}.tsv
  J. pairwise_segment_complementarity_{seg}_{scale}.tsv
  K. pairwise_chr_complementarity_{seg}_{scale}.tsv
  L. pairwise_genome_complementarity_{seg}_{scale}.tsv
  + adaptive bin derivation + empirical ROH length distribution TSV

Usage:
  python3 07_build_breeding_tables.py \\
    --tracts-bed roh_tracts_all.bed \\
    --callable-bed callable.bed \\
    --chrom-sizes chrom_sizes.tsv \\
    --theta-dir 02_heterozygosity/03_theta/multiscale \\
    --sample-list samples_qcpass.txt \\
    --out-dir 11_breeding_tables \\
    [--ancestry-labels ancestry.tsv] \\
    [--order-file covtree_order.txt] \\
    [--pruned81 pruned81.txt] \\
    [--max-pairs 0]
"""

import argparse, csv, os, sys, math, itertools
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

# ═══════════════════════════════════════════════════════════════════════════
# ROH BIN DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════════
FIXED_BINS = [
    ("1-2Mb",   1_000_000,  2_000_000),
    ("2-4Mb",   2_000_000,  4_000_000),
    ("4-8Mb",   4_000_000,  8_000_000),
    ("8-16Mb",  8_000_000, 16_000_000),
    ("gt16Mb", 16_000_000,        None),
]

def assign_bin(length, bins):
    for name, lo, hi in bins:
        if length >= lo and (hi is None or length < hi):
            return name
    return None

def derive_adaptive_bins(all_lengths, out_dir):
    """Derive data-adapted ROH bins from empirical distribution."""
    if not all_lengths:
        return FIXED_BINS  # fallback
    lengths_mb = sorted([l / 1e6 for l in all_lengths if l >= 1_000_000])
    if len(lengths_mb) < 10:
        return FIXED_BINS

    # Write empirical distribution
    with open(os.path.join(out_dir, "empirical_roh_length_distribution.tsv"), "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["length_bp", "length_Mb"])
        for l in sorted(all_lengths):
            if l >= 1_000_000:
                w.writerow([l, f"{l/1e6:.3f}"])

    # Compute quantiles to find natural breaks
    import statistics
    q25 = statistics.quantiles(lengths_mb, n=4)[0]
    q50 = statistics.median(lengths_mb)
    q75 = statistics.quantiles(lengths_mb, n=4)[2]
    q90 = sorted(lengths_mb)[int(0.9 * len(lengths_mb))]
    mx = max(lengths_mb)

    # Build adaptive bins: aim for 5-6 bins with meaningful boundaries
    # Round to nearest Mb for interpretability
    def r(v): return max(1, round(v))

    boundaries = sorted(set([1, r(q25), r(q50), r(q75), r(q90)]))
    # Remove duplicates and ensure ascending
    boundaries = sorted(set(boundaries))
    if boundaries[0] != 1:
        boundaries = [1] + boundaries

    bins = []
    for i in range(len(boundaries)):
        lo = boundaries[i]
        hi = boundaries[i+1] if i+1 < len(boundaries) else None
        if hi is not None:
            name = f"{lo}-{hi}Mb"
        else:
            name = f"gt{lo}Mb"
        bins.append((name, lo * 1_000_000, hi * 1_000_000 if hi else None))

    # Ensure last bin is unbounded
    if bins[-1][2] is not None:
        last_lo = bins[-1][1]
        bins[-1] = (f"gt{int(last_lo/1e6)}Mb", last_lo, None)

    # Write bin definitions
    with open(os.path.join(out_dir, "adaptive_bin_definitions.tsv"), "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["bin_label", "lo_bp", "hi_bp", "q25_Mb", "q50_Mb", "q75_Mb", "q90_Mb", "max_Mb"])
        for name, lo, hi in bins:
            w.writerow([name, lo, hi if hi else "Inf",
                        f"{q25:.2f}", f"{q50:.2f}", f"{q75:.2f}", f"{q90:.2f}", f"{mx:.2f}"])

    return bins

# ═══════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════
def load_tracts(path):
    """Load ROH tracts BED: chr start end sample [length]"""
    tracts = []
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            p = ln.strip().split()
            if len(p) < 4: continue
            c, s, e, samp = p[0], int(p[1]), int(p[2]), p[3]
            tracts.append((c, s, e, samp, e - s))
    return tracts

def load_chrom_sizes(path):
    sizes = {}
    with open(path) as f:
        for ln in f:
            p = ln.strip().split()
            if len(p) >= 2:
                sizes[p[0]] = int(p[1])
    return sizes

def load_callable_per_chr(path):
    d = defaultdict(int)
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            p = ln.strip().split()
            if len(p) >= 3:
                d[p[0]] += int(p[2]) - int(p[1])
    return dict(d)

def load_sample_list(path):
    with open(path) as f:
        return [ln.strip().split()[0] for ln in f if ln.strip()]

def load_order(path):
    if not path or not os.path.isfile(path): return {}
    with open(path) as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    return {s: i for i, s in enumerate(lines)}

def load_ancestry(path):
    if not path or not os.path.isfile(path): return {}
    d = {}
    with open(path) as f:
        for ln in f:
            p = ln.strip().split()
            if len(p) >= 2:
                d[p[0]] = p[1]
    return d

def load_pruned(path):
    if not path or not os.path.isfile(path): return set()
    with open(path) as f:
        return {ln.strip().split()[0] for ln in f if ln.strip()}

# ═══════════════════════════════════════════════════════════════════════════
# THETA LOADING
# ═══════════════════════════════════════════════════════════════════════════
def find_theta_files(theta_dir, samples, win, step):
    """Find thetaStat pestPG files for given win/step."""
    found = {}
    for samp in samples:
        pat = f"{samp}.win{win}.step{step}.pestPG"
        fp = os.path.join(theta_dir, pat)
        if os.path.isfile(fp):
            found[samp] = fp
    return found

def parse_theta_file(path):
    """Parse thetaStat do_stat output. Returns list of (chrom, center, tP, nSites)."""
    rows = []
    with open(path) as f:
        header = f.readline().strip().split()
        # Detect columns
        cn = [h.lower() for h in header]
        # thetaStat format varies; common: (indexStart indexStop) Chromo WinCenter tW tP ... nSites
        # Or: Chr WinCenter tW tP ...
        ci = next((i for i, h in enumerate(cn) if h in ("chr","chromo","chrom")), None)
        wi = next((i for i, h in enumerate(cn) if h in ("wincenter","midpos")), None)
        ti = next((i for i, h in enumerate(cn) if h in ("tp","thetapi","pi")), None)
        ni = next((i for i, h in enumerate(cn) if h.startswith("nsite")), None)

        if ci is None or wi is None or ti is None:
            # Positional fallback
            if len(header) >= 5:
                ci, wi, ti = 0, 1, 3
                ni = len(header) - 1

        if ci is None: return rows

        for ln in f:
            p = ln.strip().split()
            try:
                chrom = p[ci]
                center = int(float(p[wi]))
                tp = float(p[ti])
                ns = int(p[ni]) if ni is not None and ni < len(p) else 1
            except (ValueError, IndexError):
                continue
            if ns > 0:
                rows.append((chrom, center, tp, ns))
    return rows

# ═══════════════════════════════════════════════════════════════════════════
# TABLE BUILDERS
# ═══════════════════════════════════════════════════════════════════════════

def build_roh_bins_wide(tracts, callable_total, bins, suffix, out_dir):
    """Table A: per_sample_roh_bins_wide_{suffix}.tsv"""
    per_sample = defaultdict(lambda: defaultdict(int))
    per_sample_total = defaultdict(int)
    for c, s, e, samp, length in tracts:
        if length < 1_000_000: continue
        b = assign_bin(length, bins)
        if b:
            per_sample[samp][b] += length
        per_sample_total[samp] += length

    path = os.path.join(out_dir, f"per_sample_roh_bins_wide_{suffix}.tsv")
    bin_names = [b[0] for b in bins]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        header = ["sample"]
        for bn in bin_names:
            header += [f"roh_bp_{bn}", f"pct_genome_{bn}"]
        header += ["total_roh_bp", "froh"]
        w.writerow(header)
        for samp in sorted(per_sample.keys() | per_sample_total.keys()):
            row = [samp]
            for bn in bin_names:
                bp = per_sample[samp].get(bn, 0)
                pct = 100 * bp / callable_total if callable_total > 0 else 0
                row += [bp, f"{pct:.4f}"]
            tot = per_sample_total.get(samp, 0)
            froh = tot / callable_total if callable_total > 0 else 0
            row += [tot, f"{froh:.6f}"]
            w.writerow(row)
    return path

def build_genome_roh_segments(tracts, order, out_dir):
    """Table B: per_sample_genome_roh_segments.tsv"""
    path = os.path.join(out_dir, "per_sample_genome_roh_segments.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "chrom", "start", "end", "roh_length_bp",
                     "roh_bin_fixed", "order_index"])
        for c, s, e, samp, length in sorted(tracts, key=lambda x: (x[3], x[0], x[1])):
            if length < 1_000_000: continue
            bn = assign_bin(length, FIXED_BINS) or "lt1Mb"
            oi = order.get(samp, -1)
            w.writerow([samp, c, s, e, length, bn, oi])
    return path

def build_recurrence_windows(tracts, chrom_sizes, samples, win_bp, label, out_dir,
                             ancestry, pruned):
    """Table D: recurrence_roh_windows_{label}.tsv — multi-group"""
    # Build group assignments
    groups = {}
    for s in samples:
        groups[s] = ancestry.get(s, "all")

    # Build windows
    windows = []
    for chrom in sorted(chrom_sizes.keys()):
        clen = chrom_sizes[chrom]
        for start in range(0, clen, win_bp):
            end = min(start + win_bp, clen)
            windows.append((chrom, start, end))

    # For each window, count per-group how many samples have ROH overlapping
    # Index tracts by chrom for speed
    tracts_by_chr = defaultdict(list)
    for c, s, e, samp, length in tracts:
        tracts_by_chr[c].append((s, e, samp))

    group_names = sorted(set(groups.values()))
    group_samples = {g: [s for s in samples if groups.get(s) == g] for g in group_names}

    path = os.path.join(out_dir, f"recurrence_roh_windows_{label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        header = ["chrom", "start", "end"]
        for g in group_names:
            header += [f"n_with_roh_{g}", f"pct_with_roh_{g}", f"n_total_{g}"]
        header += ["n_with_roh_all", "pct_with_roh_all", "n_total_all"]
        w.writerow(header)

        for chrom, wstart, wend in windows:
            row = [chrom, wstart, wend]
            ct = tracts_by_chr.get(chrom, [])
            # Find samples with ROH overlapping this window
            samps_with_roh = set()
            for ts, te, tsamp in ct:
                if ts < wend and te > wstart:
                    samps_with_roh.add(tsamp)

            for g in group_names:
                gs = group_samples[g]
                n = sum(1 for s in gs if s in samps_with_roh)
                pct = 100 * n / len(gs) if gs else 0
                row += [n, f"{pct:.1f}", len(gs)]

            n_all = len(samps_with_roh & set(samples))
            pct_all = 100 * n_all / len(samples) if samples else 0
            row += [n_all, f"{pct_all:.1f}", len(samples)]
            w.writerow(row)
    return path

def build_per_chr_theta_summary(theta_data, tracts, scale_label, out_dir):
    """Table F: per_chr_theta_summary_{scale}.tsv"""
    # theta_data: dict sample -> list of (chrom, center, tP, nSites)
    # For ROH overlap: need tracts indexed
    tracts_by_sample_chr = defaultdict(list)
    for c, s, e, samp, length in tracts:
        tracts_by_sample_chr[(samp, c)].append((s, e))

    path = os.path.join(out_dir, f"per_chr_theta_summary_{scale_label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "chrom", "n_windows", "mean_theta", "median_theta",
                     "sd_theta", "min_theta", "max_theta",
                     "mean_theta_in_roh", "mean_theta_out_roh", "ratio_theta_in_out"])

        for samp in sorted(theta_data.keys()):
            # Group by chrom
            by_chr = defaultdict(list)
            for chrom, center, tp, ns in theta_data[samp]:
                tps = tp / ns if ns > 0 else 0
                by_chr[chrom].append((center, tps))

            for chrom in sorted(by_chr.keys()):
                vals = [v for _, v in by_chr[chrom]]
                centers = [c for c, _ in by_chr[chrom]]
                if not vals: continue
                import statistics
                mn = statistics.mean(vals)
                md = statistics.median(vals)
                sd = statistics.stdev(vals) if len(vals) > 1 else 0
                mi, mx = min(vals), max(vals)

                # Split in/out ROH
                roh_intervals = tracts_by_sample_chr.get((samp, chrom), [])
                in_vals, out_vals = [], []
                for center, tps in zip(centers, vals):
                    in_roh = any(rs <= center <= re for rs, re in roh_intervals)
                    if in_roh:
                        in_vals.append(tps)
                    else:
                        out_vals.append(tps)

                mn_in = statistics.mean(in_vals) if in_vals else None
                mn_out = statistics.mean(out_vals) if out_vals else None
                ratio = mn_in / mn_out if mn_in is not None and mn_out and mn_out > 0 else None

                w.writerow([samp, chrom, len(vals),
                            f"{mn:.8e}", f"{md:.8e}", f"{sd:.8e}", f"{mi:.8e}", f"{mx:.8e}",
                            f"{mn_in:.8e}" if mn_in is not None else "NA",
                            f"{mn_out:.8e}" if mn_out is not None else "NA",
                            f"{ratio:.6f}" if ratio is not None else "NA"])
    return path

def build_chr_breeding_features(theta_data, tracts, callable_per_chr, scale_label, out_dir):
    """Table H: per_sample_chr_breeding_features_{scale}.tsv"""
    tracts_by_sc = defaultdict(list)
    for c, s, e, samp, length in tracts:
        tracts_by_sc[(samp, c)].append((s, e, length))

    theta_by_sc = defaultdict(list)
    for samp, rows in theta_data.items():
        for chrom, center, tp, ns in rows:
            tps = tp / ns if ns > 0 else 0
            theta_by_sc[(samp, chrom)].append((center, tps))

    all_keys = set(tracts_by_sc.keys()) | set(theta_by_sc.keys())
    path = os.path.join(out_dir, f"per_sample_chr_breeding_features_{scale_label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "chrom", "callable_bp_chr",
                     "mean_theta_chr", "median_theta_chr", "frac_low_theta_windows_chr",
                     "roh_bp_chr", "froh_chr", "n_roh_chr", "longest_roh_chr",
                     "mean_theta_in_roh_chr", "mean_theta_out_roh_chr", "ratio_theta_in_out_chr"])

        import statistics
        for samp, chrom in sorted(all_keys):
            cbp = callable_per_chr.get(chrom, 0)
            roh_list = tracts_by_sc.get((samp, chrom), [])
            roh_bp = sum(l for _, _, l in roh_list)
            n_roh = len(roh_list)
            longest = max((l for _, _, l in roh_list), default=0)
            froh = roh_bp / cbp if cbp > 0 else 0

            theta_wins = theta_by_sc.get((samp, chrom), [])
            vals = [v for _, v in theta_wins]

            if vals:
                mn = statistics.mean(vals)
                md = statistics.median(vals)
                # Low theta = below 10th percentile of all sample's theta
                all_theta_vals = [v for _, v in theta_by_sc.get((samp, chrom), [])]
                threshold = sorted(all_theta_vals)[max(0, int(0.1 * len(all_theta_vals)))] if all_theta_vals else 0
                frac_low = sum(1 for v in vals if v <= threshold) / len(vals)
            else:
                mn = md = frac_low = None

            roh_intervals = [(s, e) for s, e, _ in roh_list]
            in_vals = [v for c, v in theta_wins if any(rs <= c <= re for rs, re in roh_intervals)]
            out_vals = [v for c, v in theta_wins if not any(rs <= c <= re for rs, re in roh_intervals)]
            mn_in = statistics.mean(in_vals) if in_vals else None
            mn_out = statistics.mean(out_vals) if out_vals else None
            ratio = mn_in / mn_out if mn_in is not None and mn_out and mn_out > 0 else None

            w.writerow([samp, chrom, cbp,
                        f"{mn:.8e}" if mn is not None else "NA",
                        f"{md:.8e}" if md is not None else "NA",
                        f"{frac_low:.4f}" if frac_low is not None else "NA",
                        roh_bp, f"{froh:.6f}", n_roh, longest,
                        f"{mn_in:.8e}" if mn_in is not None else "NA",
                        f"{mn_out:.8e}" if mn_out is not None else "NA",
                        f"{ratio:.6f}" if ratio is not None else "NA"])
    return path

def build_segment_breeding_features(theta_data, tracts, callable_per_chr, chrom_sizes,
                                     seg_bp, seg_label, scale_label, out_dir):
    """Table I: per_sample_segment_breeding_features_{seg}_{scale}.tsv"""
    # Build segments
    segments = []
    for chrom in sorted(chrom_sizes.keys()):
        clen = chrom_sizes[chrom]
        for i, start in enumerate(range(0, clen, seg_bp)):
            end = min(start + seg_bp, clen)
            sid = f"{chrom}_seg{i+1}"
            segments.append((chrom, start, end, sid))

    tracts_by_sc = defaultdict(list)
    for c, s, e, samp, length in tracts:
        tracts_by_sc[(samp, c)].append((s, e))

    theta_by_sc = defaultdict(list)
    for samp, rows in theta_data.items():
        for chrom, center, tp, ns in rows:
            tps = tp / ns if ns > 0 else 0
            theta_by_sc[(samp, chrom)].append((center, tps))

    samples = sorted(theta_data.keys() | {t[3] for t in tracts})
    path = os.path.join(out_dir, f"per_sample_segment_breeding_features_{seg_label}_{scale_label}.tsv")

    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "chrom", "segment_id", "start", "end",
                     "callable_bp_segment", "mean_theta_segment", "median_theta_segment",
                     "frac_low_theta_windows_segment", "roh_bp_segment", "froh_segment",
                     "theta_in_roh_segment", "theta_out_roh_segment", "ratio_theta_in_out_segment",
                     "low_div_flag", "high_roh_flag"])

        import statistics
        for samp in samples:
            for chrom, seg_s, seg_e, sid in segments:
                seg_len = seg_e - seg_s
                # ROH in segment
                roh_intervals = tracts_by_sc.get((samp, chrom), [])
                roh_bp = 0
                for rs, re in roh_intervals:
                    ov_s = max(rs, seg_s)
                    ov_e = min(re, seg_e)
                    if ov_e > ov_s:
                        roh_bp += ov_e - ov_s
                froh = roh_bp / seg_len if seg_len > 0 else 0

                # Theta in segment
                theta_wins = theta_by_sc.get((samp, chrom), [])
                seg_vals = [(c, v) for c, v in theta_wins if seg_s <= c < seg_e]
                vals = [v for _, v in seg_vals]

                if vals:
                    mn = statistics.mean(vals)
                    md = statistics.median(vals)
                    # Global threshold: 10th percentile of ALL this sample's theta
                    all_vals = []
                    for ch_key in theta_by_sc:
                        if ch_key[0] == samp:
                            all_vals.extend(v for _, v in theta_by_sc[ch_key])
                    thr = sorted(all_vals)[max(0, int(0.1 * len(all_vals)))] if all_vals else 0
                    frac_low = sum(1 for v in vals if v <= thr) / len(vals)
                else:
                    mn = md = frac_low = None

                # In/out ROH within segment
                in_v = [v for c, v in seg_vals if any(rs <= c <= re for rs, re in roh_intervals)]
                out_v = [v for c, v in seg_vals if not any(rs <= c <= re for rs, re in roh_intervals)]
                mn_in = statistics.mean(in_v) if in_v else None
                mn_out = statistics.mean(out_v) if out_v else None
                ratio = mn_in / mn_out if mn_in is not None and mn_out and mn_out > 0 else None

                low_div = 1 if frac_low is not None and frac_low > 0.5 else 0
                high_roh = 1 if froh > 0.5 else 0

                w.writerow([samp, chrom, sid, seg_s, seg_e,
                            seg_len,
                            f"{mn:.8e}" if mn is not None else "NA",
                            f"{md:.8e}" if md is not None else "NA",
                            f"{frac_low:.4f}" if frac_low is not None else "NA",
                            roh_bp, f"{froh:.6f}",
                            f"{mn_in:.8e}" if mn_in is not None else "NA",
                            f"{mn_out:.8e}" if mn_out is not None else "NA",
                            f"{ratio:.6f}" if ratio is not None else "NA",
                            low_div, high_roh])
    return path

def build_pairwise_complementarity(seg_features_path, seg_label, scale_label, out_dir, max_pairs=0):
    """Tables J, K, L: pairwise complementarity at segment, chr, genome level."""
    import csv
    # Load segment features
    rows_by_sample = defaultdict(list)
    with open(seg_features_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows_by_sample[row["sample"]].append(row)

    samples = sorted(rows_by_sample.keys())
    pairs = list(itertools.combinations(samples, 2))
    if max_pairs > 0 and len(pairs) > max_pairs:
        print(f"  Limiting to {max_pairs} pairs (of {len(pairs)} total)")
        pairs = pairs[:max_pairs]

    # Segment-level
    seg_path = os.path.join(out_dir, f"pairwise_segment_complementarity_{seg_label}_{scale_label}.tsv")
    chr_agg = defaultdict(lambda: {"n_seg": 0, "n_both_bad": 0, "n_comp": 0, "scores": []})
    genome_agg = defaultdict(lambda: {"n_chr": set(), "n_both_bad": 0, "n_comp": 0, "scores": []})

    with open(seg_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_A", "sample_B", "chrom", "segment_id", "start", "end",
                     "A_mean_theta_segment", "B_mean_theta_segment",
                     "A_roh_bp_segment", "B_roh_bp_segment",
                     "A_low_div_flag", "B_low_div_flag",
                     "A_high_roh_flag", "B_high_roh_flag",
                     "both_bad_flag", "complementarity_flag", "complementarity_score_simple"])

        for sa, sb in pairs:
            rows_a = {(r["chrom"], r["segment_id"]): r for r in rows_by_sample[sa]}
            rows_b = {(r["chrom"], r["segment_id"]): r for r in rows_by_sample[sb]}
            all_segs = sorted(set(rows_a.keys()) | set(rows_b.keys()))

            for key in all_segs:
                ra = rows_a.get(key, {})
                rb = rows_b.get(key, {})
                chrom = ra.get("chrom", rb.get("chrom", ""))
                sid = ra.get("segment_id", rb.get("segment_id", ""))
                start = ra.get("start", rb.get("start", ""))
                end = ra.get("end", rb.get("end", ""))

                a_theta = ra.get("mean_theta_segment", "NA")
                b_theta = rb.get("mean_theta_segment", "NA")
                a_roh = int(ra.get("roh_bp_segment", 0))
                b_roh = int(rb.get("roh_bp_segment", 0))
                a_low = int(ra.get("low_div_flag", 0))
                b_low = int(rb.get("low_div_flag", 0))
                a_high = int(ra.get("high_roh_flag", 0))
                b_high = int(rb.get("high_roh_flag", 0))

                a_bad = a_low or a_high
                b_bad = b_low or b_high
                both_bad = 1 if a_bad and b_bad else 0
                comp = 1 if (a_bad and not b_bad) or (b_bad and not a_bad) else 0

                # Simple score: theta difference (higher = more complementary)
                try:
                    at = float(a_theta); bt = float(b_theta)
                    score = abs(at - bt)
                except (ValueError, TypeError):
                    score = 0

                w.writerow([sa, sb, chrom, sid, start, end,
                            a_theta, b_theta, a_roh, b_roh,
                            a_low, b_low, a_high, b_high,
                            both_bad, comp, f"{score:.8e}"])

                # Aggregate
                chr_agg[(sa, sb, chrom)]["n_seg"] += 1
                chr_agg[(sa, sb, chrom)]["n_both_bad"] += both_bad
                chr_agg[(sa, sb, chrom)]["n_comp"] += comp
                chr_agg[(sa, sb, chrom)]["scores"].append(score)
                genome_agg[(sa, sb)]["n_chr"].add(chrom)
                genome_agg[(sa, sb)]["n_both_bad"] += both_bad
                genome_agg[(sa, sb)]["n_comp"] += comp
                genome_agg[(sa, sb)]["scores"].append(score)

    # Chr-level
    chr_path = os.path.join(out_dir, f"pairwise_chr_complementarity_{seg_label}_{scale_label}.tsv")
    with open(chr_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_A", "sample_B", "chrom", "n_segments",
                     "n_both_bad_segments", "n_complementary_segments",
                     "mean_complementarity_score", "chr_pair_risk_score", "chr_pair_rescue_score"])
        for (sa, sb, chrom), d in sorted(chr_agg.items()):
            mn_sc = sum(d["scores"]) / len(d["scores"]) if d["scores"] else 0
            risk = d["n_both_bad"] / d["n_seg"] if d["n_seg"] > 0 else 0
            rescue = d["n_comp"] / d["n_seg"] if d["n_seg"] > 0 else 0
            w.writerow([sa, sb, chrom, d["n_seg"], d["n_both_bad"], d["n_comp"],
                        f"{mn_sc:.8e}", f"{risk:.4f}", f"{rescue:.4f}"])

    # Genome-level
    gen_path = os.path.join(out_dir, f"pairwise_genome_complementarity_{seg_label}_{scale_label}.tsv")
    with open(gen_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_A", "sample_B", "n_chr",
                     "total_both_bad_segments", "total_complementary_segments",
                     "mean_genome_complementarity_score",
                     "genome_pair_risk_score", "genome_pair_rescue_score"])
        for (sa, sb), d in sorted(genome_agg.items()):
            n_chr = len(d["n_chr"])
            total_seg = len(d["scores"])
            mn_sc = sum(d["scores"]) / total_seg if total_seg > 0 else 0
            risk = d["n_both_bad"] / total_seg if total_seg > 0 else 0
            rescue = d["n_comp"] / total_seg if total_seg > 0 else 0
            w.writerow([sa, sb, n_chr, d["n_both_bad"], d["n_comp"],
                        f"{mn_sc:.8e}", f"{risk:.4f}", f"{rescue:.4f}"])

    return seg_path, chr_path, gen_path

# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracts-bed", required=True)
    ap.add_argument("--callable-bed", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--theta-dir", required=True, help="Directory with multiscale theta files")
    ap.add_argument("--theta-main-dir", default="", help="Directory with main 500kb theta files")
    ap.add_argument("--sample-list", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--ancestry-labels", default="")
    ap.add_argument("--order-file", default="")
    ap.add_argument("--pruned81", default="")
    ap.add_argument("--max-pairs", type=int, default=0, help="Limit pairwise comparisons (0=all)")

    # Scale definitions (can be overridden)
    ap.add_argument("--theta-scales", nargs="+", default=["5000_1000", "10000_2000", "50000_10000"])
    ap.add_argument("--theta-labels", nargs="+", default=["5kb_1kb", "10kb_2kb", "50kb_10kb"])
    ap.add_argument("--segment-sizes", nargs="+", type=int, default=[500000, 1000000, 2000000, 5000000])
    ap.add_argument("--segment-labels", nargs="+", default=["500kb", "1Mb", "2Mb", "5Mb"])
    ap.add_argument("--recurrence-sizes", nargs="+", type=int, default=[100000, 250000, 500000, 1000000])
    ap.add_argument("--recurrence-labels", nargs="+", default=["100kb", "250kb", "500kb", "1Mb"])

    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("=== Loading data ===")
    tracts = load_tracts(args.tracts_bed)
    chrom_sizes = load_chrom_sizes(args.chrom_sizes)
    callable_per_chr = load_callable_per_chr(args.callable_bed)
    callable_total = sum(callable_per_chr.values())
    samples = load_sample_list(args.sample_list)
    order = load_order(args.order_file)
    ancestry = load_ancestry(args.ancestry_labels)
    pruned = load_pruned(args.pruned81)

    print(f"  {len(tracts)} tracts, {len(samples)} samples, {len(chrom_sizes)} chroms")
    print(f"  Callable: {callable_total:,} bp")

    # ── A. ROH bins (dual) ───────────────────────────────────────────────
    print("\n=== A. ROH bins (fixed + adaptive) ===")
    all_lengths = [l for _, _, _, _, l in tracts if l >= 1_000_000]
    adaptive_bins = derive_adaptive_bins(all_lengths, args.out_dir)
    print(f"  Adaptive bins: {[b[0] for b in adaptive_bins]}")

    build_roh_bins_wide(tracts, callable_total, FIXED_BINS, "fixedBins", args.out_dir)
    build_roh_bins_wide(tracts, callable_total, adaptive_bins, "adaptiveBins", args.out_dir)

    # ── B. Genome ROH segments ───────────────────────────────────────────
    print("\n=== B. Genome ROH segments ===")
    build_genome_roh_segments(tracts, order, args.out_dir)

    # ── D. Recurrence windows ────────────────────────────────────────────
    print("\n=== D. Recurrence windows ===")
    for win_bp, label in zip(args.recurrence_sizes, args.recurrence_labels):
        print(f"  {label}...")
        build_recurrence_windows(tracts, chrom_sizes, samples, win_bp, label,
                                 args.out_dir, ancestry, pruned)

    # ── Theta-dependent tables (loop over scales) ────────────────────────
    for scale, scale_label in zip(args.theta_scales, args.theta_labels):
        win, step = scale.split("_")
        print(f"\n=== Theta scale: {scale_label} (win={win}, step={step}) ===")

        theta_files = find_theta_files(args.theta_dir, samples, win, step)
        if not theta_files:
            # Also try main dir
            if args.theta_main_dir:
                theta_files = find_theta_files(args.theta_main_dir, samples, win, step)
        if not theta_files:
            print(f"  No theta files found for {scale_label}, skipping")
            continue

        print(f"  Found {len(theta_files)} theta files")
        theta_data = {}
        for samp, fp in theta_files.items():
            theta_data[samp] = parse_theta_file(fp)

        # ── F. Per-chr theta summary ─────────────────────────────────────
        print(f"  F. per_chr_theta_summary_{scale_label}")
        build_per_chr_theta_summary(theta_data, tracts, scale_label, args.out_dir)

        # ── H. Chr breeding features ─────────────────────────────────────
        print(f"  H. per_sample_chr_breeding_features_{scale_label}")
        build_chr_breeding_features(theta_data, tracts, callable_per_chr, scale_label, args.out_dir)

        # ── Loop over segment sizes ──────────────────────────────────────
        for seg_bp, seg_label in zip(args.segment_sizes, args.segment_labels):
            print(f"  I. segment features {seg_label}_{scale_label}")
            seg_path = build_segment_breeding_features(
                theta_data, tracts, callable_per_chr, chrom_sizes,
                seg_bp, seg_label, scale_label, args.out_dir)

            # ── J,K,L. Pairwise complementarity ─────────────────────────
            if args.max_pairs != -1:  # -1 to skip entirely
                print(f"  J,K,L. pairwise complementarity {seg_label}_{scale_label}")
                build_pairwise_complementarity(
                    seg_path, seg_label, scale_label, args.out_dir,
                    max_pairs=args.max_pairs)

    print("\n=== Breeding tables complete ===")
    print(f"Output: {args.out_dir}")


if __name__ == "__main__":
    main()
