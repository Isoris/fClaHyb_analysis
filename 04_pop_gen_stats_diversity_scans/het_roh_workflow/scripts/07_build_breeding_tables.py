#!/usr/bin/env python3
# v2 note: added tables C/E/G; overlap-based ROH classification; per-chr shared
#          low-theta thresholds; pruned81 support; per_chr_theta_thresholds output
"""
07_build_breeding_tables.py — Build multiscale breeding-ready derived tables.

Outputs (all implemented):
  A. per_sample_roh_bins_wide_{fixedBins,adaptiveBins}.tsv
  B. per_sample_genome_roh_segments.tsv
  C. per_chr_roh_group_summary.tsv
  D. recurrence_roh_windows_{scale}.tsv
  E. per_sample_window_theta_{scale}.tsv
  F. per_chr_theta_summary_{scale}.tsv
  G. per_chr_segment_theta_summary_{seg}_{scale}.tsv
  H. per_sample_chr_breeding_features_{scale}.tsv
  I. per_sample_segment_breeding_features_{seg}_{scale}.tsv
  J. pairwise_segment_complementarity_{seg}_{scale}.tsv
  K. pairwise_chr_complementarity_{seg}_{scale}.tsv
  L. pairwise_genome_complementarity_{seg}_{scale}.tsv
  + adaptive_bin_definitions.tsv, empirical_roh_length_distribution.tsv
  + per_chr_theta_thresholds_{scale}.tsv
"""

import argparse, csv, os, sys, math, itertools, statistics
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Set

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
    if not all_lengths or len(all_lengths) < 10:
        return FIXED_BINS
    lengths_mb = sorted([l / 1e6 for l in all_lengths])
    with open(os.path.join(out_dir, "empirical_roh_length_distribution.tsv"), "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["length_bp", "length_Mb"])
        for l in sorted(all_lengths):
            w.writerow([l, f"{l/1e6:.3f}"])
    q25 = statistics.quantiles(lengths_mb, n=4)[0]
    q50 = statistics.median(lengths_mb)
    q75 = statistics.quantiles(lengths_mb, n=4)[2]
    q90 = sorted(lengths_mb)[int(0.9 * len(lengths_mb))]
    mx = max(lengths_mb)
    def r(v): return max(1, round(v))
    boundaries = sorted(set([1, r(q25), r(q50), r(q75), r(q90)]))
    if boundaries[0] != 1: boundaries = [1] + boundaries
    bins = []
    for i in range(len(boundaries)):
        lo = boundaries[i]
        hi = boundaries[i+1] if i+1 < len(boundaries) else None
        name = f"{lo}-{hi}Mb" if hi else f"gt{lo}Mb"
        bins.append((name, lo*1_000_000, hi*1_000_000 if hi else None))
    if bins[-1][2] is not None:
        last_lo = bins[-1][1]
        bins[-1] = (f"gt{int(last_lo/1e6)}Mb", last_lo, None)
    with open(os.path.join(out_dir, "adaptive_bin_definitions.tsv"), "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["bin_label","lo_bp","hi_bp","q25_Mb","q50_Mb","q75_Mb","q90_Mb","max_Mb"])
        for name, lo, hi in bins:
            w.writerow([name, lo, hi or "Inf", f"{q25:.2f}", f"{q50:.2f}", f"{q75:.2f}", f"{q90:.2f}", f"{mx:.2f}"])
    return bins

# ═══════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════
def load_tracts(path):
    tracts = []
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            p = ln.strip().split()
            if len(p) < 4: continue
            tracts.append((p[0], int(p[1]), int(p[2]), p[3], int(p[2])-int(p[1])))
    return tracts

def load_chrom_sizes(path):
    d = {}
    with open(path) as f:
        for ln in f:
            p = ln.strip().split()
            if len(p) >= 2: d[p[0]] = int(p[1])
    return d

def load_callable_per_chr(path):
    d = defaultdict(int)
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            p = ln.strip().split()
            if len(p) >= 3: d[p[0]] += int(p[2]) - int(p[1])
    return dict(d)

def load_list(path):
    if not path or not os.path.isfile(path): return []
    with open(path) as f:
        return [ln.strip().split()[0] for ln in f if ln.strip()]

def load_map(path):
    if not path or not os.path.isfile(path): return {}
    d = {}
    with open(path) as f:
        for i, ln in enumerate([x.strip() for x in f if x.strip()]):
            d[ln.split()[0]] = i
    return d

def load_ancestry(path):
    if not path or not os.path.isfile(path): return {}
    d = {}
    with open(path) as f:
        for ln in f:
            p = ln.strip().split()
            if len(p) >= 2: d[p[0]] = p[1]
    return d

# ═══════════════════════════════════════════════════════════════════════════
# THETA LOADING + ROH OVERLAP HELPERS
# ═══════════════════════════════════════════════════════════════════════════
def find_theta_files(theta_dir, samples, win, step):
    found = {}
    for samp in samples:
        fp = os.path.join(theta_dir, f"{samp}.win{win}.step{step}.pestPG")
        if os.path.isfile(fp): found[samp] = fp
    return found

def parse_theta_file(path):
    """Returns list of (chrom, center, tP_sum, nSites)."""
    rows = []
    with open(path) as f:
        header = f.readline().strip().split()
        cn = [h.lower() for h in header]
        ci = next((i for i,h in enumerate(cn) if h in ("chr","chromo","chrom")), None)
        wi = next((i for i,h in enumerate(cn) if h in ("wincenter","midpos")), None)
        ti = next((i for i,h in enumerate(cn) if h in ("tp","thetapi","pi")), None)
        ni = next((i for i,h in enumerate(cn) if h.startswith("nsite")), None)
        if ci is None and len(header) >= 5: ci, wi, ti = 0, 1, 3; ni = len(header)-1
        if ci is None: return rows
        for ln in f:
            p = ln.strip().split()
            try:
                chrom = p[ci]; center = int(float(p[wi]))
                tp = float(p[ti]); ns = int(p[ni]) if ni is not None and ni < len(p) else 1
            except (ValueError, IndexError): continue
            if ns > 0: rows.append((chrom, center, tp, ns))
    return rows

def compute_roh_overlap(win_start, win_end, roh_intervals):
    """Compute fraction of window [win_start, win_end) overlapping ROH."""
    win_len = win_end - win_start
    if win_len <= 0: return 0.0
    overlap = 0
    for rs, re in roh_intervals:
        ov = min(win_end, re) - max(win_start, rs)
        if ov > 0: overlap += ov
    return min(overlap / win_len, 1.0)

def index_tracts_by_sample_chr(tracts):
    d = defaultdict(list)
    for c, s, e, samp, length in tracts:
        d[(samp, c)].append((s, e))
    return d

# ═══════════════════════════════════════════════════════════════════════════
# CHROMOSOME-SHARED LOW-THETA THRESHOLDS
# ═══════════════════════════════════════════════════════════════════════════
def compute_chr_theta_thresholds(all_theta_data, win_size, scale_label, out_dir):
    """Compute 10th-percentile theta threshold per chromosome across all samples.
    Returns: dict chrom -> threshold_value"""
    chr_vals = defaultdict(list)
    half_w = win_size // 2
    for samp, rows in all_theta_data.items():
        for chrom, center, tp, ns in rows:
            tps = tp / ns if ns > 0 else 0
            chr_vals[chrom].append(tps)

    thresholds = {}
    path = os.path.join(out_dir, f"per_chr_theta_thresholds_{scale_label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "scale", "low_theta_threshold_p10", "n_windows_total"])
        for chrom in sorted(chr_vals.keys()):
            vals = sorted(chr_vals[chrom])
            idx = max(0, int(0.1 * len(vals)))
            thr = vals[idx] if vals else 0
            thresholds[chrom] = thr
            w.writerow([chrom, scale_label, f"{thr:.8e}", len(vals)])
    return thresholds

# ═══════════════════════════════════════════════════════════════════════════
# TABLE BUILDERS
# ═══════════════════════════════════════════════════════════════════════════

def build_roh_bins_wide(tracts, callable_total, bins, suffix, out_dir):
    per_sample = defaultdict(lambda: defaultdict(int))
    per_sample_total = defaultdict(int)
    for c, s, e, samp, length in tracts:
        if length < 1_000_000: continue
        b = assign_bin(length, bins)
        if b: per_sample[samp][b] += length
        per_sample_total[samp] += length
    bin_names = [b[0] for b in bins]
    path = os.path.join(out_dir, f"per_sample_roh_bins_wide_{suffix}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        header = ["sample"]
        for bn in bin_names: header += [f"roh_bp_{bn}", f"pct_genome_{bn}"]
        header += ["total_roh_bp", "froh"]
        w.writerow(header)
        for samp in sorted(per_sample.keys() | per_sample_total.keys()):
            row = [samp]
            for bn in bin_names:
                bp = per_sample[samp].get(bn, 0)
                pct = 100*bp/callable_total if callable_total > 0 else 0
                row += [bp, f"{pct:.4f}"]
            tot = per_sample_total.get(samp, 0)
            row += [tot, f"{tot/callable_total:.6f}" if callable_total > 0 else "NA"]
            w.writerow(row)

def build_genome_roh_segments(tracts, order, out_dir):
    path = os.path.join(out_dir, "per_sample_genome_roh_segments.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample","chrom","start","end","roh_length_bp","roh_bin_fixed","order_index"])
        for c, s, e, samp, length in sorted(tracts, key=lambda x: (x[3], x[0], x[1])):
            if length < 1_000_000: continue
            bn = assign_bin(length, FIXED_BINS) or "lt1Mb"
            w.writerow([samp, c, s, e, length, bn, order.get(samp, -1)])

def build_chr_roh_group_summary(tracts, callable_per_chr, samples, ancestry, pruned, out_dir):
    """Table C: per_chr_roh_group_summary.tsv"""
    groups = {}
    for s in samples:
        groups[s] = ancestry.get(s, "all")
    # Add pruned81 group
    if pruned:
        for s in samples:
            if s in pruned: groups[s + "__pruned81"] = "pruned81"

    # Per sample-chr ROH
    sc_roh = defaultdict(int)
    for c, s, e, samp, length in tracts:
        sc_roh[(samp, c)] += length

    group_names = sorted(set(groups.values()))
    all_chroms = sorted(callable_per_chr.keys())

    path = os.path.join(out_dir, "per_chr_roh_group_summary.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["group","chrom","mean_froh_chr","median_froh_chr","mean_roh_bp_chr",
                     "frac_samples_with_roh_chr","n_samples"])
        for grp in ["all"] + [g for g in group_names if g != "all"]:
            if grp == "pruned81":
                grp_samples = [s for s in samples if s in pruned]
            elif grp == "all":
                grp_samples = list(samples)
            else:
                grp_samples = [s for s in samples if ancestry.get(s) == grp]
            if not grp_samples: continue
            for chrom in all_chroms:
                cbp = callable_per_chr.get(chrom, 0)
                frohs = []
                roh_bps = []
                n_with = 0
                for s in grp_samples:
                    rb = sc_roh.get((s, chrom), 0)
                    roh_bps.append(rb)
                    frohs.append(rb / cbp if cbp > 0 else 0)
                    if rb > 0: n_with += 1
                w.writerow([grp, chrom,
                    f"{statistics.mean(frohs):.6f}" if frohs else "NA",
                    f"{statistics.median(frohs):.6f}" if frohs else "NA",
                    f"{statistics.mean(roh_bps):.0f}" if roh_bps else "NA",
                    f"{n_with/len(grp_samples):.4f}" if grp_samples else "NA",
                    len(grp_samples)])

def build_recurrence_windows(tracts, chrom_sizes, samples, win_bp, label, out_dir, ancestry, pruned):
    groups = {}
    for s in samples: groups[s] = ancestry.get(s, "all")
    tracts_by_chr = defaultdict(list)
    for c, s, e, samp, length in tracts:
        tracts_by_chr[c].append((s, e, samp))
    group_names = sorted(set(groups.values()))
    group_samples = {g: [s for s in samples if groups.get(s) == g] for g in group_names}
    # Add pruned81
    if pruned:
        group_names.append("pruned81")
        group_samples["pruned81"] = [s for s in samples if s in pruned]

    path = os.path.join(out_dir, f"recurrence_roh_windows_{label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        header = ["chrom","start","end"]
        for g in group_names: header += [f"n_with_roh_{g}", f"pct_with_roh_{g}", f"n_total_{g}"]
        header += ["n_with_roh_all","pct_with_roh_all","n_total_all"]
        w.writerow(header)
        for chrom in sorted(chrom_sizes.keys()):
            clen = chrom_sizes[chrom]
            ct = tracts_by_chr.get(chrom, [])
            for start in range(0, clen, win_bp):
                end = min(start + win_bp, clen)
                samps_with_roh = {ts for ts, te, tsamp in ct if ts < end and te > start for ts in [tsamp]}
                # fix: the comprehension above is wrong, redo properly
                samps_with = set()
                for ts, te, tsamp in ct:
                    if ts < end and te > start: samps_with.add(tsamp)
                row = [chrom, start, end]
                for g in group_names:
                    gs = group_samples[g]
                    n = sum(1 for s in gs if s in samps_with)
                    pct = 100*n/len(gs) if gs else 0
                    row += [n, f"{pct:.1f}", len(gs)]
                n_all = len(samps_with & set(samples))
                row += [n_all, f"{100*n_all/len(samples):.1f}" if samples else "0", len(samples)]
                w.writerow(row)

def build_per_sample_window_theta(theta_data, tracts_idx, win_size, scale_label,
                                   out_dir, order, ancestry, pruned):
    """Table E: per_sample_window_theta_{scale}.tsv"""
    half_w = win_size // 2
    path = os.path.join(out_dir, f"per_sample_window_theta_{scale_label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample","chrom","window_start","window_end","window_center",
                     "theta_pi_per_site","n_sites","overlap_roh_frac","in_roh_majority",
                     "order_index","ancestry","is_pruned81"])
        for samp in sorted(theta_data.keys()):
            oi = order.get(samp, -1)
            anc = ancestry.get(samp, "NA")
            pr = 1 if samp in pruned else 0
            for chrom, center, tp, ns in theta_data[samp]:
                tps = tp / ns if ns > 0 else 0
                ws = max(0, center - half_w)
                we = center + half_w
                roh_ivs = tracts_idx.get((samp, chrom), [])
                ov_frac = compute_roh_overlap(ws, we, roh_ivs)
                in_roh = 1 if ov_frac > 0.5 else 0
                w.writerow([samp, chrom, ws, we, center,
                            f"{tps:.8e}", ns, f"{ov_frac:.4f}", in_roh, oi, anc, pr])

def build_per_chr_theta_summary(theta_data, tracts_idx, win_size, thresholds,
                                 scale_label, out_dir):
    """Table F: per_chr_theta_summary_{scale}.tsv"""
    half_w = win_size // 2
    path = os.path.join(out_dir, f"per_chr_theta_summary_{scale_label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample","chrom","n_windows","mean_theta","median_theta",
                     "sd_theta","min_theta","max_theta",
                     "mean_theta_in_roh","mean_theta_out_roh","ratio_theta_in_out",
                     "frac_low_theta_windows"])
        for samp in sorted(theta_data.keys()):
            by_chr = defaultdict(list)
            for chrom, center, tp, ns in theta_data[samp]:
                tps = tp / ns if ns > 0 else 0
                ws = max(0, center - half_w); we = center + half_w
                roh_ivs = tracts_idx.get((samp, chrom), [])
                ov = compute_roh_overlap(ws, we, roh_ivs)
                by_chr[chrom].append((tps, ov))
            for chrom in sorted(by_chr.keys()):
                vals = [v for v, _ in by_chr[chrom]]
                ovs = [o for _, o in by_chr[chrom]]
                if not vals: continue
                mn = statistics.mean(vals); md = statistics.median(vals)
                sd = statistics.stdev(vals) if len(vals) > 1 else 0
                in_v = [v for v, o in zip(vals, ovs) if o > 0.5]
                out_v = [v for v, o in zip(vals, ovs) if o <= 0.5]
                mn_in = statistics.mean(in_v) if in_v else None
                mn_out = statistics.mean(out_v) if out_v else None
                ratio = mn_in/mn_out if mn_in is not None and mn_out and mn_out > 0 else None
                thr = thresholds.get(chrom, 0)
                frac_low = sum(1 for v in vals if v <= thr) / len(vals)
                w.writerow([samp, chrom, len(vals),
                    f"{mn:.8e}", f"{md:.8e}", f"{sd:.8e}", f"{min(vals):.8e}", f"{max(vals):.8e}",
                    f"{mn_in:.8e}" if mn_in is not None else "NA",
                    f"{mn_out:.8e}" if mn_out is not None else "NA",
                    f"{ratio:.6f}" if ratio is not None else "NA",
                    f"{frac_low:.4f}"])

def build_chr_segment_theta_summary(theta_data, tracts_idx, chrom_sizes, win_size,
                                     thresholds, seg_bp, seg_label, scale_label, out_dir):
    """Table G: per_chr_segment_theta_summary_{seg}_{scale}.tsv"""
    half_w = win_size // 2
    segments = []
    for chrom in sorted(chrom_sizes.keys()):
        for i, start in enumerate(range(0, chrom_sizes[chrom], seg_bp)):
            segments.append((chrom, start, min(start+seg_bp, chrom_sizes[chrom]), f"{chrom}_seg{i+1}"))

    path = os.path.join(out_dir, f"per_chr_segment_theta_summary_{seg_label}_{scale_label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample","chrom","segment_id","start","end","n_windows",
                     "mean_theta","median_theta","sd_theta","min_theta","max_theta",
                     "roh_bp_segment","froh_segment","frac_low_theta_windows",
                     "mean_theta_in_roh","mean_theta_out_roh","ratio_theta_in_out"])
        for samp in sorted(theta_data.keys()):
            # Index theta by chrom
            by_chr = defaultdict(list)
            for chrom, center, tp, ns in theta_data[samp]:
                tps = tp / ns if ns > 0 else 0
                by_chr[chrom].append((center, tps))
            for chrom, seg_s, seg_e, sid in segments:
                wins = [(c, v) for c, v in by_chr.get(chrom, []) if seg_s <= c < seg_e]
                vals = [v for _, v in wins]
                # ROH in segment
                roh_ivs = tracts_idx.get((samp, chrom), [])
                roh_bp = 0
                for rs, re in roh_ivs:
                    ov = min(seg_e, re) - max(seg_s, rs)
                    if ov > 0: roh_bp += ov
                seg_len = seg_e - seg_s
                froh = roh_bp / seg_len if seg_len > 0 else 0
                if not vals:
                    w.writerow([samp, chrom, sid, seg_s, seg_e, 0,
                        "NA","NA","NA","NA","NA", roh_bp, f"{froh:.6f}",
                        "NA","NA","NA","NA"])
                    continue
                mn = statistics.mean(vals); md = statistics.median(vals)
                sd = statistics.stdev(vals) if len(vals) > 1 else 0
                thr = thresholds.get(chrom, 0)
                frac_low = sum(1 for v in vals if v <= thr) / len(vals)
                in_v = []; out_v = []
                for c, v in wins:
                    ws = max(0, c - half_w); we = c + half_w
                    ov_frac = compute_roh_overlap(ws, we, roh_ivs)
                    (in_v if ov_frac > 0.5 else out_v).append(v)
                mn_in = statistics.mean(in_v) if in_v else None
                mn_out = statistics.mean(out_v) if out_v else None
                ratio = mn_in/mn_out if mn_in is not None and mn_out and mn_out > 0 else None
                w.writerow([samp, chrom, sid, seg_s, seg_e, len(vals),
                    f"{mn:.8e}", f"{md:.8e}", f"{sd:.8e}", f"{min(vals):.8e}", f"{max(vals):.8e}",
                    roh_bp, f"{froh:.6f}", f"{frac_low:.4f}",
                    f"{mn_in:.8e}" if mn_in is not None else "NA",
                    f"{mn_out:.8e}" if mn_out is not None else "NA",
                    f"{ratio:.6f}" if ratio is not None else "NA"])

def build_chr_breeding_features(theta_data, tracts_idx, callable_per_chr, win_size,
                                 thresholds, scale_label, out_dir):
    """Table H"""
    half_w = win_size // 2
    all_keys = set()
    for samp, rows in theta_data.items():
        for chrom, *_ in rows: all_keys.add((samp, chrom))
    for (samp, chrom) in tracts_idx: all_keys.add((samp, chrom))

    path = os.path.join(out_dir, f"per_sample_chr_breeding_features_{scale_label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample","chrom","callable_bp_chr",
                     "mean_theta_chr","median_theta_chr","frac_low_theta_windows_chr",
                     "roh_bp_chr","froh_chr","n_roh_chr","longest_roh_chr",
                     "mean_theta_in_roh_chr","mean_theta_out_roh_chr","ratio_theta_in_out_chr"])
        for samp, chrom in sorted(all_keys):
            cbp = callable_per_chr.get(chrom, 0)
            roh_ivs = tracts_idx.get((samp, chrom), [])
            roh_bp = sum(max(0, e-s) for s, e in roh_ivs)
            n_roh = len(roh_ivs)
            longest = max((e-s for s, e in roh_ivs), default=0)
            froh = roh_bp / cbp if cbp > 0 else 0
            theta_wins = [(c, tp/ns if ns > 0 else 0)
                          for c_, c, tp, ns in [(chrom, center, tp, ns)
                          for center, tp, ns in [(r[1], r[2], r[3])
                          for r in theta_data.get(samp, []) if r[0] == chrom]]]
            # Simpler rebuild
            theta_wins = []
            for r_chrom, center, tp, ns in theta_data.get(samp, []):
                if r_chrom == chrom:
                    theta_wins.append((center, tp/ns if ns > 0 else 0))
            vals = [v for _, v in theta_wins]
            thr = thresholds.get(chrom, 0)
            if vals:
                mn = statistics.mean(vals); md = statistics.median(vals)
                frac_low = sum(1 for v in vals if v <= thr) / len(vals)
            else:
                mn = md = frac_low = None
            in_v = []; out_v = []
            for c, v in theta_wins:
                ws = max(0, c - half_w); we = c + half_w
                ov = compute_roh_overlap(ws, we, roh_ivs)
                (in_v if ov > 0.5 else out_v).append(v)
            mn_in = statistics.mean(in_v) if in_v else None
            mn_out = statistics.mean(out_v) if out_v else None
            ratio = mn_in/mn_out if mn_in is not None and mn_out and mn_out > 0 else None
            w.writerow([samp, chrom, cbp,
                f"{mn:.8e}" if mn is not None else "NA",
                f"{md:.8e}" if md is not None else "NA",
                f"{frac_low:.4f}" if frac_low is not None else "NA",
                roh_bp, f"{froh:.6f}", n_roh, longest,
                f"{mn_in:.8e}" if mn_in is not None else "NA",
                f"{mn_out:.8e}" if mn_out is not None else "NA",
                f"{ratio:.6f}" if ratio is not None else "NA"])

def build_segment_breeding_features(theta_data, tracts_idx, chrom_sizes, win_size,
                                     thresholds, seg_bp, seg_label, scale_label, out_dir):
    """Table I"""
    half_w = win_size // 2
    segments = []
    for chrom in sorted(chrom_sizes.keys()):
        for i, start in enumerate(range(0, chrom_sizes[chrom], seg_bp)):
            segments.append((chrom, start, min(start+seg_bp, chrom_sizes[chrom]), f"{chrom}_seg{i+1}"))
    samples = sorted(theta_data.keys() | {k[0] for k in tracts_idx})
    path = os.path.join(out_dir, f"per_sample_segment_breeding_features_{seg_label}_{scale_label}.tsv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample","chrom","segment_id","start","end","callable_bp_segment",
                     "mean_theta_segment","median_theta_segment","frac_low_theta_windows_segment",
                     "roh_bp_segment","froh_segment","theta_in_roh_segment","theta_out_roh_segment",
                     "ratio_theta_in_out_segment","low_div_flag","high_roh_flag"])
        for samp in samples:
            by_chr = defaultdict(list)
            for r_chrom, center, tp, ns in theta_data.get(samp, []):
                by_chr[r_chrom].append((center, tp/ns if ns > 0 else 0))
            for chrom, seg_s, seg_e, sid in segments:
                seg_len = seg_e - seg_s
                roh_ivs = tracts_idx.get((samp, chrom), [])
                roh_bp = sum(max(0, min(seg_e,re)-max(seg_s,rs)) for rs, re in roh_ivs if min(seg_e,re)>max(seg_s,rs))
                froh = roh_bp / seg_len if seg_len > 0 else 0
                wins = [(c, v) for c, v in by_chr.get(chrom, []) if seg_s <= c < seg_e]
                vals = [v for _, v in wins]
                thr = thresholds.get(chrom, 0)
                if vals:
                    mn = statistics.mean(vals); md = statistics.median(vals)
                    frac_low = sum(1 for v in vals if v <= thr) / len(vals)
                else: mn = md = frac_low = None
                in_v = []; out_v = []
                for c, v in wins:
                    ws = max(0, c - half_w); we = c + half_w
                    ov = compute_roh_overlap(ws, we, roh_ivs)
                    (in_v if ov > 0.5 else out_v).append(v)
                mn_in = statistics.mean(in_v) if in_v else None
                mn_out = statistics.mean(out_v) if out_v else None
                ratio = mn_in/mn_out if mn_in is not None and mn_out and mn_out > 0 else None
                low_div = 1 if frac_low is not None and frac_low > 0.5 else 0
                high_roh = 1 if froh > 0.5 else 0
                w.writerow([samp, chrom, sid, seg_s, seg_e, seg_len,
                    f"{mn:.8e}" if mn is not None else "NA",
                    f"{md:.8e}" if md is not None else "NA",
                    f"{frac_low:.4f}" if frac_low is not None else "NA",
                    roh_bp, f"{froh:.6f}",
                    f"{mn_in:.8e}" if mn_in is not None else "NA",
                    f"{mn_out:.8e}" if mn_out is not None else "NA",
                    f"{ratio:.6f}" if ratio is not None else "NA",
                    low_div, high_roh])
    return path

def build_pairwise_complementarity(seg_path, seg_label, scale_label, out_dir, max_pairs=0):
    """Tables J, K, L"""
    rows_by_sample = defaultdict(list)
    with open(seg_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader: rows_by_sample[row["sample"]].append(row)
    samples = sorted(rows_by_sample.keys())
    pairs = list(itertools.combinations(samples, 2))
    if max_pairs > 0 and len(pairs) > max_pairs:
        print(f"  Limiting to {max_pairs} pairs (of {len(pairs)} total)")
        pairs = pairs[:max_pairs]
    chr_agg = defaultdict(lambda: {"n_seg":0,"n_both_bad":0,"n_comp":0,"scores":[]})
    genome_agg = defaultdict(lambda: {"n_chr":set(),"n_both_bad":0,"n_comp":0,"scores":[]})
    seg_out = os.path.join(out_dir, f"pairwise_segment_complementarity_{seg_label}_{scale_label}.tsv")
    with open(seg_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_A","sample_B","chrom","segment_id","start","end",
                     "A_mean_theta","B_mean_theta","A_roh_bp","B_roh_bp",
                     "A_low_div","B_low_div","A_high_roh","B_high_roh",
                     "both_bad_flag","complementarity_flag","complementarity_score"])
        for sa, sb in pairs:
            ra = {(r["chrom"],r["segment_id"]): r for r in rows_by_sample[sa]}
            rb = {(r["chrom"],r["segment_id"]): r for r in rows_by_sample[sb]}
            for key in sorted(set(ra)|set(rb)):
                a = ra.get(key, {}); b = rb.get(key, {})
                chrom = a.get("chrom", b.get("chrom",""))
                sid = a.get("segment_id", b.get("segment_id",""))
                start = a.get("start", b.get("start","")); end = a.get("end", b.get("end",""))
                at = a.get("mean_theta_segment","NA"); bt = b.get("mean_theta_segment","NA")
                ar = int(a.get("roh_bp_segment",0)); br = int(b.get("roh_bp_segment",0))
                al = int(a.get("low_div_flag",0)); bl = int(b.get("low_div_flag",0))
                ah = int(a.get("high_roh_flag",0)); bh = int(b.get("high_roh_flag",0))
                a_bad = al or ah; b_bad = bl or bh
                both = 1 if a_bad and b_bad else 0
                comp = 1 if (a_bad and not b_bad) or (b_bad and not a_bad) else 0
                try: score = abs(float(at)-float(bt))
                except: score = 0
                w.writerow([sa,sb,chrom,sid,start,end,at,bt,ar,br,al,bl,ah,bh,both,comp,f"{score:.8e}"])
                chr_agg[(sa,sb,chrom)]["n_seg"] += 1
                chr_agg[(sa,sb,chrom)]["n_both_bad"] += both
                chr_agg[(sa,sb,chrom)]["n_comp"] += comp
                chr_agg[(sa,sb,chrom)]["scores"].append(score)
                genome_agg[(sa,sb)]["n_chr"].add(chrom)
                genome_agg[(sa,sb)]["n_both_bad"] += both
                genome_agg[(sa,sb)]["n_comp"] += comp
                genome_agg[(sa,sb)]["scores"].append(score)
    # Chr level
    chr_out = os.path.join(out_dir, f"pairwise_chr_complementarity_{seg_label}_{scale_label}.tsv")
    with open(chr_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_A","sample_B","chrom","n_segments","n_both_bad","n_complementary",
                     "mean_score","risk_score","rescue_score"])
        for (sa,sb,ch), d in sorted(chr_agg.items()):
            ms = statistics.mean(d["scores"]) if d["scores"] else 0
            risk = d["n_both_bad"]/d["n_seg"] if d["n_seg"] else 0
            rescue = d["n_comp"]/d["n_seg"] if d["n_seg"] else 0
            w.writerow([sa,sb,ch,d["n_seg"],d["n_both_bad"],d["n_comp"],f"{ms:.8e}",f"{risk:.4f}",f"{rescue:.4f}"])
    # Genome level
    gen_out = os.path.join(out_dir, f"pairwise_genome_complementarity_{seg_label}_{scale_label}.tsv")
    with open(gen_out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_A","sample_B","n_chr","total_both_bad","total_complementary",
                     "mean_score","risk_score","rescue_score"])
        for (sa,sb), d in sorted(genome_agg.items()):
            n = len(d["scores"])
            ms = statistics.mean(d["scores"]) if d["scores"] else 0
            risk = d["n_both_bad"]/n if n else 0; rescue = d["n_comp"]/n if n else 0
            w.writerow([sa,sb,len(d["n_chr"]),d["n_both_bad"],d["n_comp"],f"{ms:.8e}",f"{risk:.4f}",f"{rescue:.4f}"])

# ═══════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracts-bed", required=True)
    ap.add_argument("--callable-bed", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--theta-dir", required=True)
    ap.add_argument("--theta-main-dir", default="")
    ap.add_argument("--sample-list", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--ancestry-labels", default="")
    ap.add_argument("--order-file", default="")
    ap.add_argument("--pruned81", default="")
    ap.add_argument("--max-pairs", type=int, default=0)
    ap.add_argument("--theta-scales", nargs="+", default=["5000_1000","10000_2000","50000_10000"])
    ap.add_argument("--theta-labels", nargs="+", default=["5kb_1kb","10kb_2kb","50kb_10kb"])
    ap.add_argument("--segment-sizes", nargs="+", type=int, default=[500000,1000000,2000000,5000000])
    ap.add_argument("--segment-labels", nargs="+", default=["500kb","1Mb","2Mb","5Mb"])
    ap.add_argument("--recurrence-sizes", nargs="+", type=int, default=[100000,250000,500000,1000000])
    ap.add_argument("--recurrence-labels", nargs="+", default=["100kb","250kb","500kb","1Mb"])
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("=== Loading data ===")
    tracts = load_tracts(args.tracts_bed)
    chrom_sizes = load_chrom_sizes(args.chrom_sizes)
    callable_per_chr = load_callable_per_chr(args.callable_bed)
    callable_total = sum(callable_per_chr.values())
    samples = load_list(args.sample_list)
    order = load_map(args.order_file)
    ancestry = load_ancestry(args.ancestry_labels)
    pruned = set(load_list(args.pruned81))
    tracts_idx = index_tracts_by_sample_chr(tracts)
    print(f"  {len(tracts)} tracts, {len(samples)} samples, {len(chrom_sizes)} chroms, callable={callable_total:,} bp")

    # A. ROH bins
    print("\n=== A. ROH bins (fixed + adaptive) ===")
    all_lengths = [l for *_, l in tracts if l >= 1_000_000]
    adaptive_bins = derive_adaptive_bins(all_lengths, args.out_dir)
    print(f"  Adaptive: {[b[0] for b in adaptive_bins]}")
    build_roh_bins_wide(tracts, callable_total, FIXED_BINS, "fixedBins", args.out_dir)
    build_roh_bins_wide(tracts, callable_total, adaptive_bins, "adaptiveBins", args.out_dir)

    # B. Genome ROH segments
    print("\n=== B. Genome ROH segments ===")
    build_genome_roh_segments(tracts, order, args.out_dir)

    # C. Per-chr ROH group summary
    print("\n=== C. Per-chr ROH group summary ===")
    build_chr_roh_group_summary(tracts, callable_per_chr, samples, ancestry, pruned, args.out_dir)

    # D. Recurrence windows
    print("\n=== D. Recurrence windows ===")
    for win_bp, label in zip(args.recurrence_sizes, args.recurrence_labels):
        print(f"  {label}...")
        build_recurrence_windows(tracts, chrom_sizes, samples, win_bp, label, args.out_dir, ancestry, pruned)

    # Theta-dependent tables
    for scale, scale_label in zip(args.theta_scales, args.theta_labels):
        win_str, step_str = scale.split("_")
        win_size = int(win_str)
        print(f"\n=== Theta scale: {scale_label} ===")
        theta_files = find_theta_files(args.theta_dir, samples, win_str, step_str)
        if not theta_files and args.theta_main_dir:
            theta_files = find_theta_files(args.theta_main_dir, samples, win_str, step_str)
        if not theta_files:
            print(f"  No theta files for {scale_label}, skipping"); continue
        print(f"  Found {len(theta_files)} theta files")
        theta_data = {s: parse_theta_file(fp) for s, fp in theta_files.items()}

        # Shared thresholds
        thresholds = compute_chr_theta_thresholds(theta_data, win_size, scale_label, args.out_dir)

        # E. Per-sample window theta
        print(f"  E. per_sample_window_theta_{scale_label}")
        build_per_sample_window_theta(theta_data, tracts_idx, win_size, scale_label,
                                       args.out_dir, order, ancestry, pruned)

        # F. Per-chr theta summary
        print(f"  F. per_chr_theta_summary_{scale_label}")
        build_per_chr_theta_summary(theta_data, tracts_idx, win_size, thresholds, scale_label, args.out_dir)

        # H. Chr breeding features
        print(f"  H. per_sample_chr_breeding_features_{scale_label}")
        build_chr_breeding_features(theta_data, tracts_idx, callable_per_chr, win_size,
                                     thresholds, scale_label, args.out_dir)

        for seg_bp, seg_label in zip(args.segment_sizes, args.segment_labels):
            # G. Segment theta summary
            print(f"  G. segment_theta_summary {seg_label}_{scale_label}")
            build_chr_segment_theta_summary(theta_data, tracts_idx, chrom_sizes, win_size,
                                             thresholds, seg_bp, seg_label, scale_label, args.out_dir)

            # I. Segment breeding features
            print(f"  I. segment_breeding_features {seg_label}_{scale_label}")
            seg_path = build_segment_breeding_features(theta_data, tracts_idx, chrom_sizes, win_size,
                                                        thresholds, seg_bp, seg_label, scale_label, args.out_dir)

            # J,K,L. Pairwise
            if args.max_pairs != -1:
                print(f"  J,K,L. pairwise {seg_label}_{scale_label}")
                build_pairwise_complementarity(seg_path, seg_label, scale_label, args.out_dir, args.max_pairs)

    print("\n=== Breeding tables complete ===")

if __name__ == "__main__":
    main()
