#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
summarize_ibd_roh.py

Summarize ngsF-HMM IBD tracts (ROH-like) into per-sample:
- total ROH bp, FROH (normalized by callable ∩ nonrepeat span)
- length-bin contributions like paper: 1–2, 2–4, 4–8, 8–16, >16 Mb
- optional ancestry label from Q-matrix
- optional sample ordering from an order file (e.g., cov-tree leaf order)
- optional relatedness annotation from kinship/relatedness table

Inputs:
  1) IBD tracts BED from convert_ibd.pl:
     chr   start   end   sample   length(optional)
  2) callable BED, repeats BED (or directly an "accessible bed" = callable minus repeats)
  3) (optional) Q matrix: sample + K columns, or K columns with separate sample list
  4) (optional) sample order file (one sample per line)
  5) (optional) kinship table (pairwise)

Outputs:
  - per_sample_roh.tsv
  - per_sample_roh_bins_long.tsv (plot-ready long format for stacked bars)

Notes:
  - This script does NOT call ROH from genotypes; it only summarizes tract intervals.
  - For masking: it will compute "accessible_nonrepeat" = callable minus repeats
    and (optionally) intersect tracts with it via bedtools if you choose.
"""

import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional


# -----------------------------
# Helpers
# -----------------------------

def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def run(cmd: List[str]) -> None:
    """Run a command, raise on failure."""
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        die(f"Command not found: {cmd[0]}")
    except subprocess.CalledProcessError as e:
        die(f"Command failed ({e.returncode}): {' '.join(cmd)}")


def bed_sum_length(bed_path: str) -> int:
    """Sum interval lengths (end-start) in a BED-like file (chr start end ...)."""
    total = 0
    with open(bed_path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 3:
                continue
            s = int(parts[1])
            e = int(parts[2])
            if e > s:
                total += (e - s)
    return total


def ensure_sorted_bed(in_bed: str, out_bed: str) -> None:
    """Sort a BED file using bedtools sort."""
    run(["bedtools", "sort", "-i", in_bed])
    # bedtools sort prints to stdout; re-run piping is annoying without shell.
    # Use Unix sort fallback if bedtools sort can't write:
    # We'll just do a safe GNU sort here.
    with open(out_bed, "w") as out:
        run_with_stdout(["sort", "-k1,1", "-k2,2n", in_bed], out)


def run_with_stdout(cmd: List[str], out_handle) -> None:
    try:
        subprocess.run(cmd, check=True, stdout=out_handle)
    except FileNotFoundError:
        die(f"Command not found: {cmd[0]}")
    except subprocess.CalledProcessError as e:
        die(f"Command failed ({e.returncode}): {' '.join(cmd)}")


def build_accessible_nonrepeat(callable_bed: str, repeats_bed: str, out_bed: str) -> None:
    """
    accessible = callable - repeats, then sort+merge
    """
    tmp_sub = out_bed + ".tmp.sub.bed"
    tmp_sort = out_bed + ".tmp.sort.bed"

    # bedtools subtract
    with open(tmp_sub, "w") as out:
        run_with_stdout(["bedtools", "subtract", "-a", callable_bed, "-b", repeats_bed], out)

    # sort
    with open(tmp_sort, "w") as out:
        run_with_stdout(["sort", "-k1,1", "-k2,2n", tmp_sub], out)

    # merge
    with open(out_bed, "w") as out:
        run_with_stdout(["bedtools", "merge", "-i", tmp_sort], out)

    for p in (tmp_sub, tmp_sort):
        try:
            os.remove(p)
        except OSError:
            pass


def intersect_tracts_with_accessible(tracts_bed: str, accessible_bed: str, out_bed: str) -> None:
    """
    Intersect tracts with accessible regions. Keeps tract sample IDs.
    Input tracts bed: chr start end sample [length]
    Output: chr start end sample length
    """
    tmp = out_bed + ".tmp.intersect.bed"

    # -wa keeps A features, but the resulting intervals are A intervals, not clipped.
    # We want clipped intersections, so use -wa -wb and then compute intersection start/end.
    # bedtools intersect with -wao yields extra, but easiest is:
    # bedtools intersect -a tracts -b accessible
    # This returns clipped A intervals (i.e., the overlap region) by default, good.
    # However, it prints original A columns, but start/end are overlap.
    # Great: "bedtools intersect" outputs the intersection intervals.

    with open(tmp, "w") as out:
        run_with_stdout(["bedtools", "intersect", "-a", tracts_bed, "-b", accessible_bed], out)

    # Now normalize to 5 columns with length
    with open(out_bed, "w") as out, open(tmp, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 4:
                continue
            chr_, s, e, sample = parts[0], int(parts[1]), int(parts[2]), parts[3]
            length = e - s
            if length <= 0:
                continue
            out.write(f"{chr_}\t{s}\t{e}\t{sample}\t{length}\n")

    try:
        os.remove(tmp)
    except OSError:
        pass


@dataclass
class QInfo:
    K: int
    qvals: Dict[str, List[float]]          # sample -> [q1..qK]
    label: Dict[str, str]                   # sample -> label string
    label_k: Dict[str, int]                 # sample -> argmax cluster index (1..K)


def read_q_matrix(q_path: str,
                  q_has_sample_col: bool,
                  sample_list: Optional[str],
                  label_prefix: str = "anc") -> QInfo:
    """
    Reads a Q matrix.
    Supported formats:
      A) sample q1 q2 ... qK        (q_has_sample_col=True)
      B) q1 q2 ... qK              (q_has_sample_col=False) with sample_list provided

    Returns per-sample q-values and ancestry label = argmax cluster.
    """
    samples: List[str] = []
    if not q_has_sample_col:
        if not sample_list:
            die("If --q-has-sample-col is false, you must provide --q-samples")
        with open(sample_list, "r") as f:
            samples = [ln.strip().split()[0] for ln in f if ln.strip()]

    qvals: Dict[str, List[float]] = {}
    with open(q_path, "r") as f:
        rows = [ln.strip().split() for ln in f if ln.strip()]

    if not rows:
        die(f"Empty Q file: {q_path}")

    if q_has_sample_col:
        for r in rows:
            if len(r) < 2:
                continue
            s = r[0]
            qs = [float(x) for x in r[1:]]
            qvals[s] = qs
    else:
        if len(rows) != len(samples):
            die(f"Q rows ({len(rows)}) != samples ({len(samples)}). Check --q-samples.")
        for s, r in zip(samples, rows):
            qs = [float(x) for x in r]
            qvals[s] = qs

    # Determine K
    K = len(next(iter(qvals.values())))
    label: Dict[str, str] = {}
    label_k: Dict[str, int] = {}
    for s, qs in qvals.items():
        if len(qs) != K:
            die(f"Inconsistent K for sample {s} in Q matrix.")
        k = max(range(K), key=lambda i: qs[i])  # 0-based
        label_k[s] = k + 1
        label[s] = f"{label_prefix}{k+1}"
    return QInfo(K=K, qvals=qvals, label=label, label_k=label_k)


def read_order_file(order_path: str) -> Dict[str, int]:
    """
    One sample per line; returns sample -> order index (0..n-1).
    """
    order: Dict[str, int] = {}
    with open(order_path, "r") as f:
        for i, ln in enumerate([x.strip() for x in f if x.strip()]):
            s = ln.split()[0]
            order[s] = i
    return order


@dataclass
class KinshipInfo:
    max_kin: Dict[str, float]
    n_close: Dict[str, int]


def read_kinship(kin_path: str, threshold: float) -> KinshipInfo:
    """
    Reads a pairwise kinship table with at least 3 cols: sample1 sample2 kinship
    Returns per-sample:
      - max kinship to any other sample
      - count of pairs above threshold
    """
    max_kin: Dict[str, float] = defaultdict(lambda: 0.0)
    n_close: Dict[str, int] = defaultdict(int)

    with open(kin_path, "r") as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.strip().split()
            if len(parts) < 3:
                continue
            a, b = parts[0], parts[1]
            try:
                k = float(parts[2])
            except ValueError:
                continue
            if a != b:
                if k > max_kin[a]:
                    max_kin[a] = k
                if k > max_kin[b]:
                    max_kin[b] = k
                if k >= threshold:
                    n_close[a] += 1
                    n_close[b] += 1

    return KinshipInfo(max_kin=dict(max_kin), n_close=dict(n_close))


# -----------------------------
# ROH binning
# -----------------------------

BINS = [
    ("1-2Mb", 1_000_000, 2_000_000),
    ("2-4Mb", 2_000_000, 4_000_000),
    ("4-8Mb", 4_000_000, 8_000_000),
    ("8-16Mb", 8_000_000, 16_000_000),
    (">16Mb", 16_000_000, None),
]


def bin_length_bp(length: int) -> Optional[str]:
    for name, lo, hi in BINS:
        if length >= lo and (hi is None or length < hi):
            return name
    return None  # <1Mb


def summarize_tracts(tracts_bed: str, min_roh_bp: int) -> Tuple[Dict[str, int], Dict[str, Dict[str, int]]]:
    """
    Returns:
      total_bp[sample] = total ROH bp (>=min_roh_bp applied at tract level)
      binned_bp[sample][binname] = summed bp in each bin
    Input bed must have at least: chr start end sample length
    """
    total_bp: Dict[str, int] = defaultdict(int)
    binned_bp: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    with open(tracts_bed, "r") as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.strip().split()
            if len(parts) < 5:
                # allow 4 columns and compute length
                if len(parts) < 4:
                    continue
                s = int(parts[1]); e = int(parts[2])
                length = e - s
                sample = parts[3]
            else:
                sample = parts[3]
                length = int(float(parts[4]))
            if length <= 0:
                continue
            if length < min_roh_bp:
                continue

            total_bp[sample] += length
            b = bin_length_bp(length)
            if b is not None:
                binned_bp[sample][b] += length

    # Ensure bins exist as 0
    for sample in list(total_bp.keys()):
        for b, _, _ in BINS:
            _ = binned_bp[sample][b]  # touch to create
    return dict(total_bp), {s: dict(d) for s, d in binned_bp.items()}


# -----------------------------
# Main
# -----------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracts-bed", required=True,
                    help="BED tracts from convert_ibd.pl: chr start end sample [len]")
    ap.add_argument("--callable-bed", required=True,
                    help="Callable regions BED")
    ap.add_argument("--repeats-bed", required=True,
                    help="Repeats BED to exclude")
    ap.add_argument("--out-prefix", default="roh",
                    help="Output prefix (default: roh)")

    ap.add_argument("--min-roh-bp", type=int, default=1_000_000,
                    help="Minimum tract length to keep (default 1,000,000)")
    ap.add_argument("--use-bedtools-intersect", action="store_true",
                    help="Intersect tracts with accessible_nonrepeat using bedtools (recommended)")
    ap.add_argument("--merge-gap", type=int, default=0,
                    help="Optional: merge tracts within this gap after intersect (bp). 0=off.")

    # Q matrix / ancestry
    ap.add_argument("--q-matrix", default=None,
                    help="Q matrix file (optional)")
    ap.add_argument("--q-has-sample-col", action="store_true",
                    help="Q file has sample name as first column")
    ap.add_argument("--q-samples", default=None,
                    help="Sample list (one per line) if Q matrix has no sample column")
    ap.add_argument("--label-prefix", default="anc",
                    help="Prefix for ancestry label (default 'anc')")

    # Ordering
    ap.add_argument("--order-file", default=None,
                    help="Sample order file (one sample per line), e.g. cov-tree leaf order")

    # Relatedness
    ap.add_argument("--kinship", default=None,
                    help="Pairwise kinship table: sample1 sample2 kinship")
    ap.add_argument("--kin-threshold", type=float, default=0.0884,
                    help="Kinship threshold to count 'close relatives' (default 0.0884 ~ 2nd degree-ish)")

    args = ap.parse_args()

    # Prepare accessible bed
    out_accessible = f"{args.out_prefix}.accessible_nonrepeat.bed"
    build_accessible_nonrepeat(args.callable_bed, args.repeats_bed, out_accessible)
    callable_bp = bed_sum_length(out_accessible)
    if callable_bp <= 0:
        die("Callable_nonrepeat span is 0. Check your BED inputs.")

    # Optionally intersect tracts with accessible bed
    tracts_for_summary = args.tracts_bed
    if args.use_bedtools_intersect:
        out_inter = f"{args.out_prefix}.tracts.accessible.bed"
        intersect_tracts_with_accessible(args.tracts_bed, out_accessible, out_inter)
        tracts_for_summary = out_inter

        # Optional merge within gap
        if args.merge_gap and args.merge_gap > 0:
            tmp_sort = f"{args.out_prefix}.tracts.accessible.sorted.bed"
            tmp_merge = f"{args.out_prefix}.tracts.accessible.merged.bed"

            # sort by chr,start then sample; merge per sample requires grouping.
            # simplest: split by sample later. For now, do per-sample merge in python-style.
            # We'll do a safe per-sample merge without bedtools grouping.
            tracts_for_summary = merge_tracts_per_sample(tmp_in=out_inter,
                                                         out_path=tmp_merge,
                                                         max_gap=args.merge_gap)

    # Read Q matrix if provided
    qinfo: Optional[QInfo] = None
    if args.q_matrix:
        qinfo = read_q_matrix(args.q_matrix, args.q_has_sample_col, args.q_samples, args.label_prefix)

    # Read order file if provided
    order: Dict[str, int] = {}
    if args.order_file:
        order = read_order_file(args.order_file)

    # Read kinship if provided
    kin: Optional[KinshipInfo] = None
    if args.kinship:
        kin = read_kinship(args.kinship, args.kin_threshold)

    # Summarize tracts
    total_bp, binned_bp = summarize_tracts(tracts_for_summary, args.min_roh_bp)

    # Collect all sample IDs seen in tracts OR metadata
    all_samples = set(total_bp.keys()) | set(binned_bp.keys())
    if qinfo:
        all_samples |= set(qinfo.qvals.keys())
    if order:
        all_samples |= set(order.keys())
    if kin:
        all_samples |= set(kin.max_kin.keys()) | set(kin.n_close.keys())

    # Prepare output rows
    out_per = f"{args.out_prefix}.per_sample_roh.tsv"
    out_long = f"{args.out_prefix}.per_sample_roh_bins_long.tsv"

    # Define header
    q_cols: List[str] = []
    if qinfo:
        q_cols = [f"Q{k+1}" for k in range(qinfo.K)]

    header = [
        "sample",
        "order_index",
        "ancestry_label",
        *q_cols,
        "max_kinship",
        "n_close_rel",
        "n_roh_bins_ge1Mb",  # count of tracts >= min_roh_bp
        "roh_total_bp",
        "callable_nonrepeat_bp",
        "FROH",
    ] + [f"roh_bp_{b[0]}" for b in BINS] + [f"pct_{b[0]}" for b in BINS]

    # Count tracts per sample (>=min) by scanning bed once
    tract_count = count_tracts_per_sample(tracts_for_summary, args.min_roh_bp)

    # Sort samples: order_index if provided, else alphabetical
    def sort_key(s: str):
        return (order.get(s, 10**12), s)

    samples_sorted = sorted(all_samples, key=sort_key)

    with open(out_per, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)

        for s in samples_sorted:
            roh_bp = total_bp.get(s, 0)
            froh = roh_bp / callable_bp if callable_bp > 0 else "NA"

            # ancestry
            anc_label = qinfo.label.get(s, "NA") if qinfo else "NA"
            qvals = qinfo.qvals.get(s, []) if qinfo else []
            if qinfo and not qvals:
                qvals = [float("nan")] * qinfo.K

            # relatedness
            maxk = kin.max_kin.get(s, 0.0) if kin else 0.0
            nclose = kin.n_close.get(s, 0) if kin else 0

            # bins
            perbin = binned_bp.get(s, {})
            roh_bp_bins = []
            pct_bins = []
            for bname, _, _ in BINS:
                x = int(perbin.get(bname, 0))
                roh_bp_bins.append(x)
                pct_bins.append(100.0 * x / callable_bp if callable_bp > 0 else float("nan"))

            row = [
                s,
                order.get(s, -1),
                anc_label,
                *qvals,
                maxk,
                nclose,
                tract_count.get(s, 0),
                roh_bp,
                callable_bp,
                froh,
                *roh_bp_bins,
                *pct_bins
            ]
            w.writerow(row)

    # Long format for plotting stacked bars
    with open(out_long, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "order_index", "ancestry_label", "bin", "roh_bp", "pct_genome_in_bin"])

        for s in samples_sorted:
            anc_label = qinfo.label.get(s, "NA") if qinfo else "NA"
            for bname, _, _ in BINS:
                x = int(binned_bp.get(s, {}).get(bname, 0))
                pct = 100.0 * x / callable_bp if callable_bp > 0 else float("nan")
                w.writerow([s, order.get(s, -1), anc_label, bname, x, pct])

    print("Wrote:")
    print(f"  {out_accessible}")
    if args.use_bedtools_intersect:
        print(f"  {args.out_prefix}.tracts.accessible.bed")
    print(f"  {out_per}")
    print(f"  {out_long}")


def count_tracts_per_sample(tracts_bed: str, min_len: int) -> Dict[str, int]:
    c = defaultdict(int)
    with open(tracts_bed, "r") as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.strip().split()
            if len(parts) < 4:
                continue
            sample = parts[3]
            if len(parts) >= 5:
                try:
                    length = int(float(parts[4]))
                except ValueError:
                    length = int(parts[2]) - int(parts[1])
            else:
                length = int(parts[2]) - int(parts[1])
            if length >= min_len:
                c[sample] += 1
    return dict(c)


def merge_tracts_per_sample(tmp_in: str, out_path: str, max_gap: int) -> str:
    """
    Merge tracts per sample if separated by <= max_gap, per chromosome.
    Keeps 5-col output with recomputed length.
    """
    # load per sample
    per = defaultdict(list)  # sample -> list of (chr,s,e)
    with open(tmp_in, "r") as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"):
                continue
            parts = ln.strip().split()
            if len(parts) < 4:
                continue
            chr_ = parts[0]
            s = int(parts[1]); e = int(parts[2])
            sample = parts[3]
            if e <= s:
                continue
            per[sample].append((chr_, s, e))

    with open(out_path, "w") as out:
        for sample, segs in per.items():
            # sort by chr then start
            segs.sort(key=lambda x: (x[0], x[1], x[2]))
            curr_chr = None
            curr_s = None
            curr_e = None
            for chr_, s, e in segs:
                if curr_chr is None:
                    curr_chr, curr_s, curr_e = chr_, s, e
                    continue
                if chr_ != curr_chr or s > curr_e + max_gap:
                    # flush
                    out.write(f"{curr_chr}\t{curr_s}\t{curr_e}\t{sample}\t{curr_e-curr_s}\n")
                    curr_chr, curr_s, curr_e = chr_, s, e
                else:
                    # merge/extend
                    if e > curr_e:
                        curr_e = e
            if curr_chr is not None:
                out.write(f"{curr_chr}\t{curr_s}\t{curr_e}\t{sample}\t{curr_e-curr_s}\n")

    return out_path


if __name__ == "__main__":
    main()
