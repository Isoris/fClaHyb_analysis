###
# pairwise_first_degree_curation.py does this:
# Reads pairwise relatedness data (sample1, sample2, theta), flags pairs above
# a chosen theta threshold, and greedily decides which samples to remove using
# keep/remove preferences, QC scores, graph degree, and deterministic tie-breaks.
# It writes keep/remove lists plus flagged-pair, decision, and summary tables.
###
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
from collections import Counter

def load_set(path):
    s = set()
    if not path:
        return s
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                s.add(line.split()[0])
    return s

def load_qc(path):
    """
    TSV/space-delimited file:
    sample<TAB>score

    Higher score = better sample to KEEP
    """
    qc = {}
    if not path:
        return qc
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            sample = parts[0]
            try:
                score = float(parts[1])
            except ValueError:
                continue
            qc[sample] = score
    return qc

def read_pairs(path):
    """
    Input format:
    sample1  sample2  theta
    """
    pairs = []
    samples = set()
    with open(path) as f:
        for i, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"Line {i}: expected at least 3 columns, got: {line}")
            a, b = parts[0], parts[1]
            try:
                theta = float(parts[2])
            except ValueError:
                raise ValueError(f"Line {i}: third column is not numeric: {parts[2]}")
            pairs.append((a, b, theta))
            samples.add(a)
            samples.add(b)
    return pairs, samples

def choose_remove(a, b, degree, qc, prefer_keep, prefer_remove):
    """
    Decide which sample to remove from a pair.
    Returns: (remove_sample, keep_sample, reason)
    Priority:
      1) prefer_keep list
      2) prefer_remove list
      3) QC score (higher score kept)
      4) degree among flagged pairs (higher degree removed)
      5) lexical tie-break
    """

    # 1) prefer_keep list
    a_keep = a in prefer_keep
    b_keep = b in prefer_keep
    if a_keep and not b_keep:
        return b, a, "prefer_keep"
    if b_keep and not a_keep:
        return a, b, "prefer_keep"

    # 2) prefer_remove list
    a_drop = a in prefer_remove
    b_drop = b in prefer_remove
    if a_drop and not b_drop:
        return a, b, "prefer_remove"
    if b_drop and not a_drop:
        return b, a, "prefer_remove"

    # 3) QC score (higher = better, so remove lower)
    a_qc = qc.get(a, None)
    b_qc = qc.get(b, None)
    if a_qc is not None and b_qc is not None and a_qc != b_qc:
        if a_qc < b_qc:
            return a, b, "lower_qc"
        else:
            return b, a, "lower_qc"
    elif a_qc is not None and b_qc is None:
        return b, a, "missing_qc"
    elif b_qc is not None and a_qc is None:
        return a, b, "missing_qc"

    # 4) remove hub if tied (more connected among flagged pairs)
    if degree[a] != degree[b]:
        if degree[a] > degree[b]:
            return a, b, "higher_degree_hub"
        else:
            return b, a, "higher_degree_hub"

    # 5) deterministic lexical fallback
    # keep lexicographically smaller, remove larger
    if a <= b:
        return b, a, "lexical_tie"
    else:
        return a, b, "lexical_tie"

def main():
    ap = argparse.ArgumentParser(
        description="Pairwise first-degree curation from sample1 sample2 theta input."
    )
    ap.add_argument("-i", "--input", required=True,
                    help="Input file: sample1 sample2 theta")
    ap.add_argument("-o", "--output-prefix", required=True,
                    help="Prefix for output files")
    ap.add_argument("-t", "--threshold", type=float, default=0.177,
                    help="Theta threshold to flag pairs (default: 0.177 for first-degree)")
    ap.add_argument("--qc", default=None,
                    help="Optional QC file: sample score (higher score kept)")
    ap.add_argument("--prefer-keep", default=None,
                    help="Optional file with sample IDs that should be protected if possible")
    ap.add_argument("--prefer-remove", default=None,
                    help="Optional file with sample IDs to sacrifice first if possible")
    args = ap.parse_args()

    pairs, all_samples = read_pairs(args.input)
    qc = load_qc(args.qc)
    prefer_keep = load_set(args.prefer_keep)
    prefer_remove = load_set(args.prefer_remove)

    # Flag first-degree-or-closer pairs
    flagged = [(a, b, theta) for (a, b, theta) in pairs if theta >= args.threshold]

    # Sort highest theta first
    flagged.sort(key=lambda x: (-x[2], x[0], x[1]))

    # Degree among flagged pairs only
    degree = Counter()
    for a, b, theta in flagged:
        degree[a] += 1
        degree[b] += 1

    kept = set(all_samples)
    removed = set()
    decisions = []
    skipped_pairs = 0

    # Greedy pairwise curation
    for a, b, theta in flagged:
        # If the pair is already resolved, skip
        if a not in kept or b not in kept:
            skipped_pairs += 1
            continue

        rm, kp, reason = choose_remove(a, b, degree, qc, prefer_keep, prefer_remove)
        kept.remove(rm)
        removed.add(rm)

        decisions.append({
            "sample1": a,
            "sample2": b,
            "theta": theta,
            "removed": rm,
            "kept": kp,
            "reason": reason,
            "degree_removed": degree[rm],
            "degree_kept": degree[kp],
            "qc_removed": qc.get(rm, ""),
            "qc_kept": qc.get(kp, "")
        })

    # Outputs
    remove_path = args.output_prefix + "_toRemove.txt"
    keep_path = args.output_prefix + "_toKeep.txt"
    decisions_path = args.output_prefix + "_decisions.tsv"
    summary_path = args.output_prefix + "_summary.tsv"
    flagged_path = args.output_prefix + "_flagged_pairs.tsv"

    with open(remove_path, "w") as f:
        for s in sorted(removed):
            f.write(s + "\n")

    with open(keep_path, "w") as f:
        for s in sorted(kept):
            f.write(s + "\n")

    with open(flagged_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample1", "sample2", "theta"])
        for a, b, theta in flagged:
            w.writerow([a, b, f"{theta:.6f}"])

    with open(decisions_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "sample1", "sample2", "theta",
            "removed", "kept", "reason",
            "degree_removed", "degree_kept",
            "qc_removed", "qc_kept"
        ])
        for d in decisions:
            w.writerow([
                d["sample1"], d["sample2"], f'{d["theta"]:.6f}',
                d["removed"], d["kept"], d["reason"],
                d["degree_removed"], d["degree_kept"],
                d["qc_removed"], d["qc_kept"]
            ])

    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["threshold", args.threshold])
        w.writerow(["total_samples", len(all_samples)])
        w.writerow(["total_pairs", len(pairs)])
        w.writerow(["flagged_pairs_theta_ge_threshold", len(flagged)])
        w.writerow(["removed_samples", len(removed)])
        w.writerow(["retained_samples", len(kept)])
        w.writerow(["skipped_pairs_already_resolved", skipped_pairs])

    print("[DONE]")
    print(f"Threshold: {args.threshold}")
    print(f"Total samples: {len(all_samples)}")
    print(f"Total pairs: {len(pairs)}")
    print(f"Flagged pairs (theta >= threshold): {len(flagged)}")
    print(f"Removed samples: {len(removed)}")
    print(f"Retained samples: {len(kept)}")
    print(f"Wrote: {remove_path}")
    print(f"Wrote: {keep_path}")
    print(f"Wrote: {flagged_path}")
    print(f"Wrote: {decisions_path}")
    print(f"Wrote: {summary_path}")

if __name__ == "__main__":
    main()
