#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_parse_roh_and_het.py — Convert .ibd → BED ROH, compute FROH, het in/out ROH

Wraps:
  1) convert_ibd.pl  → per-sample ROH BED tracts
  2) summarize_ibd_roh.py → ROH/FROH summary with bins
  3) Intersect local theta with ROH intervals → het inside/outside ROH
  4) Merge into master summary table

Usage:
  python3 04_parse_roh_and_het.py --config 00_config.sh

Or with explicit paths:
  python3 04_parse_roh_and_het.py \\
    --ibd best.ibd \\
    --pos main_qcpass.pos \\
    --ind samples.ind \\
    --callable-bed callable.bed \\
    --repeats-bed repeats.bed \\
    --het-summary genomewide_heterozygosity.tsv \\
    --theta-dir 02_heterozygosity/03_theta \\
    --out-dir 04_roh_summary \\
    --tables-dir 09_final_tables \\
    --convert-ibd-pl /path/to/convert_ibd.pl \\
    --summarize-roh-py /path/to/summarize_ibd_roh.py
"""

import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def die(msg, code=1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def run_cmd(cmd, desc=""):
    """Run a shell command, print and check."""
    print(f"  Running: {' '.join(cmd)}" + (f"  ({desc})" if desc else ""))
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        die(f"Command not found: {cmd[0]}")
    except subprocess.CalledProcessError as e:
        die(f"Command failed ({e.returncode}): {' '.join(cmd)}")


def run_cmd_stdout(cmd, outpath, desc=""):
    """Run command and capture stdout to file."""
    print(f"  Running: {' '.join(cmd)} > {outpath}" + (f"  ({desc})" if desc else ""))
    with open(outpath, "w") as out:
        try:
            subprocess.run(cmd, check=True, stdout=out)
        except FileNotFoundError:
            die(f"Command not found: {cmd[0]}")
        except subprocess.CalledProcessError as e:
            die(f"Command failed ({e.returncode}): {' '.join(cmd)}")


def parse_het_summary(het_tsv):
    """Read genomewide_heterozygosity.tsv → dict of sample→het."""
    het = {}
    with open(het_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            s = row.get("sample", "")
            h = row.get("het_genomewide", row.get("het_saf1", "NA"))
            try:
                het[s] = float(h)
            except (ValueError, TypeError):
                het[s] = None
    return het


def compute_per_chr_roh(tracts_bed, callable_bed):
    """
    Compute per-sample per-chromosome ROH summary.
    Returns: list of dicts with sample, chrom, roh_bp, n_tracts, longest_tract, callable_bp, froh_chr
    """
    # Read callable spans per chrom
    callable_per_chr = defaultdict(int)
    with open(callable_bed) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            c = parts[0]
            callable_per_chr[c] += int(parts[2]) - int(parts[1])

    # Read tracts: chr start end sample [length]
    per_sample_chr = defaultdict(lambda: defaultdict(list))
    with open(tracts_bed) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            c, s, e, sample = parts[0], int(parts[1]), int(parts[2]), parts[3]
            length = e - s
            if length > 0:
                per_sample_chr[sample][c].append(length)

    rows = []
    for sample in sorted(per_sample_chr.keys()):
        for chrom in sorted(per_sample_chr[sample].keys()):
            tracts = per_sample_chr[sample][chrom]
            roh_bp = sum(tracts)
            n = len(tracts)
            longest = max(tracts) if tracts else 0
            cbp = callable_per_chr.get(chrom, 0)
            froh = roh_bp / cbp if cbp > 0 else 0
            rows.append({
                "sample": sample,
                "chrom": chrom,
                "roh_bp": roh_bp,
                "n_tracts": n,
                "longest_tract": longest,
                "callable_bp_chr": cbp,
                "froh_chr": froh,
            })
    return rows


def compute_het_in_out_roh(theta_dir, tracts_bed, samples, win_size=500000):
    """
    For each sample, intersect local theta windows with ROH BED to compute
    mean theta (diversity proxy) inside ROH vs outside ROH.

    Uses thetaStat window outputs already computed.
    Returns: dict of sample → {het_in_roh, het_out_roh, n_win_in, n_win_out}
    """
    # Read all ROH tracts per sample
    roh_per_sample = defaultdict(list)  # sample -> [(chr, start, end), ...]
    with open(tracts_bed) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            c, s, e, samp = parts[0], int(parts[1]), int(parts[2]), parts[3]
            roh_per_sample[samp].append((c, s, e))

    results = {}
    for sample in samples:
        # Find theta window file
        # Pattern: {sample}.win500000.step500000.pestPG
        # thetaStat do_stat output is named differently per ANGSD version
        theta_candidates = [
            os.path.join(theta_dir, f"{sample}.win{win_size}.step{win_size}.pestPG"),
            os.path.join(theta_dir, f"{sample}.win{win_size}.step{win_size}.gz.pestPG"),
        ]
        # Also try globbing
        import glob
        theta_glob = glob.glob(os.path.join(theta_dir, f"{sample}.win*step*.pestPG"))
        theta_candidates.extend(theta_glob)

        theta_file = None
        for tc in theta_candidates:
            if os.path.isfile(tc):
                theta_file = tc
                break

        if theta_file is None:
            results[sample] = {
                "het_proxy_in_roh": None,
                "het_proxy_out_roh": None,
                "n_win_in": 0,
                "n_win_out": 0,
            }
            continue

        # Read theta windows
        windows = []
        with open(theta_file) as f:
            header = f.readline().strip().split()
            # Find columns: (chr|Chromo), WinCenter, tW, tP, nSites
            chr_idx = next((i for i, h in enumerate(header) if h.lower() in ("chr", "chromo", "chrom")), None)
            center_idx = next((i for i, h in enumerate(header) if h.lower() in ("wincenter", "midpos")), None)
            tp_idx = next((i for i, h in enumerate(header) if h.lower() in ("tp", "thetapi", "pi")), None)
            nsites_idx = next((i for i, h in enumerate(header) if h.lower() in ("nsites", "nsite", "nsit")), None)

            if any(x is None for x in [chr_idx, center_idx, tp_idx]):
                # Try positional fallback for thetaStat do_stat standard output
                # Standard columns: (indexStart indexStop) Chromo WinCenter tW tP tF tH tL Tajima fuf fud fayh zeng nSites
                # or: Chr WinCenter tW tP ...
                if len(header) >= 5:
                    chr_idx = 0
                    center_idx = 1
                    tp_idx = 3  # tP
                    nsites_idx = len(header) - 1 if header[-1].lower().startswith("nsite") else None

            if chr_idx is None or center_idx is None or tp_idx is None:
                results[sample] = {
                    "het_proxy_in_roh": None,
                    "het_proxy_out_roh": None,
                    "n_win_in": 0,
                    "n_win_out": 0,
                }
                continue

            for line in f:
                parts = line.strip().split()
                if len(parts) <= max(chr_idx, center_idx, tp_idx):
                    continue
                try:
                    chrom = parts[chr_idx]
                    center = int(float(parts[center_idx]))
                    tp = float(parts[tp_idx])
                    nsites = int(parts[nsites_idx]) if nsites_idx is not None and nsites_idx < len(parts) else 1
                except (ValueError, IndexError):
                    continue
                if nsites > 0:
                    windows.append((chrom, center, tp / nsites, nsites))

        # Classify windows as in/out ROH
        roh_intervals = roh_per_sample.get(sample, [])
        in_vals = []
        out_vals = []
        half_win = win_size // 2

        for chrom, center, theta_per_site, nsites in windows:
            wstart = center - half_win
            wend = center + half_win
            in_roh = False
            for rc, rs, re in roh_intervals:
                if rc == chrom and rs < wend and re > wstart:
                    # Overlap
                    overlap = min(wend, re) - max(wstart, rs)
                    if overlap > half_win:  # majority overlap
                        in_roh = True
                        break
            if in_roh:
                in_vals.append(theta_per_site)
            else:
                out_vals.append(theta_per_site)

        results[sample] = {
            "het_proxy_in_roh": sum(in_vals) / len(in_vals) if in_vals else None,
            "het_proxy_out_roh": sum(out_vals) / len(out_vals) if out_vals else None,
            "n_win_in": len(in_vals),
            "n_win_out": len(out_vals),
        }

    return results


def main():
    ap = argparse.ArgumentParser(description="Parse ROH from ngsF-HMM, compute FROH, het in/out ROH")

    ap.add_argument("--ibd", required=True, help="Best .ibd file from ngsF-HMM")
    ap.add_argument("--pos", required=True, help="Position file (chrom\\tpos)")
    ap.add_argument("--ind", required=True, help="samples.ind file")
    ap.add_argument("--callable-bed", required=True, help="Callable BED mask")
    ap.add_argument("--repeats-bed", default="", help="Repeats BED (optional)")
    ap.add_argument("--het-summary", required=True, help="genomewide_heterozygosity.tsv")
    ap.add_argument("--theta-dir", required=True, help="Directory with thetaStat window outputs")
    ap.add_argument("--out-dir", required=True, help="Output directory for ROH summaries")
    ap.add_argument("--tables-dir", required=True, help="Final tables output directory")
    ap.add_argument("--convert-ibd-pl", required=True, help="Path to convert_ibd.pl")
    ap.add_argument("--summarize-roh-py", required=True, help="Path to summarize_ibd_roh.py")
    ap.add_argument("--win-size", type=int, default=500000, help="Window size for theta (must match thetaStat)")

    # Optional metadata for summarize_ibd_roh.py
    ap.add_argument("--q-matrix", default="", help="Q matrix file")
    ap.add_argument("--q-has-sample-col", action="store_true")
    ap.add_argument("--order-file", default="", help="Sample order file")
    ap.add_argument("--kinship", default="", help="Kinship/relatedness table")

    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.tables_dir, exist_ok=True)

    # ── Step 1: convert .ibd → BED tracts using convert_ibd.pl ────────────
    print("=== Step 1: Convert .ibd to BED tracts ===")
    tracts_bed = os.path.join(args.out_dir, "roh_tracts_all.bed")
    run_cmd_stdout(
        ["perl", args.convert_ibd_pl,
         "--ind", args.ind,
         "--pos", args.pos,
         "--ibd_pos", args.ibd],
        tracts_bed,
        "convert_ibd.pl"
    )
    n_tracts = sum(1 for _ in open(tracts_bed) if _.strip())
    print(f"  Total tracts: {n_tracts}")

    # ── Step 2: summarize_ibd_roh.py ──────────────────────────────────────
    print("\n=== Step 2: Summarize ROH/FROH ===")
    roh_prefix = os.path.join(args.out_dir, "catfish_roh")

    summarize_cmd = [
        "python3", args.summarize_roh_py,
        "--tracts-bed", tracts_bed,
        "--callable-bed", args.callable_bed,
        "--repeats-bed", args.repeats_bed if args.repeats_bed else args.callable_bed,
        "--use-bedtools-intersect",
        "--out-prefix", roh_prefix,
    ]
    if args.q_matrix:
        summarize_cmd += ["--q-matrix", args.q_matrix]
        if args.q_has_sample_col:
            summarize_cmd += ["--q-has-sample-col"]
    if args.order_file:
        summarize_cmd += ["--order-file", args.order_file]
    if args.kinship:
        summarize_cmd += ["--kinship", args.kinship]

    run_cmd(summarize_cmd, "summarize_ibd_roh.py")

    # ── Step 3: Per-chromosome ROH summary ────────────────────────────────
    print("\n=== Step 3: Per-chromosome ROH summary ===")
    chr_rows = compute_per_chr_roh(tracts_bed, args.callable_bed)

    per_chr_roh = os.path.join(args.tables_dir, "per_chr_roh_summary.tsv")
    with open(per_chr_roh, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "chrom", "roh_bp", "n_tracts",
                                           "longest_tract", "callable_bp_chr", "froh_chr"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(chr_rows)
    print(f"  Wrote {per_chr_roh} ({len(chr_rows)} rows)")

    # ── Step 4: Per-sample ROH details (from tracts BED) ─────────────────
    print("\n=== Step 4: Per-sample ROH detail table ===")
    sample_roh = defaultdict(lambda: {"roh_bp": 0, "n_tracts": 0, "longest": 0, "lengths": []})
    with open(tracts_bed) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            s_name = parts[3]
            length = int(parts[2]) - int(parts[1])
            if length <= 0:
                continue
            d = sample_roh[s_name]
            d["roh_bp"] += length
            d["n_tracts"] += 1
            d["lengths"].append(length)
            if length > d["longest"]:
                d["longest"] = length

    # Read callable total
    callable_bp = 0
    with open(args.callable_bed) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                callable_bp += int(parts[2]) - int(parts[1])

    per_sample_roh_tsv = os.path.join(args.tables_dir, "per_sample_roh.tsv")
    with open(per_sample_roh_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "roh_total_bp", "n_tracts", "longest_roh",
                     "mean_roh_length", "median_roh_length", "froh",
                     "callable_bp"])
        for sname in sorted(sample_roh.keys()):
            d = sample_roh[sname]
            lengths = sorted(d["lengths"])
            mean_l = sum(lengths) / len(lengths) if lengths else 0
            median_l = lengths[len(lengths) // 2] if lengths else 0
            froh = d["roh_bp"] / callable_bp if callable_bp > 0 else 0
            w.writerow([sname, d["roh_bp"], d["n_tracts"], d["longest"],
                        f"{mean_l:.0f}", f"{median_l:.0f}", f"{froh:.8f}",
                        callable_bp])
    print(f"  Wrote {per_sample_roh_tsv}")

    # ── Step 5: Het in/out ROH ────────────────────────────────────────────
    print("\n=== Step 5: Heterozygosity (diversity proxy) in/out ROH ===")
    samples = []
    with open(args.ind) as f:
        for line in f:
            s = line.strip().split()[0] if line.strip() else ""
            if s:
                samples.append(s)

    het_inout = compute_het_in_out_roh(args.theta_dir, tracts_bed, samples, args.win_size)

    het_inout_tsv = os.path.join(args.tables_dir, "per_sample_het_in_out_roh.tsv")
    with open(het_inout_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "theta_proxy_in_roh", "theta_proxy_out_roh",
                     "n_windows_in_roh", "n_windows_out_roh", "ratio_in_out"])
        for s in sorted(het_inout.keys()):
            d = het_inout[s]
            hin = d["het_proxy_in_roh"]
            hout = d["het_proxy_out_roh"]
            ratio = (hin / hout) if (hin is not None and hout is not None and hout > 0) else None
            w.writerow([
                s,
                f"{hin:.8e}" if hin is not None else "NA",
                f"{hout:.8e}" if hout is not None else "NA",
                d["n_win_in"],
                d["n_win_out"],
                f"{ratio:.6f}" if ratio is not None else "NA",
            ])
    print(f"  Wrote {het_inout_tsv}")

    # ── Step 6: Merge into master summary ─────────────────────────────────
    print("\n=== Step 6: Master summary table ===")
    het_gw = parse_het_summary(args.het_summary)

    master_tsv = os.path.join(args.tables_dir, "master_summary.tsv")
    with open(master_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "sample", "het_genomewide", "roh_total_bp", "froh",
            "n_tracts", "longest_roh", "mean_roh_length",
            "theta_proxy_in_roh", "theta_proxy_out_roh",
            "ratio_theta_in_out",
            "callable_bp"
        ])
        all_samples = sorted(set(list(het_gw.keys()) + list(sample_roh.keys())))
        for s in all_samples:
            hg = het_gw.get(s)
            sr = sample_roh.get(s, {"roh_bp": 0, "n_tracts": 0, "longest": 0, "lengths": []})
            hi = het_inout.get(s, {})
            froh = sr["roh_bp"] / callable_bp if callable_bp > 0 else 0
            mean_l = sum(sr["lengths"]) / len(sr["lengths"]) if sr["lengths"] else 0
            hin = hi.get("het_proxy_in_roh")
            hout = hi.get("het_proxy_out_roh")
            ratio = (hin / hout) if (hin is not None and hout is not None and hout > 0) else None

            w.writerow([
                s,
                f"{hg:.8e}" if hg is not None else "NA",
                sr["roh_bp"],
                f"{froh:.8f}",
                sr["n_tracts"],
                sr["longest"],
                f"{mean_l:.0f}",
                f"{hin:.8e}" if hin is not None else "NA",
                f"{hout:.8e}" if hout is not None else "NA",
                f"{ratio:.6f}" if ratio is not None else "NA",
                callable_bp,
            ])
    print(f"  Wrote {master_tsv}")

    # ── Copy key files to tables dir ──────────────────────────────────────
    import shutil
    for src_name in ["catfish_roh.per_sample_roh.tsv", "catfish_roh.per_sample_roh_bins_long.tsv"]:
        src = os.path.join(args.out_dir, src_name)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(args.tables_dir, src_name))

    # Copy tracts BED
    shutil.copy2(tracts_bed, os.path.join(args.tables_dir, "roh_tracts_all.bed"))

    print("\n=== Step 4 (parse_roh_and_het) complete ===")
    print(f"  ROH tracts: {tracts_bed}")
    print(f"  Per-sample ROH: {per_sample_roh_tsv}")
    print(f"  Per-chr ROH: {per_chr_roh}")
    print(f"  Het in/out ROH: {het_inout_tsv}")
    print(f"  Master summary: {master_tsv}")


if __name__ == "__main__":
    main()
