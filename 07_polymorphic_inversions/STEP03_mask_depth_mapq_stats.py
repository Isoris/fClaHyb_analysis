#!/usr/bin/env python3
"""
STEP03_mask_depth_mapq_stats.py

Compute stats for one or more BED masks:
  1) mask size in bp
  2) percent of genome covered
  3) number of intervals
  4) mean interval length
  5) depth stats over those regions from BAMs
  6) depth stats over those regions using only reads with MAPQ >= chosen threshold
"""

import argparse
import os
import statistics
import subprocess


def run_cmd(cmd):
    p = subprocess.run(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
        )
    return p.stdout


def read_fai_total_bp(fai_path):
    total = 0
    with open(fai_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            total += int(fields[1])
    return total


def read_bed_stats(bed_path):
    total_bp = 0
    n_intervals = 0
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            start = int(fields[1])
            end = int(fields[2])
            length = end - start
            if length <= 0:
                continue
            total_bp += length
            n_intervals += 1
    mean_len = (total_bp / n_intervals) if n_intervals else 0.0
    return total_bp, n_intervals, mean_len


def read_bamlist(bamlist_path):
    bams = []
    with open(bamlist_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            bams.append(line.split()[0])
    return bams


def parse_bedcov_output(text):
    depths = []
    for line in text.strip().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        start = int(fields[1])
        end = int(fields[2])
        covered_bases = float(fields[3])
        length = end - start
        if length <= 0:
            continue
        depths.append(covered_bases / length)
    return depths


def summarize(values):
    if not values:
        return {"mean": "NA", "median": "NA", "sd": "NA", "min": "NA", "max": "NA"}
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    return {
        "mean": f"{statistics.mean(values):.6f}",
        "median": f"{statistics.median(values):.6f}",
        "sd": f"{sd:.6f}",
        "min": f"{min(values):.6f}",
        "max": f"{max(values):.6f}",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fai", help="Reference FASTA .fai file")
    ap.add_argument("--genome-size", type=int, help="Total genome size in bp")
    ap.add_argument("--bamlist", required=True, help="One BAM path per line")
    ap.add_argument("--mapq", type=int, default=60, help="MAPQ threshold")
    ap.add_argument("--mask", action="append", required=True, help="Mask definition as name=path/to/file.bed")
    ap.add_argument("--out", required=True, help="Output TSV")
    args = ap.parse_args()

    if args.fai:
        genome_bp = read_fai_total_bp(args.fai)
    elif args.genome_size:
        genome_bp = args.genome_size
    else:
        raise SystemExit("ERROR: provide either --fai or --genome-size")

    bams = read_bamlist(args.bamlist)
    if not bams:
        raise SystemExit("ERROR: no BAMs found in bamlist")

    mask_defs = []
    for item in args.mask:
        if "=" not in item:
            raise SystemExit(f"ERROR: --mask must be name=path, got: {item}")
        name, path = item.split("=", 1)
        if not os.path.exists(path):
            raise SystemExit(f"ERROR: mask BED not found: {path}")
        mask_defs.append((name, path))

    header = [
        "mask_name", "mask_bed", "genome_bp", "mask_bp", "pct_genome",
        "n_intervals", "mean_interval_bp", "n_bams",
        "depth_all_mean_of_samples", "depth_all_median_of_samples",
        "depth_all_sd_of_samples", "depth_all_min_of_samples", "depth_all_max_of_samples",
        f"depth_mapq{args.mapq}_mean_of_samples", f"depth_mapq{args.mapq}_median_of_samples",
        f"depth_mapq{args.mapq}_sd_of_samples", f"depth_mapq{args.mapq}_min_of_samples",
        f"depth_mapq{args.mapq}_max_of_samples",
    ]

    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")

        for mask_name, bed_path in mask_defs:
            mask_bp, n_intervals, mean_interval_bp = read_bed_stats(bed_path)
            pct_genome = (100.0 * mask_bp / genome_bp) if genome_bp else float("nan")

            depth_all_per_sample = []
            depth_mapq_per_sample = []

            for bam in bams:
                out_all = run_cmd(["samtools", "bedcov", bed_path, bam])
                depths_all = parse_bedcov_output(out_all)
                depth_all_per_sample.append(statistics.mean(depths_all) if depths_all else 0.0)

                cmd_mapq = f"samtools view -b -q {args.mapq} {bam} | samtools bedcov {bed_path} -"
                p = subprocess.run(cmd_mapq, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if p.returncode != 0:
                    raise RuntimeError(f"Command failed: {cmd_mapq}\nSTDERR:\n{p.stderr}")
                depths_mapq = parse_bedcov_output(p.stdout)
                depth_mapq_per_sample.append(statistics.mean(depths_mapq) if depths_mapq else 0.0)

            s_all = summarize(depth_all_per_sample)
            s_mapq = summarize(depth_mapq_per_sample)

            row = [
                mask_name, bed_path, str(genome_bp), str(mask_bp), f"{pct_genome:.6f}",
                str(n_intervals), f"{mean_interval_bp:.2f}", str(len(bams)),
                s_all["mean"], s_all["median"], s_all["sd"], s_all["min"], s_all["max"],
                s_mapq["mean"], s_mapq["median"], s_mapq["sd"], s_mapq["min"], s_mapq["max"],
            ]
            out.write("\t".join(row) + "\n")

    print(f"[INFO] Wrote {args.out}")


if __name__ == "__main__":
    main()
