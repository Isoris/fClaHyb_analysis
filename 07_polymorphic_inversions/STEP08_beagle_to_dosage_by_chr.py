#!/usr/bin/env python3
"""
STEP08_beagle_to_dosage_by_chr.py

Read an ANGSD .beagle.gz file, compute expected genotype dosage per sample,
and split output by chromosome.

For each site:
  dosage = P(AB) + 2 * P(BB)

Outputs per chromosome:
  <outdir>/<chrom>.sites.tsv.gz
  <outdir>/<chrom>.dosage.tsv.gz
"""

import argparse
import gzip
import os
import re
import sys
from typing import Dict, List, Tuple, TextIO


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beagle", required=True, help="Input ANGSD .beagle.gz file")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--digits", type=int, default=6, help="Decimal places")
    ap.add_argument("--prefix", default="", help="Optional output prefix")
    return ap.parse_args()


def infer_chr_pos(marker: str) -> Tuple[str, int]:
    marker2 = marker.replace(":", "_")
    parts = marker2.split("_")
    if len(parts) < 2:
        raise ValueError(f"Could not parse marker '{marker}'")

    pos_idx = None
    for i in range(len(parts) - 1, -1, -1):
        if re.fullmatch(r"\d+", parts[i]):
            pos_idx = i
            break

    if pos_idx is None:
        raise ValueError(f"Could not find numeric position in marker '{marker}'")

    pos = int(parts[pos_idx])
    chrom_parts = parts[:pos_idx]
    if not chrom_parts:
        raise ValueError(f"Could not infer chromosome from marker '{marker}'")

    chrom = "_".join(chrom_parts)
    return chrom, pos


def open_gzip_text_writer(path: str) -> TextIO:
    return gzip.open(path, "wt")


class ChromWriters:
    def __init__(self, outdir: str, prefix: str, samples: List[str]):
        self.outdir = outdir
        self.prefix = prefix
        self.samples = samples
        self.handles: Dict[str, Tuple[TextIO, TextIO]] = {}

    def _make_paths(self, chrom: str) -> Tuple[str, str]:
        base = f"{self.prefix}{chrom}" if self.prefix else chrom
        sites_path = os.path.join(self.outdir, f"{base}.sites.tsv.gz")
        dosage_path = os.path.join(self.outdir, f"{base}.dosage.tsv.gz")
        return sites_path, dosage_path

    def get(self, chrom: str) -> Tuple[TextIO, TextIO]:
        if chrom not in self.handles:
            sites_path, dosage_path = self._make_paths(chrom)
            sites_fh = open_gzip_text_writer(sites_path)
            dosage_fh = open_gzip_text_writer(dosage_path)

            sites_fh.write("marker\tchrom\tpos\tallele1\tallele2\n")
            dosage_fh.write("marker\t" + "\t".join(self.samples) + "\n")

            self.handles[chrom] = (sites_fh, dosage_fh)
        return self.handles[chrom]

    def close_all(self) -> None:
        for sites_fh, dosage_fh in self.handles.values():
            sites_fh.close()
            dosage_fh.close()


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    n_sites = 0
    chrom_counts: Dict[str, int] = {}

    with gzip.open(args.beagle, "rt") as fh:
        header = fh.readline().rstrip("\n").split()
        if len(header) < 6:
            raise RuntimeError("BEAGLE header too short")

        if (len(header) - 3) % 3 != 0:
            raise RuntimeError("Header genotype-likelihood columns after first 3 fields are not divisible by 3")

        raw_sample_cols = header[3:]
        n_samples = len(raw_sample_cols) // 3

        samples: List[str] = []
        for i in range(n_samples):
            s0 = raw_sample_cols[i * 3]
            s1 = raw_sample_cols[i * 3 + 1]
            s2 = raw_sample_cols[i * 3 + 2]
            if not (s0 == s1 == s2):
                raise RuntimeError(f"Expected repeated sample names in groups of 3, got: {s0}, {s1}, {s2}")
            samples.append(s0)

        writers = ChromWriters(args.outdir, args.prefix, samples)
        fmt = "{:." + str(args.digits) + "f}"

        try:
            for line_num, line in enumerate(fh, start=2):
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split()

                expected_len = 3 + 3 * n_samples
                if len(fields) != expected_len:
                    raise RuntimeError(f"Line {line_num}: expected {expected_len} columns, found {len(fields)}")

                marker = fields[0]
                allele1 = fields[1]
                allele2 = fields[2]

                chrom, pos = infer_chr_pos(marker)
                sites_fh, dosage_fh = writers.get(chrom)

                sites_fh.write(f"{marker}\t{chrom}\t{pos}\t{allele1}\t{allele2}\n")

                dosages: List[str] = [marker]
                gl = fields[3:]

                for i in range(n_samples):
                    p_ab = float(gl[i * 3 + 1])
                    p_bb = float(gl[i * 3 + 2])
                    dosage = p_ab + 2.0 * p_bb
                    dosages.append(fmt.format(dosage))

                dosage_fh.write("\t".join(dosages) + "\n")

                n_sites += 1
                chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

        finally:
            writers.close_all()

    sys.stderr.write(f"[INFO] Samples: {n_samples}\n")
    sys.stderr.write(f"[INFO] Sites processed: {n_sites}\n")
    for chrom in sorted(chrom_counts):
        sys.stderr.write(f"[INFO] {chrom}: {chrom_counts[chrom]} sites\n")


if __name__ == "__main__":
    main()
