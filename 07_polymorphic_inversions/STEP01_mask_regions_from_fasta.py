#!/usr/bin/env python3
"""
STEP01_mask_regions_from_fasta.py

Extract BED intervals + per-contig stats for:
  1) hard-masked regions:            N (and n)
  2) soft-masked regions:            a/c/g/t
  3) normal uppercase regions:       A/C/G/T
  4) soft+hard masked regions:       a/c/g/t + N/n
  5) inversion-callable regions:     A/C/G/T + a/c/g/t   (all ACGT, any case)

Outputs (0-based, half-open BED):
  <prefix>.hardN.bed
  <prefix>.softacgt.bed
  <prefix>.normalACGT.bed
  <prefix>.masked_acgtN.bed
  <prefix>.inversion_acgt_allcase.bed
  <prefix>.stats.tsv
"""

import argparse
from dataclasses import dataclass
from typing import Optional, TextIO


@dataclass
class RunState:
    start: Optional[int] = None


@dataclass
class TrackStats:
    bp: int = 0
    n_intervals: int = 0
    max_len: int = 0


def close_run(chrom: str, pos: int, state: RunState, out: TextIO, label: str, tstats: TrackStats):
    if state.start is None:
        return
    s = state.start
    e = pos
    if e > s:
        out.write(f"{chrom}\t{s}\t{e}\t{label}\n")
        length = e - s
        tstats.n_intervals += 1
        tstats.bp += length
        if length > tstats.max_len:
            tstats.max_len = length
    state.start = None


def open_run(pos: int, state: RunState):
    if state.start is None:
        state.start = pos


def flush_contig(
    chrom, pos,
    hard_state, soft_state, norm_state, mask_state, inversion_state,
    hard_out, soft_out, norm_out, mask_out, inversion_out,
    hard_t, soft_t, norm_t, mask_t, inversion_t,
):
    close_run(chrom, pos, hard_state, hard_out, "hardN", hard_t)
    close_run(chrom, pos, soft_state, soft_out, "softacgt", soft_t)
    close_run(chrom, pos, norm_state, norm_out, "normalACGT", norm_t)
    close_run(chrom, pos, mask_state, mask_out, "masked_acgtN", mask_t)
    close_run(chrom, pos, inversion_state, inversion_out, "inversion_acgt_allcase", inversion_t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--fasta", required=True, help="Input FASTA")
    ap.add_argument("-p", "--prefix", required=True, help="Output prefix")
    args = ap.parse_args()

    hard_bed_path = f"{args.prefix}.hardN.bed"
    soft_bed_path = f"{args.prefix}.softacgt.bed"
    norm_bed_path = f"{args.prefix}.normalACGT.bed"
    mask_bed_path = f"{args.prefix}.masked_acgtN.bed"
    inversion_bed_path = f"{args.prefix}.inversion_acgt_allcase.bed"
    stats_path = f"{args.prefix}.stats.tsv"

    HARD = set("Nn")
    SOFT = set("acgt")
    NORM = set("ACGT")
    INVERSION = set("ACGTacgt")

    with open(hard_bed_path, "w") as hard_out, \
         open(soft_bed_path, "w") as soft_out, \
         open(norm_bed_path, "w") as norm_out, \
         open(mask_bed_path, "w") as mask_out, \
         open(inversion_bed_path, "w") as inversion_out, \
         open(stats_path, "w") as stats_out:

        stats_out.write(
            "contig\tlength_bp\t"
            "bp_hardN\tbp_softacgt\tbp_normalACGT\tbp_masked_acgtN\tbp_inversion_acgt_allcase\tbp_other\t"
            "pct_hardN\tpct_softacgt\tpct_normalACGT\tpct_masked_acgtN\tpct_inversion_acgt_allcase\tpct_other\t"
            "nint_hardN\tnint_softacgt\tnint_normalACGT\tnint_masked_acgtN\tnint_inversion_acgt_allcase\t"
            "maxrun_hardN\tmaxrun_softacgt\tmaxrun_normalACGT\tmaxrun_masked_acgtN\tmaxrun_inversion_acgt_allcase\n"
        )

        hard_state = RunState()
        soft_state = RunState()
        norm_state = RunState()
        mask_state = RunState()
        inversion_state = RunState()

        chrom = None
        pos = 0
        bp_hard = bp_soft = bp_norm = bp_mask = bp_inversion = bp_other = 0

        hard_t = TrackStats()
        soft_t = TrackStats()
        norm_t = TrackStats()
        mask_t = TrackStats()
        inversion_t = TrackStats()

        tot_len = 0
        tot_hard = tot_soft = tot_norm = tot_mask = tot_inversion = tot_other = 0

        tot_hard_t = TrackStats()
        tot_soft_t = TrackStats()
        tot_norm_t = TrackStats()
        tot_mask_t = TrackStats()
        tot_inversion_t = TrackStats()

        def write_contig_stats(contig_name, length, bh, bs, bn, bm, bi, bo, ht, st, nt, mt, it):
            if length == 0:
                return

            def pct(x):
                return (100.0 * x / length) if length else 0.0

            stats_out.write(
                f"{contig_name}\t{length}\t"
                f"{bh}\t{bs}\t{bn}\t{bm}\t{bi}\t{bo}\t"
                f"{pct(bh):.6f}\t{pct(bs):.6f}\t{pct(bn):.6f}\t{pct(bm):.6f}\t{pct(bi):.6f}\t{pct(bo):.6f}\t"
                f"{ht.n_intervals}\t{st.n_intervals}\t{nt.n_intervals}\t{mt.n_intervals}\t{it.n_intervals}\t"
                f"{ht.max_len}\t{st.max_len}\t{nt.max_len}\t{mt.max_len}\t{it.max_len}\n"
            )

        with open(args.fasta, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue

                if line.startswith(">"):
                    if chrom is not None:
                        flush_contig(
                            chrom, pos,
                            hard_state, soft_state, norm_state, mask_state, inversion_state,
                            hard_out, soft_out, norm_out, mask_out, inversion_out,
                            hard_t, soft_t, norm_t, mask_t, inversion_t
                        )
                        write_contig_stats(
                            chrom, pos,
                            bp_hard, bp_soft, bp_norm, bp_mask, bp_inversion, bp_other,
                            hard_t, soft_t, norm_t, mask_t, inversion_t
                        )
                        tot_len += pos
                        tot_hard += bp_hard; tot_soft += bp_soft; tot_norm += bp_norm
                        tot_mask += bp_mask; tot_inversion += bp_inversion; tot_other += bp_other
                        for src, dst in [(hard_t, tot_hard_t), (soft_t, tot_soft_t),
                                         (norm_t, tot_norm_t), (mask_t, tot_mask_t),
                                         (inversion_t, tot_inversion_t)]:
                            dst.bp += src.bp
                            dst.n_intervals += src.n_intervals
                            dst.max_len = max(dst.max_len, src.max_len)

                    chrom = line[1:].split()[0]
                    pos = 0
                    hard_state = RunState(); soft_state = RunState(); norm_state = RunState()
                    mask_state = RunState(); inversion_state = RunState()
                    bp_hard = bp_soft = bp_norm = bp_mask = bp_inversion = bp_other = 0
                    hard_t = TrackStats(); soft_t = TrackStats(); norm_t = TrackStats()
                    mask_t = TrackStats(); inversion_t = TrackStats()
                    continue

                for c in line:
                    in_hard = c in HARD
                    in_soft = c in SOFT
                    in_norm = c in NORM
                    in_mask = in_hard or in_soft
                    in_inversion = c in INVERSION

                    if in_hard:     bp_hard += 1
                    elif in_soft:   bp_soft += 1
                    elif in_norm:   bp_norm += 1
                    else:           bp_other += 1
                    if in_mask:     bp_mask += 1
                    if in_inversion: bp_inversion += 1

                    if in_hard:     open_run(pos, hard_state)
                    else:           close_run(chrom, pos, hard_state, hard_out, "hardN", hard_t)
                    if in_soft:     open_run(pos, soft_state)
                    else:           close_run(chrom, pos, soft_state, soft_out, "softacgt", soft_t)
                    if in_norm:     open_run(pos, norm_state)
                    else:           close_run(chrom, pos, norm_state, norm_out, "normalACGT", norm_t)
                    if in_mask:     open_run(pos, mask_state)
                    else:           close_run(chrom, pos, mask_state, mask_out, "masked_acgtN", mask_t)
                    if in_inversion: open_run(pos, inversion_state)
                    else:           close_run(chrom, pos, inversion_state, inversion_out, "inversion_acgt_allcase", inversion_t)

                    pos += 1

        # flush last contig
        if chrom is not None:
            flush_contig(
                chrom, pos,
                hard_state, soft_state, norm_state, mask_state, inversion_state,
                hard_out, soft_out, norm_out, mask_out, inversion_out,
                hard_t, soft_t, norm_t, mask_t, inversion_t
            )
            write_contig_stats(
                chrom, pos,
                bp_hard, bp_soft, bp_norm, bp_mask, bp_inversion, bp_other,
                hard_t, soft_t, norm_t, mask_t, inversion_t
            )
            tot_len += pos
            tot_hard += bp_hard; tot_soft += bp_soft; tot_norm += bp_norm
            tot_mask += bp_mask; tot_inversion += bp_inversion; tot_other += bp_other
            for src, dst in [(hard_t, tot_hard_t), (soft_t, tot_soft_t),
                             (norm_t, tot_norm_t), (mask_t, tot_mask_t),
                             (inversion_t, tot_inversion_t)]:
                dst.bp += src.bp
                dst.n_intervals += src.n_intervals
                dst.max_len = max(dst.max_len, src.max_len)

        def pct_total(x):
            return (100.0 * x / tot_len) if tot_len else 0.0

        stats_out.write(
            f"__TOTAL__\t{tot_len}\t"
            f"{tot_hard}\t{tot_soft}\t{tot_norm}\t{tot_mask}\t{tot_inversion}\t{tot_other}\t"
            f"{pct_total(tot_hard):.6f}\t{pct_total(tot_soft):.6f}\t{pct_total(tot_norm):.6f}\t"
            f"{pct_total(tot_mask):.6f}\t{pct_total(tot_inversion):.6f}\t{pct_total(tot_other):.6f}\t"
            f"{tot_hard_t.n_intervals}\t{tot_soft_t.n_intervals}\t{tot_norm_t.n_intervals}\t"
            f"{tot_mask_t.n_intervals}\t{tot_inversion_t.n_intervals}\t"
            f"{tot_hard_t.max_len}\t{tot_soft_t.max_len}\t{tot_norm_t.max_len}\t"
            f"{tot_mask_t.max_len}\t{tot_inversion_t.max_len}\n"
        )

    print("Wrote:")
    print(" ", hard_bed_path)
    print(" ", soft_bed_path)
    print(" ", norm_bed_path)
    print(" ", mask_bed_path)
    print(" ", inversion_bed_path)
    print(" ", stats_path)


if __name__ == "__main__":
    main()
