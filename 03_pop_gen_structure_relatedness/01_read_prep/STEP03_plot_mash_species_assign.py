#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per_cga", required=True, help="results.per_CGA.tsv")
    ap.add_argument("--per_run", default=None, help="results.tsv (optional)")
    ap.add_argument("--outprefix", default="mash_species_assign", help="output prefix")
    ap.add_argument("--margin", type=float, default=0.002, help="same as your MARGIN (for reference lines)")
    ap.add_argument("--label_top", type=int, default=15, help="label N most suspicious CGAs on scatter")
    args = ap.parse_args()

    df = pd.read_csv(args.per_cga, sep="\t")
    # expected columns:
    # CGA_sample, n_runs, median_dist_gar, median_dist_mac, majority_call, calls_breakdown

    df["median_dist_gar"] = pd.to_numeric(df["median_dist_gar"], errors="coerce")
    df["median_dist_mac"] = pd.to_numeric(df["median_dist_mac"], errors="coerce")
    df["delta_mac_minus_gar"] = df["median_dist_mac"] - df["median_dist_gar"]
    df["absdiff"] = (df["median_dist_mac"] - df["median_dist_gar"]).abs()

    # ---- 1) Scatter: dist_gar vs dist_mac ----
    x = df["median_dist_gar"].to_numpy()
    y = df["median_dist_mac"].to_numpy()

    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]; y = y[finite]
    df_sc = df.loc[finite].copy()

    plt.figure()
    plt.scatter(x, y, s=18)
    # y=x diagonal
    lo = min(x.min(), y.min())
    hi = max(x.max(), y.max())
    plt.plot([lo, hi], [lo, hi])
    plt.xlabel("median_dist_gar")
    plt.ylabel("median_dist_mac")
    plt.title("Mash distance: Gar vs Mac (per CGA)")

    # label the most suspicious = smallest delta (closest to diagonal / wrong)
    df_sc = df_sc.sort_values("delta_mac_minus_gar", ascending=True)
    top = df_sc.head(args.label_top)
    for _, r in top.iterrows():
        plt.annotate(
            r["CGA_sample"],
            (r["median_dist_gar"], r["median_dist_mac"]),
            fontsize=8,
            xytext=(3, 3),
            textcoords="offset points",
        )

    plt.tight_layout()
    plt.savefig(f"{args.outprefix}.scatter_perCGA.png", dpi=200)
    plt.savefig(f"{args.outprefix}.scatter_perCGA.pdf")
    plt.close()

    # ---- 2) Histogram of delta ----
    plt.figure()
    d = df["delta_mac_minus_gar"].to_numpy()
    d = d[np.isfinite(d)]
    plt.hist(d, bins=40)
    plt.xlabel("delta = median_dist_mac - median_dist_gar")
    plt.ylabel("count")
    plt.title("Separation between references (per CGA)")
    plt.tight_layout()
    plt.savefig(f"{args.outprefix}.hist_delta_perCGA.png", dpi=200)
    plt.savefig(f"{args.outprefix}.hist_delta_perCGA.pdf")
    plt.close()

    # ---- 3) Ranked “most suspicious” table + plot ----
    sus = df.sort_values("delta_mac_minus_gar", ascending=True).copy()
    sus_cols = ["CGA_sample","n_runs","median_dist_gar","median_dist_mac","delta_mac_minus_gar","majority_call","calls_breakdown"]
    sus[sus_cols].to_csv(f"{args.outprefix}.most_suspicious_perCGA.tsv", sep="\t", index=False)

    topn = min(30, len(sus))
    plt.figure(figsize=(10, 6))
    plt.plot(np.arange(topn), sus["delta_mac_minus_gar"].head(topn).to_numpy(), marker="o", linestyle="-")
    plt.xticks(np.arange(topn), sus["CGA_sample"].head(topn).to_list(), rotation=90)
    plt.ylabel("delta = dist_mac - dist_gar (smaller = more suspicious)")
    plt.title(f"Top {topn} most suspicious CGAs by delta")
    plt.tight_layout()
    plt.savefig(f"{args.outprefix}.top_suspicious_perCGA.png", dpi=200)
    plt.savefig(f"{args.outprefix}.top_suspicious_perCGA.pdf")
    plt.close()

    # ---- Optional: per-run scatter if provided ----
    if args.per_run:
        dr = pd.read_csv(args.per_run, sep="\t")
        dr["dist_gar"] = pd.to_numeric(dr["dist_gar"], errors="coerce")
        dr["dist_mac"] = pd.to_numeric(dr["dist_mac"], errors="coerce")
        dr["delta_mac_minus_gar"] = dr["dist_mac"] - dr["dist_gar"]

        xr = dr["dist_gar"].to_numpy()
        yr = dr["dist_mac"].to_numpy()
        fr = np.isfinite(xr) & np.isfinite(yr)
        xr = xr[fr]; yr = yr[fr]

        plt.figure()
        plt.scatter(xr, yr, s=10)
        lo = min(xr.min(), yr.min())
        hi = max(xr.max(), yr.max())
        plt.plot([lo, hi], [lo, hi])
        plt.xlabel("dist_gar (per run)")
        plt.ylabel("dist_mac (per run)")
        plt.title("Mash distance: Gar vs Mac (per run-lane unit)")
        plt.tight_layout()
        plt.savefig(f"{args.outprefix}.scatter_perRUN.png", dpi=200)
        plt.savefig(f"{args.outprefix}.scatter_perRUN.pdf")
        plt.close()

        # export most suspicious runs too
        dr.sort_values("delta_mac_minus_gar", ascending=True).to_csv(
            f"{args.outprefix}.most_suspicious_perRUN.tsv", sep="\t", index=False
        )

    # quick console summary
    print("Wrote:")
    print(f"  {args.outprefix}.scatter_perCGA.png/.pdf")
    print(f"  {args.outprefix}.hist_delta_perCGA.png/.pdf")
    print(f"  {args.outprefix}.top_suspicious_perCGA.png/.pdf")
    print(f"  {args.outprefix}.most_suspicious_perCGA.tsv")
    if args.per_run:
        print(f"  {args.outprefix}.scatter_perRUN.png/.pdf")
        print(f"  {args.outprefix}.most_suspicious_perRUN.tsv")

if __name__ == "__main__":
    main()
