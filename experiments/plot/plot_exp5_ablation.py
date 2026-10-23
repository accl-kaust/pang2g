#!/usr/bin/env python3
"""
plot_exp5_ablation.py - Ablation study figure (Figure 5 in paper).

Shows the throughput contribution of each design decision.

Usage:
    python experiments/plot/plot_exp5_ablation.py \\
        --csv-main results/phase5/exp5/exp5_ablation_main_*.csv \\
        --csv-a4   results/phase5/exp5/exp5_ablation_anchor_*.csv \\
        --csv-a5   results/phase5/exp5/exp5_ablation_matching_*.csv \\
        --out      results/phase5/figures/
"""

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


ABLATION_LABELS = {
    "A1": "Cache-oblivious\nlayout",
    "A2": "Work-stealing\nscheduler",
    "A3": "Steal-largest\npolicy",
    "A4": "Minimizer\nanchor index",
    "A5": "Hungarian\nmatching",
}

VARIANT_DISPLAY = {
    "veb_layout":       ("vEB layout", "#2166ac"),
    "bfs_layout":       ("BFS layout", "#a6cee3"),
    "work_stealing":    ("Work-stealing", "#2166ac"),
    "static_partition": ("Static partition", "#a6cee3"),
    "steal_largest":    ("Steal-largest", "#2166ac"),
    "steal_random":     ("Steal-random", "#a6cee3"),
    "minimizer_index":  ("Minimizer index", "#2166ac"),
    "naive_anchor":     ("Naïve O(|VA||VB|)", "#a6cee3"),
    "hungarian":        ("Hungarian O(k³)", "#2166ac"),
    "brute_force":      ("Brute-force", "#a6cee3"),
}


def load_csvs(paths):
    dfs = []
    for p in paths:
        expanded = glob.glob(p) if "*" in p else [p]
        for f in expanded:
            try:
                dfs.append(pd.read_csv(f))
            except Exception as e:
                print(f"WARNING: {f}: {e}", file=sys.stderr)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def plot_ablation(df_main: pd.DataFrame, df_a4: pd.DataFrame,
                  df_a5: pd.DataFrame, out_dir: str) -> None:

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # -- Panel A: A1–A3 grouped bar chart ---------------------
    ax = axes[0]
    ablation_groups = ["A1", "A2", "A3"]
    x = np.arange(len(ablation_groups))
    width = 0.35

    bars_opt, bars_abl = [], []
    for i, ab in enumerate(ablation_groups):
        sub = df_main[df_main["ablation"] == ab]
        variants = sorted(sub["variant"].unique())
        if len(variants) < 2:
            continue

        # Production variant (first alphabetically that isn't "abl_")
        prod_v = [v for v in variants if not v.startswith("abl_") and
                  "bfs" not in v and "random" not in v and
                  "static" not in v and "naive" not in v and
                  "brute" not in v]
        abl_v = [v for v in variants if v not in prod_v]

        prod_ms = sub[sub["variant"] == prod_v[0]]["wall_ms"].mean() if prod_v else 0
        abl_ms  = sub[sub["variant"] == abl_v[0]]["wall_ms"].mean() if abl_v else 0

        prod_std = sub[sub["variant"] == prod_v[0]]["wall_ms"].std() if prod_v else 0
        abl_std  = sub[sub["variant"] == abl_v[0]]["wall_ms"].std() if abl_v else 0

        b1 = ax.bar(i - width/2, prod_ms/1000, width, color="#2166ac",
                    yerr=prod_std/1000, capsize=4)
        b2 = ax.bar(i + width/2, abl_ms/1000, width, color="#d6604d",
                    yerr=abl_std/1000, capsize=4)
        bars_opt.append(b1)
        bars_abl.append(b2)

        # Speedup annotation
        if abl_ms > 0:
            speedup = abl_ms / prod_ms
            ax.text(i, max(prod_ms, abl_ms)/1000 + 0.02, f"{speedup:.1f}×",
                    ha="center", fontsize=9, color="#d73027", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([ABLATION_LABELS[a] for a in ablation_groups], fontsize=9)
    ax.set_ylabel("Wall time (seconds)")
    ax.set_title("(A) Scheduler & Layout Ablation\n(n=10K snarls)")
    legend_elements = [
        mpatches.Patch(color="#2166ac", label="Production (optimized)"),
        mpatches.Patch(color="#d6604d", label="Ablated variant"),
    ]
    ax.legend(handles=legend_elements, fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")

    # -- Panel B: A4 - Anchor detection -----------------------
    ax = axes[1]
    if not df_a4.empty:
        for variant in df_a4["variant"].unique():
            sub = df_a4[df_a4["variant"] == variant].sort_values("n_snarls")
            label, color = VARIANT_DISPLAY.get(variant, (variant, "#888"))
            x_vals = sub.groupby("n_snarls")["wall_ms"].mean().index.values
            y_vals = sub.groupby("n_snarls")["wall_ms"].mean().values / 1000
            ax.plot(x_vals, y_vals, label=label, color=color,
                    marker="o", linewidth=2.5, markersize=7)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Graph size (snarls)")
        ax.set_ylabel("Wall time (seconds)")
        ax.set_title("(B) Anchor Index vs. Naïve Scan\n(A4)")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, which="both")

    # -- Panel C: A5 - Matching algorithm vs. k ---------------
    ax = axes[2]
    if not df_a5.empty:
        for variant in df_a5["variant"].unique():
            sub = df_a5[df_a5["variant"] == variant]
            sub = sub[sub["wall_ms"].apply(lambda x: str(x).isdigit())]
            sub["wall_ms"] = sub["wall_ms"].astype(int)
            if sub.empty: continue
            grp = sub.groupby("k")["wall_ms"].mean()
            label, color = VARIANT_DISPLAY.get(variant, (variant, "#888"))
            ax.plot(grp.index, grp.values, label=label, color=color,
                    marker="o", linewidth=2.5, markersize=7)

        ax.set_yscale("log")
        ax.set_xlabel("Alleles per snarl (k)")
        ax.set_ylabel("Wall time per run (ms)")
        ax.set_title("(C) Matching Algorithm vs. k\n(A5)")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, which="both")

        # Mark k=7 inflection point
        ax.axvline(x=7, color="orange", linestyle="--", alpha=0.7, label="k=7 inflection")
        ax.text(7.1, ax.get_ylim()[0] * 2, "k=7\ninflection",
                fontsize=8, color="darkorange", va="bottom")

    fig.suptitle("G2G Ablation Study - Experiment 5", fontsize=13, fontweight="bold")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        path = os.path.join(out_dir, f"fig5_ablation.{ext}")
        dpi = 150 if ext == "pdf" else 200
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        print(f"[plot] Saved: {path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv-main", nargs="+", default=[])
    parser.add_argument("--csv-a4",   nargs="+", default=[])
    parser.add_argument("--csv-a5",   nargs="+", default=[])
    parser.add_argument("--out", default="results/phase5/figures/")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    df_main = load_csvs(args.csv_main)
    df_a4   = load_csvs(args.csv_a4)
    df_a5   = load_csvs(args.csv_a5)
    print(f"[plot] Main: {len(df_main)} rows  A4: {len(df_a4)} rows  A5: {len(df_a5)} rows")
    plot_ablation(df_main, df_a4, df_a5, args.out)


if __name__ == "__main__":
    main()
