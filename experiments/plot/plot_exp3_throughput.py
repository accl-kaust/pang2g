#!/usr/bin/env python3
"""
plot_exp3_throughput.py — Throughput scaling figure (Figure 3 in paper).

Generates:
  - Log-log plot: wall time vs. graph size for all tools
  - Memory usage vs. graph size
  - Speedup of G2G-parallel vs. G2G-sequential
  - Combined 2×2 figure

Usage:
    python experiments/plot/plot_exp3_throughput.py \\
        --csv results/phase5/exp3/exp3_throughput_*.csv \\
        --out results/phase5/figures/
"""

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

TOOL_DISPLAY = {
    "g2g_parallel":   ("G2G (parallel)",    "#2166ac", "o",   2.5),
    "g2g_sequential": ("G2G (sequential)",  "#4393c3", "s",   2.0),
    "minimap2":       ("minimap2 asm5",     "#d6604d", "^",   2.0),
    "pangraph":       ("PanGraph",          "#878787", "D",   2.0),
}


def load_csv(paths: list[str]) -> pd.DataFrame:
    dfs = []
    for p in paths:
        try:
            dfs.append(pd.read_csv(p))
        except Exception as e:
            print(f"[plot] WARNING: could not load {p}: {e}", file=sys.stderr)
    if not dfs:
        raise ValueError("No CSV files loaded")
    df = pd.concat(dfs, ignore_index=True)
    df = df[df["exit_code"] == 0]
    return df


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["n_snarls", "tool"])
        .agg(
            wall_ms_mean=("wall_ms", "mean"),
            wall_ms_std=("wall_ms", "std"),
            wall_ms_min=("wall_ms", "min"),
            wall_ms_max=("wall_ms", "max"),
            peak_kb_mean=("peak_kb", "mean"),
            n=("wall_ms", "count"),
        )
        .reset_index()
    )


def plot_throughput(agg: pd.DataFrame, out_dir: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # ── Panel A: Wall time vs. n_snarls (log-log) ──
    ax = axes[0, 0]
    for tool, (label, color, marker, lw) in TOOL_DISPLAY.items():
        sub = agg[agg["tool"] == tool].sort_values("n_snarls")
        if sub.empty:
            continue
        x = sub["n_snarls"].values
        y = sub["wall_ms_mean"].values / 1000  # → seconds
        yerr = sub["wall_ms_std"].values / 1000
        ax.errorbar(x, y, yerr=np.nan_to_num(yerr),
                    label=label, color=color, marker=marker,
                    linewidth=lw, markersize=7, capsize=3)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Graph size (snarls)")
    ax.set_ylabel("Wall time (seconds)")
    ax.set_title("(A) Throughput scaling")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which="both")

    # Fit power-law slope for G2G parallel
    sub = agg[agg["tool"] == "g2g_parallel"].sort_values("n_snarls")
    if len(sub) >= 3:
        log_x = np.log10(sub["n_snarls"].values)
        log_y = np.log10(sub["wall_ms_mean"].values)
        slope, intercept = np.polyfit(log_x, log_y, 1)
        ax.text(0.05, 0.95, f"G2G slope: O(n^{slope:.2f})",
                transform=ax.transAxes, fontsize=9, va="top",
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    # ── Panel B: Memory usage vs. n_snarls ────────
    ax = axes[0, 1]
    for tool, (label, color, marker, lw) in TOOL_DISPLAY.items():
        sub = agg[agg["tool"] == tool].sort_values("n_snarls")
        if sub.empty: continue
        x = sub["n_snarls"].values
        y = sub["peak_kb_mean"].values / (1024 ** 2)  # → GB
        ax.plot(x, y, label=label, color=color, marker=marker,
                linewidth=lw, markersize=7)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Graph size (snarls)")
    ax.set_ylabel("Peak memory (GB)")
    ax.set_title("(B) Memory usage")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which="both")

    # ── Panel C: Speedup G2G-par vs. G2G-seq ─────
    ax = axes[1, 0]
    par = agg[agg["tool"] == "g2g_parallel"].set_index("n_snarls")
    seq = agg[agg["tool"] == "g2g_sequential"].set_index("n_snarls")
    common = par.index.intersection(seq.index)
    if len(common) > 0:
        speedup = seq.loc[common, "wall_ms_mean"] / par.loc[common, "wall_ms_mean"]
        ax.plot(common, speedup.values, color="#2166ac", marker="o",
                linewidth=2, markersize=7, label="G2G parallel speedup")
        ax.axhline(y=1, color="gray", linestyle="--", alpha=0.5)
        ax.set_xscale("log")
        ax.set_xlabel("Graph size (snarls)")
        ax.set_ylabel("Speedup vs. sequential")
        ax.set_title("(C) Speedup over sequential G2G")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    # ── Panel D: G2G vs. minimap2 speedup ─────────
    ax = axes[1, 1]
    mm2 = agg[agg["tool"] == "minimap2"].set_index("n_snarls")
    par = agg[agg["tool"] == "g2g_parallel"].set_index("n_snarls")
    common = mm2.index.intersection(par.index)
    if len(common) > 0:
        ratio = mm2.loc[common, "wall_ms_mean"] / par.loc[common, "wall_ms_mean"]
        colors = ["#d73027" if r < 1 else "#1a9850" for r in ratio]
        bars = ax.bar(range(len(common)), ratio.values, color=colors)
        ax.axhline(y=1, color="gray", linestyle="--", alpha=0.7)
        ax.set_xticks(range(len(common)))
        ax.set_xticklabels([f"{n/1000:.0f}K" for n in common], rotation=45, ha="right")
        ax.set_xlabel("Graph size (snarls)")
        ax.set_ylabel("minimap2 time / G2G time")
        ax.set_title("(D) G2G vs. minimap2 time ratio\n(>1 means G2G is faster)")
        ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle("G2G Throughput Scaling — Experiment 3", fontsize=13, fontweight="bold")
    plt.tight_layout()

    out_path = os.path.join(out_dir, "fig3_throughput_scaling.pdf")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"[plot] Saved: {out_path}")

    png_path = out_path.replace(".pdf", ".png")
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    print(f"[plot] Saved: {png_path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", nargs="+", required=True,
                        help="Path(s) to exp3 CSV files (glob supported)")
    parser.add_argument("--out", default="results/phase5/figures/",
                        help="Output directory")
    args = parser.parse_args()

    # Expand globs
    csv_paths = []
    for pat in args.csv:
        expanded = glob.glob(pat)
        csv_paths.extend(expanded if expanded else [pat])

    os.makedirs(args.out, exist_ok=True)
    print(f"[plot] Loading {len(csv_paths)} CSV file(s)...")

    df = load_csv(csv_paths)
    print(f"[plot] {len(df)} rows loaded  tools={df['tool'].unique().tolist()}")

    agg = aggregate(df)
    plot_throughput(agg, args.out)

    # Also save aggregated table
    table_path = os.path.join(args.out, "exp3_summary.csv")
    agg.to_csv(table_path, index=False)
    print(f"[plot] Summary table: {table_path}")


if __name__ == "__main__":
    main()
