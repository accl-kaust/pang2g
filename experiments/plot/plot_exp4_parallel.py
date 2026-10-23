#!/usr/bin/env python3
"""
plot_exp4_parallel.py - Parallel speedup figure (Figure 4 in paper).

Generates:
  - Speedup vs. thread count (with ideal linear line)
  - Parallel efficiency vs. thread count
  - Memory vs. thread count
  - Amdahl's law fit

Usage:
    python experiments/plot/plot_exp4_parallel.py \\
        --csv results/phase5/exp4/exp4_parallel_scaling_*.csv \\
        --out results/phase5/figures/
"""

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def load_csv(paths):
    dfs = []
    for p in paths:
        try:
            dfs.append(pd.read_csv(p))
        except Exception as e:
            print(f"WARNING: could not load {p}: {e}", file=sys.stderr)
    if not dfs:
        raise ValueError("No CSV files loaded")
    return pd.concat(dfs, ignore_index=True)


def amdahl(p, f):
    """Amdahl's law: S(p) = 1 / (f + (1-f)/p)"""
    return 1.0 / (f + (1.0 - f) / p)


def plot_parallel(df: pd.DataFrame, out_dir: str) -> None:
    # Aggregate over repeats and seeds
    agg = (
        df.groupby("threads")
        .agg(
            speedup_mean=("speedup", lambda x: pd.to_numeric(x, errors="coerce").mean()),
            speedup_std=("speedup",  lambda x: pd.to_numeric(x, errors="coerce").std()),
            efficiency_mean=("efficiency", lambda x: pd.to_numeric(x, errors="coerce").mean()),
            wall_ms_mean=("wall_ms", "mean"),
            wall_ms_std=("wall_ms", "std"),
            peak_gb_mean=("peak_kb", lambda x: x.mean() / 1024**2),
            n=("wall_ms", "count"),
        )
        .reset_index()
        .sort_values("threads")
    )

    threads = agg["threads"].values
    speedup = agg["speedup_mean"].values
    efficiency = agg["efficiency_mean"].values

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # -- Panel A: Speedup --------------------------------------
    ax = axes[0]
    # Ideal linear speedup
    ax.plot(threads, threads, "k--", alpha=0.5, linewidth=1.5, label="Ideal (linear)")
    # Measured speedup
    ax.errorbar(threads, speedup,
                yerr=np.nan_to_num(agg["speedup_std"].values),
                color="#2166ac", marker="o", linewidth=2.5,
                markersize=8, capsize=4, label="G2G (measured)")

    # Amdahl fit
    try:
        valid = ~np.isnan(speedup)
        popt, _ = curve_fit(amdahl, threads[valid], speedup[valid], p0=[0.05],
                            bounds=(0, 1))
        f_serial = popt[0]
        p_fit = np.linspace(1, threads.max(), 200)
        ax.plot(p_fit, amdahl(p_fit, f_serial), "r-", alpha=0.6,
                label=f"Amdahl fit (f={f_serial:.3f})")
        max_speedup = 1.0 / f_serial
        ax.text(0.05, 0.7, f"Serial fraction: {f_serial*100:.1f}%\n"
                           f"Max speedup: {max_speedup:.1f}×",
                transform=ax.transAxes, fontsize=9,
                bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
    except Exception:
        pass

    ax.set_xlabel("Thread count")
    ax.set_ylabel("Speedup (T₁ / Tₚ)")
    ax.set_title("(A) Parallel Speedup")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x)}"))

    import matplotlib.ticker as ticker

    # -- Panel B: Parallel efficiency -------------------------
    ax = axes[1]
    ax.plot(threads, efficiency, color="#d73027", marker="s",
            linewidth=2.5, markersize=8, label="G2G efficiency")
    ax.axhline(y=1.0, color="k", linestyle="--", alpha=0.5, label="Perfect efficiency")
    ax.axhline(y=0.7, color="gray", linestyle=":", alpha=0.7, label="70% threshold")
    ax.set_xlabel("Thread count")
    ax.set_ylabel("Parallel efficiency (T₁ / (p·Tₚ))")
    ax.set_title("(B) Parallel Efficiency")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x)}"))
    ax.set_ylim(0, 1.1)

    # -- Panel C: Wall time ------------------------------------
    ax = axes[2]
    ax.errorbar(threads, agg["wall_ms_mean"].values / 1000,
                yerr=np.nan_to_num(agg["wall_ms_std"].values) / 1000,
                color="#4dac26", marker="^", linewidth=2.5,
                markersize=8, capsize=4, label="G2G wall time")
    ax.set_xlabel("Thread count")
    ax.set_ylabel("Wall time (seconds)")
    ax.set_title("(C) Wall Time vs. Thread Count")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x)}"))
    ax.set_yscale("log")

    plt.suptitle("G2G Parallel Scaling - Experiment 4", fontsize=13, fontweight="bold")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        path = os.path.join(out_dir, f"fig4_parallel_scaling.{ext}")
        dpi = 150 if ext == "pdf" else 200
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        print(f"[plot] Saved: {path}")
    plt.close(fig)

    # Print summary table
    print("\n[plot] Parallel scaling summary:")
    print(agg[["threads", "speedup_mean", "efficiency_mean", "wall_ms_mean"]].to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv", nargs="+", required=True)
    parser.add_argument("--out", default="results/phase5/figures/")
    args = parser.parse_args()

    csv_paths = []
    for pat in args.csv:
        expanded = glob.glob(pat)
        csv_paths.extend(expanded if expanded else [pat])

    os.makedirs(args.out, exist_ok=True)
    df = load_csv(csv_paths)
    print(f"[plot] {len(df)} rows  threads={sorted(df['threads'].unique())}")
    plot_parallel(df, args.out)


if __name__ == "__main__":
    main()
