#!/usr/bin/env python3
"""
plot_scaling.py - turn scaling-sweep TSVs into a log-log scaling figure.

Reads one or more result TSVs written by 05_scaling_demo.sh and/or
06_scaling_limit.sh and plots wall time and peak memory against graph size
(total nodes) on log-log axes, each with a fitted power-law exponent and an
ideal-linear guide line - so the near-linear scaling is obvious at a glance in
a writeup. A FAIL in any row's `correct` column is flagged in the title.

Usage:
    python3 demo/plot_scaling.py [TSV ...] [--out FILE] [--x nodes|bp]

With no TSV arguments it auto-discovers the demo's own outputs under demo/work/.
Both TSV layouts (05's and 06's) are accepted; overlapping sizes are de-duped,
so you can pass both to get a single curve spanning the whole range.
"""
import argparse
import csv
import os
import sys

DEFAULT_TSVS = [
    "demo/work/scaling/scale_results.tsv",
    "demo/work/scaling_limit/scale_limit_results.tsv",
]


def read_rows(paths):
    """Merge result rows from the given TSVs, de-duped by reference bp, sorted."""
    by_bp = {}
    for p in paths:
        if not os.path.exists(p):
            continue
        with open(p) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                try:
                    by_bp[int(r["ref_bp"])] = {
                        "bp":        int(r["ref_bp"]),
                        "nodes":     int(r["nodes_A"]) + int(r["nodes_B"]),
                        "wall":      float(r["wall_s"]),
                        "decompose": float(r.get("decompose_s") or "nan"),
                        "mem":       float(r["peak_MB"]),
                        "ok":        r.get("correct", "PASS"),
                    }
                except (KeyError, ValueError):
                    continue
    return [by_bp[k] for k in sorted(by_bp)]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", nargs="*", help="result TSV(s); default: demo/work/*")
    ap.add_argument("--out", default="demo/sample_outputs/scaling_curve.png")
    ap.add_argument("--x", choices=["nodes", "bp"], default="nodes",
                    help="x-axis: total nodes (default) or reference bp")
    a = ap.parse_args()

    try:
        import matplotlib
        matplotlib.use("Agg")            # headless: works over SSH / on a cluster
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as e:
        sys.exit(f"ERROR: matplotlib + numpy required ({e}).\n"
                 f"  pip install matplotlib   (or: conda install -c conda-forge matplotlib)")

    rows = read_rows(a.tsv or DEFAULT_TSVS)
    if len(rows) < 2:
        sys.exit("ERROR: need >=2 data points. Run 05/06 first, or pass TSV path(s).")

    x   = np.array([r[a.x]  for r in rows], float)
    wall = np.array([r["wall"] for r in rows], float)
    mem  = np.array([r["mem"]  for r in rows], float)
    dec  = np.array([r["decompose"] for r in rows], float)
    all_ok = all(r["ok"] == "PASS" for r in rows)
    xlabel = "total nodes (A + B)" if a.x == "nodes" else "reference length (bp)"

    def fit(y):
        slope, intercept = np.polyfit(np.log(x), np.log(y), 1)
        return slope, np.exp(intercept)

    def guide(ax, y0):                   # ideal linear (slope 1) through 1st point
        ax.loglog(x, y0 * (x / x[0]), ":", color="0.6", lw=1.3,
                  label="ideal linear (N¹)")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # -- Time ----------------------------------------------------------------
    sw, _ = fit(wall)
    ax1.loglog(x, wall, "o-", color="#1f77b4", label=f"wall  (∝ N^{sw:.2f})")
    if np.isfinite(dec).all() and (dec > 0).all():
        sd, _ = fit(dec)
        ax1.loglog(x, dec, "s--", color="#7aa8d2", ms=4,
                   label=f"decompose  (∝ N^{sd:.2f})")
    guide(ax1, wall[0])
    ax1.set(xlabel=xlabel, ylabel="wall time (s)", title="Time")
    ax1.grid(True, which="both", ls=":", alpha=0.4)
    ax1.legend(fontsize=8, loc="upper left")

    # -- Memory --------------------------------------------------------------
    sm, _ = fit(mem)
    ax2.loglog(x, mem, "o-", color="#d62728", label=f"peak RSS  (∝ N^{sm:.2f})")
    guide(ax2, mem[0])
    ax2.set(xlabel=xlabel, ylabel="peak memory (MB)", title="Memory")
    ax2.grid(True, which="both", ls=":", alpha=0.4)
    ax2.legend(fontsize=8, loc="upper left")

    big = rows[-1]
    status = "exact at every size" if all_ok else "WARNING - a point FAILED correctness"
    fig.suptitle(f"G2G  --ref  scaling  ·  {len(rows)} points  ·  "
                 f"up to {big['nodes']/1e6:.1f}M nodes / {big['bp']/1e6:.0f} Mb  ·  {status}",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=150)

    print(f"wrote {a.out}")
    print(f"  {len(rows)} points · wall ∝ N^{sw:.2f} · mem ∝ N^{sm:.2f}")
    print(f"  largest: {big['nodes']:,} nodes / {big['wall']:.1f} s / {big['mem']:.0f} MB"
          + ("" if all_ok else "  [!] correctness FAIL present"))


if __name__ == "__main__":
    main()
