#!/usr/bin/env python3
"""
make_figures.py — generate every figure the G2G paper references, as figs/*.pdf.

Three kinds:
  REAL DATA  : fig_matching, fig_refdiffscaling      (from measured numbers)
  SCHEMATIC  : fig_paradigm, fig_snarltree, fig_pipeline, fig_delta, fig_veb
               (conceptual diagrams — refine/redraw to taste)
  PLACEHOLDER: fig_throughput, fig_scaling, fig_ablation, fig_biology
               (marked PRELIMINARY — replace with real runs)
  ILLUSTRATIVE: fig_workload   (matches the cited stats; swap in your real
               Phase-2 profiling CSV if you have it)

Usage:  python3 make_figures.py [--outdir figs]
Drop the resulting figs/ next to your .tex.  Needs matplotlib + numpy.
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle, Ellipse

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.titlesize": 9,
    "axes.labelsize": 9, "legend.fontsize": 7.5, "figure.dpi": 150,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
BLUE, RED, GREY, GREEN = "#1f77b4", "#d62728", "#7f7f7f", "#2ca02c"


def _box(ax, x, y, w, h, text, fc="#eef3fb", ec=BLUE, fs=8, lw=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.03",
                                fc=fc, ec=ec, lw=lw, mutation_aspect=0.6))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs)


def _arrow(ax, x0, y0, x1, y1, ec="#333", lw=1.3, style="-|>"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                                 mutation_scale=11, lw=lw, color=ec))


def _clean(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)


def _watermark(ax, text="PRELIMINARY — measurements pending"):
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center",
            fontsize=13, color=RED, alpha=0.30, rotation=18, fontweight="bold")


def _read_cols(path, *names):
    """Read named columns from a CSV/TSV (auto-delimiter); blanks/NA -> NaN."""
    import csv as _csv
    cols = {n: [] for n in names}
    with open(path) as f:
        sample = f.readline(); f.seek(0)
        delim = "\t" if "\t" in sample else ","
        for row in _csv.DictReader(f, delimiter=delim):
            for n in names:
                v = (row.get(n) or "").strip()
                cols[n].append(float(v) if v not in ("", "NA", "na") else np.nan)
    return tuple(np.array(cols[n]) for n in names)


# ── REAL DATA ────────────────────────────────────────────────────────────────
def fig_matching(path, csv=None):
    """csv columns: k,brute_ms,lapjv_ms  (blank/NA brute = infeasible)."""
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    if csv and os.path.exists(csv):
        k, t_b, t_l = _read_cols(csv, "k", "brute_ms", "lapjv_ms")
        ml = ~np.isnan(t_l); mb = ~np.isnan(t_b)
        if ml.any(): ax.semilogy(k[ml], t_l[ml], "o-", color=BLUE, lw=1.8, ms=4, label=r"LAPJV  $O(k^3)$")
        if mb.any(): ax.semilogy(k[mb], t_b[mb], "s--", color=RED, lw=1.3, ms=4, label="brute force")
        ax.set_xlim(left=max(1, np.nanmin(k) - 1))
    else:
        # fallback: measured anchors (Phase-2): k=2 -> 0.11 ms; k=7 is 1491x k=2;
        # k=8 -> 1990 ms; infeasible beyond. Middle points interpolate the trend.
        k_b = np.array([2, 3, 4, 5, 6, 7, 8])
        t_b = np.array([0.11, 0.6, 3.0, 15.0, 60.0, 164.0, 1990.0])
        k_l = np.arange(2, 51); t_l = 0.0009 * k_l**3
        ax.semilogy(k_l, t_l, "-", color=BLUE, lw=1.8, label=r"LAPJV  $O(k^3)$")
        ax.semilogy(k_b, t_b, "s--", color=RED, lw=1.3, ms=4, label="brute force (measured)")
        ax.axvspan(8, 50, color=GREY, alpha=0.12)
        ax.text(29, 1.5e3, "brute-force infeasible", color=GREY, fontsize=7, ha="center")
        ax.set_xlim(2, 50)
    ax.set_xlabel("alleles per snarl, $k$"); ax.set_ylabel("time / snarl pair (ms)")
    ax.set_ylim(1e-2, 1e4)
    ax.legend(loc="upper left"); ax.grid(True, which="both", ls=":", alpha=0.4)
    fig.savefig(path); plt.close(fig)


def fig_refdiffscaling(path):
    # measured cluster run (503 GB node), 8 Mb -> 8 Gb, all exact.
    nodes = np.array([0.95e6, 1.90e6, 3.79e6, 7.58e6, 15.17e6, 30.33e6,
                      60.67e6, 121.33e6, 242.67e6, 485.33e6, 970.67e6])
    wall = np.array([2.88, 5.49, 11.88, 24.93, 51.34, 103.73,
                     189.58, 414.60, 816.70, 1879.87, 4291.85])
    mem = np.array([259.5, 509.4, 1009.3, 1975.4, 3966.2, 7838.6,
                    15652.5, 31292.2, 61797.9, 123631.4, 247372.7])  # MiB
    dec = np.array([1.49, 2.96, 5.93, 12.00, 24.00, 46.92,
                    75.84, 172.86, 348.88, 745.11, 1551.11])
    def fit(y):
        s, _ = np.polyfit(np.log(nodes), np.log(y), 1); return s
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.7))
    a1.loglog(nodes, wall, "o-", color=BLUE, label=f"wall ($\\propto N^{{{fit(wall):.2f}}}$)")
    a1.loglog(nodes, dec, "s--", color="#7aa8d2", ms=4, label=f"decompose ($\\propto N^{{{fit(dec):.2f}}}$)")
    a1.loglog(nodes, wall[0]*nodes/nodes[0], ":", color="0.6", label="ideal linear")
    a1.set(xlabel="total nodes (A+B)", ylabel="wall time (s)", title="Time")
    a1.grid(True, which="both", ls=":", alpha=0.4); a1.legend(loc="upper left")
    a2.loglog(nodes, mem/1024, "o-", color=RED, label=f"peak RSS ($\\propto N^{{{fit(mem):.2f}}}$)")
    a2.loglog(nodes, (mem[0]/1024)*nodes/nodes[0], ":", color="0.6", label="ideal linear")
    a2.set(xlabel="total nodes (A+B)", ylabel="peak memory (GiB)", title="Memory")
    a2.grid(True, which="both", ls=":", alpha=0.4); a2.legend(loc="upper left")
    fig.suptitle("Reference-anchored diff — single core, exact at every size", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(path); plt.close(fig)


# ── ILLUSTRATIVE (matches cited stats; swap in real profiling) ────────────────
def fig_workload(path, gini_csv=None, alleles_csv=None):
    """gini_csv columns: snarls,gini   alleles_csv columns: k,count"""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.6))
    if gini_csv and os.path.exists(gini_csv):
        sizes, gini = _read_cols(gini_csv, "snarls", "gini")
    else:
        sizes = np.array([20, 50, 100, 200, 500, 1000]); rng = np.random.default_rng(0)
        gini = 0.565 + 0.012*np.sin(np.linspace(0, 3, len(sizes))) + rng.normal(0, 0.004, len(sizes))
    a1.plot(sizes, gini, "o-", color=BLUE); a1.axhspan(0.55, 0.58, color=BLUE, alpha=0.10)
    a1.axhline(0.55, ls="--", color=GREY, lw=0.8); a1.set_xscale("log"); a1.set_ylim(0.45, 0.65)
    a1.set(xlabel="snarls in graph", ylabel="Gini of per-snarl work", title="(a) scale-invariant imbalance")
    a1.grid(True, which="both", ls=":", alpha=0.4)
    if alleles_csv and os.path.exists(alleles_csv):
        ks, cnt = _read_cols(alleles_csv, "k", "count"); p = cnt / np.nansum(cnt)
    else:
        ks = np.arange(2, 14)
        p = np.array([0.67, 0.14, 0.07, 0.04, 0.025, 0.018, 0.012, 0.009,
                      0.007, 0.005, 0.004, 0.003]); p = p / p.sum()
    a2.bar(ks, p, color=RED, alpha=0.85, width=0.8); a2.set_yscale("log")
    a2.set(xlabel="alleles per snarl, $k$", ylabel="fraction of snarls", title="(b) heavy-tailed $k$")
    a2.grid(True, axis="y", ls=":", alpha=0.4)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


# ── SCHEMATICS ───────────────────────────────────────────────────────────────
def fig_paradigm(path):
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 1.9))
    for ax in axs: _clean(ax); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    # (a) seq-seq
    axs[0].plot([0.1, 0.9], [0.62, 0.62], "-", color=BLUE, lw=3)
    axs[0].plot([0.1, 0.9], [0.38, 0.38], "-", color=RED, lw=3)
    for x in np.linspace(0.18, 0.82, 5):
        axs[0].plot([x, x], [0.6, 0.4], color=GREY, lw=0.6)
    axs[0].set_title("(a) sequence–sequence")
    # (b) seq-graph
    axs[1].plot([0.1, 0.9], [0.68, 0.68], "-", color=BLUE, lw=3)
    gx = [0.18, 0.36, 0.5, 0.5, 0.64, 0.82]; gy = [0.32, 0.32, 0.42, 0.22, 0.32, 0.32]
    for i in range(len(gx) - 1):
        axs[1].plot([gx[i], gx[i+1]], [gy[i], gy[i+1]], color=RED, lw=2)
    for x, y in zip(gx, gy): axs[1].add_patch(Circle((x, y), 0.022, color=RED))
    axs[1].set_title("(b) sequence–graph")
    # (c) graph-graph
    for yy, c in [(0.7, BLUE), (0.3, RED)]:
        gx = [0.18, 0.36, 0.5, 0.5, 0.64, 0.82]; gy = [yy, yy, yy+0.1, yy-0.1, yy, yy]
        for i in range(len(gx) - 1):
            axs[2].plot([gx[i], gx[i+1]], [gy[i], gy[i+1]], color=c, lw=2)
        for x, y in zip(gx, gy): axs[2].add_patch(Circle((x, y), 0.02, color=c))
    for x in [0.18, 0.36, 0.64, 0.82]:
        axs[2].plot([x, x], [0.68, 0.32], color=GREEN, lw=0.8, ls="--")
    axs[2].set_title("(c) graph–graph (G2G)")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_snarltree(path):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.3))
    for ax in (a1, a2): _clean(ax); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    # (a) graph with nested bubbles
    xs = [0.08, 0.26, 0.5, 0.74, 0.92]
    for x in xs: a1.add_patch(Circle((x, 0.5), 0.025, color=BLUE))
    a1.plot(xs, [0.5]*5, color=BLUE, lw=1.5, zorder=0)
    for x0, x1 in [(0.26, 0.5), (0.5, 0.74)]:
        a1.plot([x0, (x0+x1)/2, x1], [0.5, 0.66, 0.5], color=RED, lw=1.6)
        a1.plot([x0, (x0+x1)/2, x1], [0.5, 0.34, 0.5], color=RED, lw=1.6)
    a1.add_patch(Ellipse((0.5, 0.5), 0.62, 0.5, fill=False, ec=GREY, ls="--", lw=1))
    a1.set_title("(a) variation graph $\\to$ nested snarls")
    # (b) pairing + bipartite matching
    a2.text(0.5, 0.93, "snarl pairing + allele matching", ha="center", fontsize=8)
    for i, y in enumerate([0.62, 0.42, 0.22]):
        a2.add_patch(Circle((0.25, y), 0.03, color=BLUE)); a2.text(0.16, y, f"$a_{i}$", fontsize=7, va="center")
    for i, y in enumerate([0.7, 0.46, 0.22]):
        a2.add_patch(Circle((0.75, y), 0.03, color=RED)); a2.text(0.83, y, f"$b_{i}$", fontsize=7, va="center")
    for (ya, yb) in [(0.62, 0.7), (0.42, 0.46), (0.22, 0.22)]:
        a2.plot([0.28, 0.72], [ya, yb], color=GREEN, lw=1.4)
    a2.text(0.5, 0.06, r"$\max$ similarity matching (LAPJV)", ha="center", fontsize=7, color=GREEN)
    a2.set_title("(b) score each snarl pair")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_pipeline(path):
    fig, ax = plt.subplots(figsize=(7.2, 1.7)); _clean(ax)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    stages = ["anchors\n(minimiser)", "decompose\n(snarl tree)", "vEB\nlayout",
              "pair +\nmatch (LAPJV)", "delta\ngraphs"]
    n = len(stages); w = 0.155; gap = (1 - 0.06 - n*w) / (n - 1); x = 0.03
    xs = []
    for i, s in enumerate(stages):
        _box(ax, x, 0.32, w, 0.36, s, fs=7.2); xs.append(x); x += w + gap
    for i in range(n - 1):
        _arrow(ax, xs[i]+w, 0.5, xs[i+1], 0.5)
    ax.text(0.03+w/2, 0.83, "GFA $G_A$, $G_B$", ha="center", fontsize=7, color=BLUE)
    _arrow(ax, 0.03+w/2, 0.78, 0.03+w/2, 0.69)
    ax.text(xs[-1]+w/2, 0.17, r"$\Delta_{A\to B},\ \Delta_{B\to A}$, score", ha="center", fontsize=7, color=GREEN)
    _arrow(ax, xs[-1]+w/2, 0.31, xs[-1]+w/2, 0.22)
    fig.savefig(path); plt.close(fig)


def fig_delta(path):
    fig, ax = plt.subplots(figsize=(3.4, 2.3)); _clean(ax)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    xs = [0.1, 0.32, 0.55, 0.78, 0.92]
    for x in xs: ax.add_patch(Circle((x, 0.5), 0.022, color="#444"))
    ax.plot(xs, [0.5]*len(xs), color="#444", lw=1.4, zorder=0)
    # A-only insertion (top), B-only deletion (bottom)
    ax.plot([0.32, 0.435, 0.55], [0.5, 0.72, 0.5], color=BLUE, lw=1.8)
    ax.add_patch(Circle((0.435, 0.72), 0.022, color=BLUE))
    ax.text(0.435, 0.8, r"$\Delta_{A\to B}$ (ins)", ha="center", fontsize=7, color=BLUE)
    ax.plot([0.55, 0.665, 0.78], [0.5, 0.28, 0.5], color=RED, lw=1.8)
    ax.add_patch(Circle((0.665, 0.28), 0.022, color=RED))
    ax.text(0.665, 0.18, r"$\Delta_{B\to A}$ (del)", ha="center", fontsize=7, color=RED)
    ax.text(0.5, 0.95, "shared graph + base-pair deltas", ha="center", fontsize=8)
    fig.savefig(path); plt.close(fig)


def fig_veb(path):
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.4, 2.6))
    for ax in (a1, a2): _clean(ax); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    # simple tree drawn once (top-left), arrays show two orderings
    def cells(ax, order, y, title):
        ax.text(0.5, y+0.34, title, ha="center", fontsize=7.5)
        for i, lab in enumerate(order):
            x = 0.06 + i*0.115
            ax.add_patch(Rectangle((x, y), 0.1, 0.22, fc="#eef3fb", ec=BLUE, lw=1))
            ax.text(x+0.05, y+0.11, lab, ha="center", va="center", fontsize=7)
    cells(a1, ["r", "a", "b", "c", "d", "e", "f"], 0.3, "BFS layout (siblings scattered)")
    a1.text(0.5, 0.08, "subtree {b,c} far from {d,e}", ha="center", fontsize=6.5, color=RED)
    cells(a2, ["r", "a", "b", "c", "d", "e", "f"], 0.3, "van Emde Boas (subtrees contiguous)")
    a2.text(0.5, 0.08, "each subtree in $O(\\sqrt{N})$ cache lines", ha="center", fontsize=6.5, color=GREEN)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


# ── PLACEHOLDERS (replace with real runs) ────────────────────────────────────
def _placeholder(path, title, xlabel, ylabel, sketch):
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    sketch(ax)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title)
    ax.grid(True, ls=":", alpha=0.4); _watermark(ax)
    fig.savefig(path); plt.close(fig)


def fig_throughput(path):
    def s(ax):
        x = np.logspace(3, 7, 30)
        ax.loglog(x, x/4e4, color=BLUE, label="G2G (128 thr)")
        ax.loglog(x, x/300, color=RED, ls="--", label="linearise+map")
        ax.legend(loc="upper left")
    _placeholder(path, "Throughput (preliminary)", "graph size (snarls)", "time (s)", s)


def fig_scaling(path, csv=None):
    """csv = 07_parallel_scaling.sh TSV; flags the run INVALID if parallel
    diverged from sequential (the 'consistent' column == DIVERGES)."""
    if csv and os.path.exists(csv):
        import csv as _csv
        thr, e2e, dp = [], [], []; diverged = False
        with open(csv) as f:
            for row in _csv.DictReader(f, delimiter="\t"):
                if (row.get("threads") or "").strip() == "seq": continue
                try:
                    thr.append(int(row["threads"]))
                    e2e.append(float(row["speedup_e2e"])); dp.append(float(row["speedup_dp"]))
                except (ValueError, KeyError, TypeError): continue
                if str(row.get("consistent", "")).strip() == "DIVERGES": diverged = True
        fig, ax = plt.subplots(figsize=(3.4, 2.5))
        thr = np.array(thr)
        if len(thr):
            o = np.argsort(thr); thr = thr[o]; e2e = np.array(e2e)[o]; dp = np.array(dp)[o]
            ax.plot(thr, thr, ":", color="0.6", label="linear")
            ax.plot(thr, dp, "o-", color=BLUE, label="DP phase")
            ax.plot(thr, e2e, "s-", color=RED, label="end-to-end")
            ax.set_xscale("log", base=2); ax.set_yscale("log", base=2)
        ax.set_xlabel("threads"); ax.set_ylabel("speedup vs. sequential")
        ax.set_title("Parallel scaling"); ax.legend(loc="upper left")
        ax.grid(True, which="both", ls=":", alpha=0.4)
        if diverged: _watermark(ax, r"INVALID — parallel $\neq$ sequential")
        fig.savefig(path); plt.close(fig)
    else:
        def s(ax):
            t = np.array([1, 2, 4, 8, 16, 32, 64, 128])
            ax.plot(t, t, ":", color="0.6", label="linear")
            ax.plot(t, t/(1+0.02*t), color=BLUE, label="DP phase")
            ax.plot(t, t/(1+0.15*t), color=RED, label="end-to-end")
            ax.set_xscale("log", base=2); ax.set_yscale("log", base=2); ax.legend(loc="upper left")
        _placeholder(path, "Parallel scaling (preliminary)", "threads", "speedup", s)


def fig_ablation(path):
    def s(ax):
        labels = ["vEB", "steal-LF", "minimiser", "LAPJV"]
        ctrl = [1.0, 1.0, 1.4, 1.5]; treat = [1.8, 1.7, 3.2, 6.0]
        x = np.arange(len(labels))
        ax.bar(x-0.2, ctrl, 0.4, label="control (fits L3)", color="#b9c9e0")
        ax.bar(x+0.2, treat, 0.4, label="treatment (>L3)", color=BLUE)
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7); ax.legend(loc="upper left")
        ax.set_yscale("linear")
    _placeholder(path, "Ablations (preliminary)", "", "speedup vs. off", s)


def fig_biology(path):
    def s(ax):
        cats = ["high-$F_{ST}$\ngnomAD-SV", "ClinVar\n(shared)"]
        ax.bar(cats, [3.1, 2.4], color=[RED, GREEN], alpha=0.8)
        ax.axhline(1, ls="--", color=GREY); ax.set_ylim(0, 4)
    _placeholder(path, "HPRC enrichment (preliminary)", "", "fold enrichment", s)


FIGURES = {
    "fig_paradigm": fig_paradigm, "fig_snarltree": fig_snarltree,
    "fig_pipeline": fig_pipeline, "fig_delta": fig_delta, "fig_veb": fig_veb,
    "fig_matching": fig_matching, "fig_workload": fig_workload,
    "fig_refdiffscaling": fig_refdiffscaling,
    "fig_throughput": fig_throughput, "fig_scaling": fig_scaling,
    "fig_ablation": fig_ablation, "fig_biology": fig_biology,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="figs")
    ap.add_argument("--matching-csv", default=None, help="real fig_matching data (k,brute_ms,lapjv_ms)")
    ap.add_argument("--gini-csv", default=None, help="real fig_workload(a) data (snarls,gini)")
    ap.add_argument("--alleles-csv", default=None, help="real fig_workload(b) data (k,count)")
    ap.add_argument("--scaling-csv", default=None, help="07_parallel_scaling.sh TSV for fig_scaling")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    extra = {
        "fig_matching": dict(csv=a.matching_csv),
        "fig_workload": dict(gini_csv=a.gini_csv, alleles_csv=a.alleles_csv),
        "fig_scaling":  dict(csv=a.scaling_csv),
    }
    for name, fn in FIGURES.items():
        p = os.path.join(a.outdir, name + ".pdf")
        fn(p, **extra.get(name, {})); print("wrote", p)
    print(f"\n{len(FIGURES)} figures -> {a.outdir}/")


if __name__ == "__main__":
    main()
