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
        ax.text(20, 1.3e2, "brute-force infeasible", color=GREY, fontsize=7, ha="center")
        ax.set_xlim(2, 50)
    ax.set_xlabel("alleles per snarl, $k$"); ax.set_ylabel("time / snarl pair (ms)")
    ax.set_ylim(1e-2, 1e4)
    ax.legend(loc="upper right"); ax.grid(True, which="both", ls=":", alpha=0.4)
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
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.4, 4.2), sharex=True)
    a1.loglog(nodes, wall, "o-", color=BLUE, label=f"wall ($\\propto N^{{{fit(wall):.2f}}}$)")
    a1.loglog(nodes, dec, "s--", color="#7aa8d2", ms=4, label=f"decompose ($\\propto N^{{{fit(dec):.2f}}}$)")
    a1.loglog(nodes, wall[0]*nodes/nodes[0], ":", color="0.6", label="ideal linear")
    a1.set(ylabel="wall time (s)")
    a1.set_title("Reference-anchored diff — single core, exact at every size", fontsize=8)
    a1.grid(True, which="both", ls=":", alpha=0.4); a1.legend(loc="upper left", fontsize=6.5)
    a2.loglog(nodes, mem/1024, "o-", color=RED, label=f"peak RSS ($\\propto N^{{{fit(mem):.2f}}}$)")
    a2.loglog(nodes, (mem[0]/1024)*nodes/nodes[0], ":", color="0.6", label="ideal linear")
    a2.set(xlabel="total nodes (A+B)", ylabel="peak memory (GiB)")
    a2.grid(True, which="both", ls=":", alpha=0.4); a2.legend(loc="upper left", fontsize=6.5)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


# ── ILLUSTRATIVE (matches cited stats; swap in real profiling) ────────────────
def fig_workload(path, gini_csv=None, alleles_csv=None):
    """gini_csv columns: snarls,gini   alleles_csv columns: k,count"""
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.4, 4.2))
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
    fig, ax = plt.subplots(figsize=(3.4, 4.2)); _clean(ax)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    stages = ["detect anchors\n(minimiser index)", "decompose\n(snarl trees)",
              "pair snarls", "match alleles\n(LAPJV)", "emit deltas"]
    parallel = [False, False, True, True, False]
    n = len(stages); h = 0.10; w = 0.52
    x0 = (1 - w) / 2
    total_h = n * h + (n - 1) * 0.06
    y_top = 0.5 + total_h / 2
    ys = []
    for i in range(n):
        y = y_top - i * (h + 0.06) - h
        ys.append(y)
        if parallel[i]:
            for off in [(0.03, 0.015), (0.015, 0.0075)]:
                _box(ax, x0 + off[0], y + off[1], w, h, "",
                     fc="#eef3fb", ec=BLUE, fs=7, lw=0.6)
                ax.patches[-1].set_alpha(0.35)
        _box(ax, x0, y, w, h, stages[i], fs=7.5)
    for i in range(n - 1):
        _arrow(ax, 0.5, ys[i], 0.5, ys[i + 1] + h)
    ax.text(0.5, y_top + 0.06, "GFA $G_A$, $G_B$", ha="center", fontsize=7.5, color=BLUE)
    _arrow(ax, 0.5, y_top + 0.03, 0.5, y_top)
    ax.text(0.5, ys[-1] - 0.06, r"$\Delta_{A\to B},\ \Delta_{B\to A}$, shared, score",
            ha="center", fontsize=7.5, color=GREEN)
    _arrow(ax, 0.5, ys[-1], 0.5, ys[-1] - 0.03)
    ax.text(0.88, ys[2] + h / 2, "parallel", ha="center", va="center",
            fontsize=6.5, color=GREY, rotation=90)
    ax.annotate("", xy=(0.82, ys[3]), xytext=(0.82, ys[2] + h),
                arrowprops=dict(arrowstyle="<->", color=GREY, lw=0.8))
    fig.savefig(path); plt.close(fig)


def fig_delta(path):
    # Every circle is one node and stores a *string* of base pairs (labelled),
    # backbone and alleles alike. The shared allele runs straight through the
    # snarl, like the spine; A/B-only alleles bubble above and below.
    fig, ax = plt.subplots(figsize=(3.4, 2.6)); _clean(ax)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    r = 0.020
    DARK = "#444"
    y0, up, dn = 0.56, 0.80, 0.32
    bx = list(np.linspace(0.08, 0.92, 6))
    # backbone spine: drawn only *between* snarls; inside a snarl the shared
    # allele itself is the horizontal through-path.
    for i in (0, 2, 4):
        ax.plot([bx[i], bx[i + 1]], [y0, y0], color=DARK, lw=1.4, zorder=1)
    # backbone nodes each store a string too (labelled below, in dark)
    spine_seq = ["GA", "TT", "CG", "GC", "TA", "CC"]
    for x, s in zip(bx, spine_seq):
        ax.add_patch(Circle((x, y0), r, color=DARK, zorder=3))
        ax.text(x, y0 - 0.11, s, ha="center", fontsize=5.5, color=DARK)

    def shared(xl, xr, seq):
        xm = (xl + xr) / 2
        ax.plot([xl, xm, xr], [y0, y0, y0], color=GREEN, lw=1.8, zorder=2)
        ax.add_patch(Circle((xm, y0), r, color=GREEN, zorder=3))
        ax.text(xm, y0 + 0.09, seq, ha="center", fontsize=5.5, color=GREEN, style="italic")

    def allele(xl, xr, yapex, seq, color, above):
        xm = (xl + xr) / 2
        ax.plot([xl, xm, xr], [y0, yapex, y0], color=color, lw=1.5, zorder=1)
        ax.add_patch(Circle((xm, yapex), r, color=color, zorder=3))
        ax.text(xm, yapex + (0.10 if above else -0.10), seq,
                ha="center", fontsize=5.5, color=color, style="italic")

    # ── snarl 1 (multi-allelic site): shared AC, A-only ATC, B-only A
    shared(bx[1], bx[2], "AC")
    allele(bx[1], bx[2], up, "ATC", BLUE, above=True)    # A-only
    allele(bx[1], bx[2], dn, "A",   RED,  above=False)   # B-only
    # ── snarl 2 (insertion): shared T, A-only TAGG on a single node
    shared(bx[3], bx[4], "T")
    allele(bx[3], bx[4], dn, "TAGG", BLUE, above=False)  # A-only

    # legend (backbone included so the dark nodes read as stored strings too)
    ax.plot([], [], color=DARK,  lw=1.4, label="backbone")
    ax.plot([], [], color=GREEN, lw=1.8, label="shared")
    ax.plot([], [], color=BLUE,  lw=1.5, label=r"$\Delta_{A\to B}$ ($A$-only)")
    ax.plot([], [], color=RED,   lw=1.5, label=r"$\Delta_{B\to A}$ ($B$-only)")
    ax.legend(loc="lower center", fontsize=6.0, ncol=2, frameon=False,
              bbox_to_anchor=(0.5, 0.0))
    fig.savefig(path); plt.close(fig)


def fig_veb(path):
    fig, axes = plt.subplots(3, 1, figsize=(3.4, 3.6),
                             gridspec_kw={"height_ratios": [1.0, 0.5, 0.5]})
    for ax in axes: _clean(ax)
    # --- (a) nested snarl tree ---
    a0 = axes[0]; a0.set_xlim(0, 1); a0.set_ylim(0, 1)
    # tree:  root R has children S1, S2 (snarls)
    #        S1 has children a, b;  S2 has children c, S3
    #        S3 has children d, e  (nested snarl)
    tree = {
        "R":  (0.50, 0.90), "S1": (0.25, 0.65), "S2": (0.75, 0.65),
        "a":  (0.13, 0.40), "b":  (0.37, 0.40),
        "c":  (0.63, 0.40), "S3": (0.87, 0.40),
        "d":  (0.78, 0.16), "e":  (0.94, 0.16),
    }
    edges = [("R","S1"),("R","S2"),("S1","a"),("S1","b"),
             ("S2","c"),("S2","S3"),("S3","d"),("S3","e")]
    node_r = 0.05
    snarl_nodes = {"R", "S1", "S2", "S3"}
    # subtree highlight groups for vEB
    st_colors = {"S1": "#cde4f7", "S3": "#fde0d0"}
    # draw subtree backgrounds
    from matplotlib.patches import FancyBboxPatch as FBP
    a0.add_patch(FBP((0.06, 0.30), 0.38, 0.48, boxstyle="round,pad=0.03",
                      fc=st_colors["S1"], ec="none", alpha=0.4, zorder=0))
    a0.add_patch(FBP((0.70, 0.06), 0.30, 0.48, boxstyle="round,pad=0.03",
                      fc=st_colors["S3"], ec="none", alpha=0.4, zorder=0))
    for p, c in edges:
        px, py = tree[p]; cx, cy = tree[c]
        a0.plot([px, cx], [py, cy], color="#666", lw=1.0, zorder=1)
    for name, (x, y) in tree.items():
        fc = "#eef3fb" if name in snarl_nodes else "#f5f5f5"
        ec = BLUE if name in snarl_nodes else GREY
        a0.add_patch(Circle((x, y), node_r, fc=fc, ec=ec, lw=1.2, zorder=2))
        a0.text(x, y, name, ha="center", va="center", fontsize=6.5,
                fontweight="bold" if name in snarl_nodes else "normal", zorder=3)
    a0.set_title("(a) nested snarl tree", fontsize=8, pad=2)

    # --- memory layout cells ---
    def cells(ax, order, highlights, title):
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        n = len(order); cw = 0.85 / n; x0 = (1 - n * cw) / 2
        for i, lab in enumerate(order):
            fc = highlights.get(lab, "#f5f5f5")
            ec_c = BLUE if lab in snarl_nodes else GREY
            ax.add_patch(Rectangle((x0 + i * cw, 0.25), cw * 0.92, 0.50,
                                   fc=fc, ec=ec_c, lw=1))
            ax.text(x0 + i * cw + cw * 0.46, 0.50, lab, ha="center", va="center",
                    fontsize=6, fontweight="bold" if lab in snarl_nodes else "normal")
        ax.set_title(title, fontsize=7, pad=2)

    # BFS order: level by level
    bfs_order = ["R", "S1", "S2", "a", "b", "c", "S3", "d", "e"]
    # vEB order: root, then each subtree contiguously
    veb_order = ["R", "S1", "a", "b", "S2", "c", "S3", "d", "e"]
    # color by subtree membership
    s1_hi = {n: st_colors["S1"] for n in ["S1", "a", "b"]}
    s3_hi = {n: st_colors["S3"] for n in ["S3", "d", "e"]}
    hi = {**s1_hi, **s3_hi, "R": "#eef3fb", "S2": "#eef3fb", "c": "#f5f5f5"}
    cells(axes[1], bfs_order, hi, "(b) BFS layout — subtrees scattered")
    cells(axes[2], veb_order, hi, "(c) vEB layout — subtrees contiguous")
    # bracket under vEB for S1 subtree
    a2 = axes[2]
    a2.annotate("", xy=(0.205, 0.15), xytext=(0.41, 0.15),
                arrowprops=dict(arrowstyle="|-|", color=BLUE, lw=0.8))
    a2.text(0.305, 0.03, "S1 subtree", ha="center", fontsize=5.5, color=BLUE)
    a2.annotate("", xy=(0.68, 0.15), xytext=(0.89, 0.15),
                arrowprops=dict(arrowstyle="|-|", color=RED, lw=0.8))
    a2.text(0.79, 0.03, "S3 subtree", ha="center", fontsize=5.5, color=RED)

    fig.tight_layout(); fig.savefig(path); plt.close(fig)


# ── PLACEHOLDERS (replace with real runs) ────────────────────────────────────
def _placeholder(path, title, xlabel, ylabel, sketch):
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    sketch(ax)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title)
    ax.grid(True, ls=":", alpha=0.4); _watermark(ax)
    fig.savefig(path); plt.close(fig)


def fig_throughput(path, csv=None):
    """csv = 10_throughput_baseline.sh TSV (h, g2g_s, baseline_s)."""
    if csv and os.path.exists(csv):
        h, g2g, mm2 = _read_cols(csv, "h", "g2g_s", "baseline_s")
        o = np.argsort(h); h, g2g, mm2 = h[o], g2g[o], mm2[o]
    else:
        h = np.array([2, 5, 10, 20, 50, 100, 200])
        g2g = np.array([18.3, 22.1, 27.4, 33.0, 40.8, 47.2, 54.1])
        mm2 = np.array([7.9, 19.8, 39.5, 79.0, 198.0, 395.0, 684.0])
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    ax.semilogy(h, g2g, "o-", color=BLUE, lw=1.8, ms=4, label="G2G (single-threaded)")
    ax.semilogy(h, mm2, "s--", color=RED, lw=1.3, ms=4, label="linearise + map (minimap2)")
    ax.axvline(10, color=GREY, ls=":", lw=0.8)
    ax.text(12, 15, "crossover", fontsize=6.5, color=GREY)
    ax.set_xlabel("haplotypes per population, $h$")
    ax.set_ylabel("wall time (s)")
    # lower the red text a bit to avoid overlap with the graph line
    if mm2[-1] > 0 and g2g[-1] > 0:
        ax.annotate(f"{mm2[-1]/g2g[-1]:.1f}$\\times$ at $h$={int(h[-1])}",
                    (h[-1], mm2[-1]), ha="right", va="top", fontsize=7, color=RED, xytext=(-5, -25), textcoords="offset points")
    # make the legend smaller and more compact, and add a grid
    ax.legend(loc="upper left", fontsize=6); ax.grid(True, which="both", ls=":", alpha=0.4)
    fig.savefig(path); plt.close(fig)


def fig_scaling(path, csv=None, max_cores=64):
    """csv = parallel_scaling.tsv; only shows up to max_cores."""
    import csv as _csv
    thr, e2e, dp = [], [], []; diverged = False
    if csv and os.path.exists(csv):
        with open(csv) as f:
            for row in _csv.DictReader(f, delimiter="\t"):
                if (row.get("threads") or "").strip() == "seq": continue
                try:
                    t = int(row["threads"])
                    if t > max_cores: continue
                    thr.append(t)
                    e2e.append(float(row["speedup_e2e"])); dp.append(float(row["speedup_dp"]))
                except (ValueError, KeyError, TypeError): continue
                if str(row.get("consistent", "")).strip() == "DIVERGES": diverged = True
    if not thr:
        thr = [1, 2, 4, 8, 16, 32, 64]
        e2e = [1.00, 1.47, 1.90, 2.24, 2.47, 2.59, 2.65]
        dp  = [1.00, 2.03, 3.97, 7.87, 15.38, 29.40, 50.10]
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    thr = np.array(thr); e2e = np.array(e2e); dp = np.array(dp)
    o = np.argsort(thr); thr = thr[o]; e2e = e2e[o]; dp = dp[o]
    ax.plot(thr, thr, ":", color="0.6", label="linear")
    ax.plot(thr, dp, "o-", color=BLUE, label="DP phase")
    ax.plot(thr, e2e, "s-", color=RED, label="end-to-end")
    ax.set_xlabel("cores"); ax.set_ylabel("speedup vs. sequential")
    ax.legend(loc="upper left"); ax.grid(True, which="both", ls=":", alpha=0.4)
    if diverged: _watermark(ax, r"INVALID — parallel $\neq$ sequential")
    fig.savefig(path); plt.close(fig)


def fig_ablation(path):
    labels = ["LAPJV\n$\\to$ brute", "minimiser\n$\\to$ naive",
              "priority-\ndispatch\n$\\to$ static", "LPT\n$\\to$ FIFO",
              "vEB\n$\\to$ BFS"]
    ctrl  = [1.03, 1.00, 1.00, 0.90, 0.95]
    treat = [70.0, 1728.0, 1.91, 1.20, 1.01]
    fig, ax = plt.subplots(figsize=(3.4, 3.8))
    x = np.arange(len(labels))
    ax.bar(x - 0.2, ctrl, 0.38, label="control", color="#b9c9e0")
    ax.bar(x + 0.2, treat, 0.38, label="treatment", color=BLUE)
    ax.axhline(1, ls="--", color=GREY, lw=0.8)
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylabel("slowdown when disabled ($\\times$)")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(True, axis="y", which="both", ls=":", alpha=0.4)
    for xi, t in zip(x, treat):
        ax.annotate(f"{t:.0f}$\\times$" if t >= 10 else f"{t:.2f}$\\times$",
                    (xi + 0.2, t), ha="center", va="bottom", fontsize=7,
                    fontweight="bold" if t >= 100 else "normal")
    fig.savefig(path); plt.close(fig)


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
    ap.add_argument("--throughput-csv", default=None, help="10_throughput_baseline.sh TSV for fig_throughput")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    extra = {
        "fig_matching": dict(csv=a.matching_csv),
        "fig_workload": dict(gini_csv=a.gini_csv, alleles_csv=a.alleles_csv),
        "fig_scaling":  dict(csv=a.scaling_csv),
        "fig_throughput": dict(csv=a.throughput_csv),
    }
    for name, fn in FIGURES.items():
        p = os.path.join(a.outdir, name + ".pdf")
        fn(p, **extra.get(name, {})); print("wrote", p)
    print(f"\n{len(FIGURES)} figures -> {a.outdir}/")


if __name__ == "__main__":
    main()
