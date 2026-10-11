#!/usr/bin/env python3
"""
make_figures.py — regenerate every figure of the PanG2G paper (latex/figs/*.pdf).

Every data figure is drawn from a committed measurement file, so a plain

    python3 make_figures.py            # -> latex/figs/

reproduces the paper's figures exactly; pass --outdir to write elsewhere, or
point an individual --*-csv flag at a fresh run (the harnesses under demo/ and
experiments/ write the same formats).

  figure               data (default path)                            produced by
  fig_matching         results/matching_cost_vs_k_cpp.csv              build/release/bench_matching
  fig_workload         results/work_distribution.csv (Gini) +          python/profile_workload.py
                       results/allele_count_hist.csv  (k histogram)
  fig_refdiffscaling   results/refdiff_scaling.tsv                     demo/06_scaling_limit.sh
  fig_throughput       demo/sample_outputs/throughput_chr22.tsv        demo/10_throughput_baseline.sh
  fig_scaling          results/parallel_scaling.tsv                    demo/07_parallel_scaling.sh
  fig_ablation         results/accl/ablation_results.tsv               demo/08_ablation.sh (EPYC 7763 node; the
                       (results/ablation_ratios.tsv = submitted run)    submitted version's ratios kept alongside)
  fig_biology          demo/sample_outputs/biology_results.tsv         demo/09_hprc_biology.sh
  fig_fst_sweep        demo/sample_outputs/fst_sweep_chr22.tsv         demo/hprc_enrich.py --fst-sweep
  fig_pairing_error    results/accl/exp1b_pairing_error.csv            experiments/exp1_correctness/pairing_error_sweep.py
                       (results/exp1b_pairing_error.csv = laptop run)   (EPYC 7763 node)
  fig_pipeline, fig_delta, fig_veb, fig_paradigm, fig_snarltree        schematics (no data)

Needs matplotlib + numpy.
"""
import argparse
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle, Ellipse

ROOT = os.path.dirname(os.path.abspath(__file__))
R = lambda *p: os.path.join(ROOT, *p)   # noqa: E731

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.titlesize": 9,
    "axes.labelsize": 9, "legend.fontsize": 7.5, "figure.dpi": 150,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
BLUE, RED, GREY, GREEN = "#1f77b4", "#d62728", "#7f7f7f", "#2ca02c"


# ── helpers ───────────────────────────────────────────────────────────────────
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


def _rows(path):
    """Rows of a CSV/TSV (auto-delimiter) as dicts of strings."""
    with open(path, newline="") as f:
        sample = f.readline(); f.seek(0)
        delim = "\t" if "\t" in sample else ","
        return list(csv.DictReader(f, delimiter=delim))


def _col(rows, name):
    out = []
    for r in rows:
        v = (r.get(name) or "").strip()
        out.append(float(v) if v not in ("", "NA", "na") else np.nan)
    return np.array(out)


def _need(path):
    if not os.path.exists(path):
        raise SystemExit(f"missing data file: {path}")
    return path


# ── data figures ──────────────────────────────────────────────────────────────
def fig_matching(path, csv_path):
    """k, brute_ms, lapjv_ms — measured by tools/bench_matching.cpp (C++)."""
    rows = _rows(_need(csv_path))
    k, t_b, t_l = _col(rows, "k"), _col(rows, "brute_ms"), _col(rows, "lapjv_ms")
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    ml, mb = ~np.isnan(t_l), ~np.isnan(t_b)
    ax.semilogy(k[ml], t_l[ml], "o-", color=BLUE, lw=1.6, ms=3, label=r"LAPJV  $O(k^3)$")
    ax.semilogy(k[mb], t_b[mb], "s--", color=RED, lw=1.3, ms=4, label="brute force  $O(k!)$")
    kmax_b = int(np.nanmax(k[mb]))
    ax.axvspan(kmax_b, np.nanmax(k), color=GREY, alpha=0.12)
    ax.text((kmax_b + np.nanmax(k)) / 2, np.nanmax(t_b) * 0.5,
            f"brute force infeasible\n($k>{kmax_b}$)", color=GREY, fontsize=7, ha="center")
    ax.set_xlabel("alleles per snarl, $k$"); ax.set_ylabel("time / snarl pair (ms)")
    ax.set_xlim(2, np.nanmax(k))
    ax.legend(loc="upper left"); ax.grid(True, which="both", ls=":", alpha=0.4)
    fig.savefig(path); plt.close(fig)


def fig_workload(path, gini_csv, alleles_csv):
    """(a) Gini of per-snarl work vs. graph size; (b) allele-count histogram."""
    g = _rows(_need(gini_csv))
    sizes, gini = _col(g, "total_snarls"), _col(g, "gini_work")
    h = _rows(_need(alleles_csv))
    ks, cnt = _col(h, "k"), _col(h, "count"); p = cnt / np.nansum(cnt)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.4, 3.3))
    a1.plot(sizes, gini, "o-", color=BLUE); a1.axhspan(0.55, 0.58, color=BLUE, alpha=0.10)
    a1.axhline(0.55, ls="--", color=GREY, lw=0.8); a1.set_xscale("log"); a1.set_ylim(0.45, 0.65)
    a1.set(xlabel="snarls in graph pair", ylabel="Gini of per-snarl work",
           title="(a) scale-invariant imbalance")
    a1.grid(True, which="both", ls=":", alpha=0.4)
    a2.bar(ks, p, color=RED, alpha=0.85, width=0.8); a2.set_yscale("log")
    a2.set_xticks(ks)
    a2.set(xlabel="alleles per snarl, $k$", ylabel="fraction of snarls", title="(b) heavy-tailed $k$")
    for x, y in zip(ks, p):
        a2.annotate(f"{100*y:.0f}%", (x, y), ha="center", va="bottom", fontsize=6.5)
    a2.grid(True, axis="y", ls=":", alpha=0.4)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_refdiffscaling(path, tsv):
    """Single panel: wall time and decomposition time (left axis) and peak
    memory (right axis) vs. total nodes, with fitted exponents."""
    rows = _rows(_need(tsv))
    nodes = _col(rows, "nodes_A") + _col(rows, "nodes_B")
    wall, mem, dec = _col(rows, "wall_s"), _col(rows, "peak_MB"), _col(rows, "decompose_s")
    def fit(y):
        s, _ = np.polyfit(np.log(nodes), np.log(y), 1); return s
    fig, a1 = plt.subplots(figsize=(3.4, 2.4))
    l1 = a1.loglog(nodes, wall, "o-", color=BLUE, ms=4, label=f"wall ($\\propto N^{{{fit(wall):.2f}}}$)")
    l2 = a1.loglog(nodes, dec, "s--", color="#7aa8d2", ms=3.5, label=f"decompose ($\\propto N^{{{fit(dec):.2f}}}$)")
    l0 = a1.loglog(nodes, wall[0] * nodes / nodes[0], ":", color="0.6", label="ideal linear")
    a1.set(xlabel="total nodes (A+B)", ylabel="time (s)")
    a1.grid(True, which="both", ls=":", alpha=0.4)
    a2 = a1.twinx()
    l3 = a2.loglog(nodes, mem / 1024, "^-", color=RED, ms=4, label=f"peak RSS ($\\propto N^{{{fit(mem):.2f}}}$)")
    a2.set_ylabel("peak memory (GiB)", color=RED); a2.tick_params(axis="y", colors=RED)
    lines = l1 + l2 + l3 + l0
    a1.legend(lines, [l.get_label() for l in lines], loc="upper left", fontsize=6.2)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_workload_matching(path, gini_csv, matching_csv):
    """(a) Gini of per-snarl work vs. graph size (Finding 2);
    (b) matching cost vs. allele count k, brute force vs. LAPJV (Finding 1)."""
    g = _rows(_need(gini_csv))
    sizes, gini = _col(g, "total_snarls"), _col(g, "gini_work")
    m = _rows(_need(matching_csv))
    k, t_b, t_l = _col(m, "k"), _col(m, "brute_ms"), _col(m, "lapjv_ms")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(3.5, 1.9))
    a1.plot(sizes, gini, "o-", color=BLUE, ms=4); a1.axhspan(0.55, 0.58, color=BLUE, alpha=0.10)
    a1.set_xscale("log"); a1.set_ylim(0.45, 0.65)
    a1.set_xlabel("snarls in graph pair", fontsize=7.5); a1.set_ylabel("Gini of per-snarl work", fontsize=7.5)
    a1.set_title("(a) imbalance is scale-invariant", fontsize=7.5)
    a1.tick_params(labelsize=6.5); a1.grid(True, which="both", ls=":", alpha=0.4)
    ml, mb = ~np.isnan(t_l), ~np.isnan(t_b)
    a2.semilogy(k[ml], t_l[ml], "-", color=BLUE, lw=1.6, label=r"LAPJV $O(k^3)$")
    a2.semilogy(k[mb], t_b[mb], "s--", color=RED, lw=1.2, ms=3.5, label=r"brute force $O(k!)$")
    kmax_b = int(np.nanmax(k[mb]))
    a2.axvspan(kmax_b, np.nanmax(k), color=GREY, alpha=0.12)
    a2.set_xlim(2, np.nanmax(k))
    a2.set_xlabel("alleles per snarl, $k$", fontsize=7.5); a2.set_ylabel("ms per snarl pair", fontsize=7.5)
    a2.set_title("(b) matching cost", fontsize=7.5)
    a2.tick_params(labelsize=6.5); a2.legend(fontsize=5.8, loc="upper left")
    a2.grid(True, which="both", ls=":", alpha=0.4)
    fig.tight_layout(w_pad=0.6); fig.savefig(path); plt.close(fig)


def fig_perf(path, scaling_tsv, throughput_tsv, max_cores=64, submitted_tsv=None):
    """(a) parallel speedup vs. cores; (b) PanG2G vs. linearize-and-map vs. h."""
    rows = [r for r in _rows(_need(scaling_tsv)) if r["threads"].strip() != "seq"]
    thr = _col(rows, "threads"); e2e = _col(rows, "speedup_e2e"); dp = _col(rows, "speedup_dp")
    keep = thr <= max_cores; thr, e2e, dp = thr[keep], e2e[keep], dp[keep]
    t = _rows(_need(throughput_tsv))
    h, g2g, mm2 = _col(t, "h"), _col(t, "g2g_s"), _col(t, "baseline_s")
    o = np.argsort(h); h, g2g, mm2 = h[o], g2g[o], mm2[o]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(3.5, 1.9))
    a1.plot(thr, thr, ":", color="0.6", label="linear")
    if submitted_tsv and os.path.exists(submitted_tsv):
        sr = [r for r in _rows(submitted_tsv) if r["threads"].strip() != "seq"]
        a1.plot(_col(sr, "threads"), _col(sr, "speedup_dp"), "^--", color="0.45", ms=3, lw=0.9,
                label="DP, submitted run")
    a1.plot(thr, dp, "o-", color=BLUE, ms=3.5, label="DP phase")
    a1.plot(thr, e2e, "s-", color=RED, ms=3.5, label="end-to-end")
    a1.set_xlabel("cores", fontsize=7.5); a1.set_ylabel("speedup vs. 1 thread", fontsize=7.5)
    a1.set_title("(a) parallel scaling", fontsize=7.5)
    a1.tick_params(labelsize=6.5); a1.legend(fontsize=5.8, loc="upper left")
    a1.grid(True, which="both", ls=":", alpha=0.4)
    a2.semilogy(h, g2g, "o-", color=BLUE, lw=1.5, ms=3.5, label="PanG2G (1 thread)")
    a2.semilogy(h, mm2, "s--", color=RED, lw=1.2, ms=3.5, label="linearize + map")
    a2.set_xlabel("haplotypes $h$ per population", fontsize=6.5); a2.set_ylabel("wall time (s)", fontsize=7.5)
    a2.set_title("(b) vs. linearize-and-map", fontsize=7.5)
    a2.annotate(f"{mm2[-1]/g2g[-1]:.1f}$\\times$", (h[-1], mm2[-1]), ha="right", va="top",
                fontsize=6.5, color=RED, xytext=(-2, -4), textcoords="offset points")
    a2.tick_params(labelsize=6.5); a2.legend(fontsize=5.8, loc="upper left")
    a2.grid(True, which="both", ls=":", alpha=0.4)
    fig.tight_layout(w_pad=0.6); fig.savefig(path); plt.close(fig)


def fig_throughput(path, tsv):
    rows = _rows(_need(tsv))
    h, g2g, mm2 = _col(rows, "h"), _col(rows, "g2g_s"), _col(rows, "baseline_s")
    o = np.argsort(h); h, g2g, mm2 = h[o], g2g[o], mm2[o]
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    ax.semilogy(h, g2g, "o-", color=BLUE, lw=1.8, ms=4, label="PanG2G (single-threaded)")
    ax.semilogy(h, mm2, "s--", color=RED, lw=1.3, ms=4, label="linearize + map (minimap2)")
    # crossover: first h where the baseline overtakes PanG2G
    cross = next((hh for hh, a, b in zip(h, g2g, mm2) if b > a), None)
    if cross is not None:
        ax.axvline(cross, color=GREY, ls=":", lw=0.8)
        ax.text(cross * 1.15, min(g2g) * 0.8, "crossover", fontsize=6.5, color=GREY)
    ax.set_xlabel("haplotypes per population, $h$"); ax.set_ylabel("wall time (s)")
    ax.annotate(f"{mm2[-1]/g2g[-1]:.1f}$\\times$ at $h$={int(h[-1])}", (h[-1], mm2[-1]),
                ha="right", va="top", fontsize=7, color=RED, xytext=(-5, -25), textcoords="offset points")
    ax.legend(loc="upper left", fontsize=6); ax.grid(True, which="both", ls=":", alpha=0.4)
    fig.savefig(path); plt.close(fig)


def fig_scaling(path, tsv, max_cores=64):
    rows = [r for r in _rows(_need(tsv)) if r["threads"].strip() != "seq"]
    thr = _col(rows, "threads"); e2e = _col(rows, "speedup_e2e"); dp = _col(rows, "speedup_dp")
    keep = thr <= max_cores; thr, e2e, dp = thr[keep], e2e[keep], dp[keep]
    diverged = any(str(r.get("consistent", "")).strip() not in ("ok", "") for r in rows)
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    ax.plot(thr, thr, ":", color="0.6", label="linear")
    ax.plot(thr, dp, "o-", color=BLUE, label="DP phase")
    ax.plot(thr, e2e, "s-", color=RED, label="end-to-end")
    ax.set_xlabel("cores"); ax.set_ylabel("speedup vs. sequential")
    ax.legend(loc="upper left"); ax.grid(True, which="both", ls=":", alpha=0.4)
    if diverged:
        ax.text(0.5, 0.5, r"INVALID — parallel $\neq$ sequential", transform=ax.transAxes,
                ha="center", va="center", fontsize=13, color=RED, alpha=0.3, rotation=18)
    fig.savefig(path); plt.close(fig)


def fig_ablation(path, tsv):
    """Either the raw 08_ablation.sh TSV (ctrl_base, ctrl_disabled, trt_base,
    trt_disabled) or the committed ratio table (ctrl_ratio, trt_ratio)."""
    rows = _rows(_need(tsv))
    if "ctrl_ratio" in rows[0]:
        ctrl, treat = _col(rows, "ctrl_ratio"), _col(rows, "trt_ratio")
    else:
        ctrl = _col(rows, "ctrl_disabled") / _col(rows, "ctrl_base")
        treat = _col(rows, "trt_disabled") / _col(rows, "trt_base")
    labels = ["LAPJV\n$\\to$ brute", "minimizer\n$\\to$ naive",
              "priority-\ndispatch\n$\\to$ static", "LPT\n$\\to$ FIFO", "vEB\n$\\to$ BFS"][:len(rows)]
    fig, ax = plt.subplots(figsize=(3.4, 2.9))
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


def fig_biology(path, tsv):
    rows = _rows(_need(tsv))
    def short(cat):
        c = cat.lower()
        if "f_st" in c or "fst" in c or "gnomad" in c: return "gnomAD-SV\nstructural"
        if "clinvar" in c: return "ClinVar\npathogenic"
        return cat
    cats = [short(r["category"]) for r in rows]
    folds, ps = _col(rows, "fold"), _col(rows, "p")
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    x = np.arange(len(cats))
    ax.bar(x, folds, 0.6, color=[RED if p < 0.05 else GREY for p in ps], alpha=0.85)
    ax.axhline(1, ls="--", color=GREY, lw=0.9, label="null (no enrichment)")
    ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=7.5)
    ax.set_ylabel("fold enrichment (obs / permuted)")
    ax.set_title("Enrichment near SV & clinical sites")
    ax.set_ylim(0, max(folds.max(), 1.0) * 1.42)
    for xi, fo, p in zip(x, folds, ps):
        tag = f"{fo:.2f}×\n$p$={p:.1e}" if p < 1e-3 else f"{fo:.2f}×\n$p$={p:.3f}"
        ax.annotate(tag, (xi, fo), ha="center", va="bottom", fontsize=6.8)
    ax.legend(fontsize=7, loc="upper right"); ax.grid(True, axis="y", ls=":", alpha=0.4)
    fig.savefig(path); plt.close(fig)


def fig_fst_sweep(path, tsv):
    rows = _rows(_need(tsv))
    t, fold, p = _col(rows, "fst_threshold"), _col(rows, "fold"), _col(rows, "p")
    o = np.argsort(t); t, fold, p = t[o], fold[o], p[o]
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    ax.plot(t, fold, "-", color=BLUE, lw=1.6, zorder=1)
    sig = p < 0.05
    ax.scatter(t[sig], fold[sig], c=RED, s=28, zorder=2, label="$p<0.05$")
    if (~sig).any():
        ax.scatter(t[~sig], fold[~sig], c=GREY, s=28, zorder=2, label="n.s.")
    ax.axhline(1, ls="--", color=GREY, lw=0.9)
    ax.set_xlabel("gnomAD-SV $F_{ST}$ threshold"); ax.set_ylabel("fold enrichment (obs / permuted)")
    ax.set_title("Delta enrichment vs. SV differentiation")
    ax.set_ylim(bottom=min(0.9, np.nanmin(fold) * 0.95))
    ax.grid(True, ls=":", alpha=0.4); ax.legend(fontsize=7, loc="upper left")
    fig.savefig(path); plt.close(fig)


def fig_pairing_error(path, csv_path):
    """Snarl-pairing error vs. graph size (Exp1b): 1 - F1 of the recovered
    delta, for the boundary-key rule and the anchor-interval DP, on planted
    graphs with short boundary nodes and 10% of snarls private to one graph.
    Lines: mean over seeds; markers at each measured size."""
    rows = _rows(_need(csv_path))
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    styles = {("boundary", 1): ("s--", RED, "boundary-key, 1-bp boundaries"),
              ("boundary", 2): ("^--", "#e8a0a0", "boundary-key, 2-bp boundaries"),
              ("boundary", 8): ("v--", GREY, "boundary-key, 8-bp boundaries"),
              ("interval", 1): ("o-", BLUE, "anchor-interval DP, 1-bp boundaries")}
    drew = False
    for (mode, L), (st, col, lab) in styles.items():
        sel = [r for r in rows if r["pairing"] == mode and int(r["anchor_len"]) == L
               and float(r["private_frac"]) > 0]
        if not sel:
            continue
        by_n = {}
        for r in sel:
            f1 = 0.5 * (float(r["f1_a"]) + float(r["f1_b"]))
            by_n.setdefault(int(r["n_snarls"]), []).append(1.0 - f1)
        n = np.array(sorted(by_n)); err = np.array([np.mean(by_n[k]) for k in n])
        ax.plot(n, 100 * err, st, color=col, ms=4, lw=1.4, label=lab); drew = True
    if not drew:
        raise SystemExit("pairing-error CSV has no private-snarl rows")
    ax.set_xscale("log")
    ax.set_xlabel("snarls per graph"); ax.set_ylabel("delta error, $1-F_1$ (%)")
    ax.set_title("Pairing error vs. size (10% private snarls)", fontsize=8)
    ax.grid(True, which="both", ls=":", alpha=0.4); ax.legend(fontsize=6, loc="center left")
    fig.savefig(path); plt.close(fig)


# ── schematics ────────────────────────────────────────────────────────────────
def fig_paradigm(path):
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 1.9))
    for ax in axs: _clean(ax); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    axs[0].plot([0.1, 0.9], [0.62, 0.62], "-", color=BLUE, lw=3)
    axs[0].plot([0.1, 0.9], [0.38, 0.38], "-", color=RED, lw=3)
    for x in np.linspace(0.18, 0.82, 5):
        axs[0].plot([x, x], [0.6, 0.4], color=GREY, lw=0.6)
    axs[0].set_title("(a) sequence–sequence")
    axs[1].plot([0.1, 0.9], [0.68, 0.68], "-", color=BLUE, lw=3)
    gx = [0.18, 0.36, 0.5, 0.5, 0.64, 0.82]; gy = [0.32, 0.32, 0.42, 0.22, 0.32, 0.32]
    for i in range(len(gx) - 1):
        axs[1].plot([gx[i], gx[i+1]], [gy[i], gy[i+1]], color=RED, lw=2)
    for x, y in zip(gx, gy): axs[1].add_patch(Circle((x, y), 0.022, color=RED))
    axs[1].set_title("(b) sequence–graph")
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
    xs = [0.08, 0.26, 0.5, 0.74, 0.92]
    for x in xs: a1.add_patch(Circle((x, 0.5), 0.025, color=BLUE))
    a1.plot(xs, [0.5]*5, color=BLUE, lw=1.5, zorder=0)
    for x0, x1 in [(0.26, 0.5), (0.5, 0.74)]:
        a1.plot([x0, (x0+x1)/2, x1], [0.5, 0.66, 0.5], color=RED, lw=1.6)
        a1.plot([x0, (x0+x1)/2, x1], [0.5, 0.34, 0.5], color=RED, lw=1.6)
    a1.add_patch(Ellipse((0.5, 0.5), 0.62, 0.5, fill=False, ec=GREY, ls="--", lw=1))
    a1.set_title("(a) variation graph $\\to$ nested snarls")
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
    stages = ["detect anchors\n(minimizer index)", "decompose\n(snarl trees)",
              "pair snarls", "match alleles\n(LAPJV)", "emit deltas"]
    parallel = [False, False, False, True, False]
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
                _box(ax, x0 + off[0], y + off[1], w, h, "", fc="#eef3fb", ec=BLUE, fs=7, lw=0.6)
                ax.patches[-1].set_alpha(0.35)
        _box(ax, x0, y, w, h, stages[i], fs=7.5)
    for i in range(n - 1):
        _arrow(ax, 0.5, ys[i], 0.5, ys[i + 1] + h)
    ax.text(0.5, y_top + 0.06, "GFA $G_A$, $G_B$", ha="center", fontsize=7.5, color=BLUE)
    _arrow(ax, 0.5, y_top + 0.03, 0.5, y_top)
    ax.text(0.5, ys[-1] - 0.06, r"$\Delta_{A\to B},\ \Delta_{B\to A}$, alignment (GAF), score",
            ha="center", fontsize=7.5, color=GREEN)
    _arrow(ax, 0.5, ys[-1], 0.5, ys[-1] - 0.03)
    ax.text(0.88, ys[3] + h / 2, "parallel", ha="center", va="center", fontsize=6.5, color=GREY, rotation=90)
    ax.annotate("", xy=(0.82, ys[3]), xytext=(0.82, ys[3] + h),
                arrowprops=dict(arrowstyle="<->", color=GREY, lw=0.8))
    fig.savefig(path); plt.close(fig)


def fig_delta(path):
    DARK = "#444"
    y0, up, dn = 0.50, 0.82, 0.18
    bx = list(np.linspace(0.10, 0.90, 6))
    spine_seq = ["GA", "TT", "CG", "GC", "TA", "CC"]
    fig, axes = plt.subplots(3, 1, figsize=(3.4, 3.6))
    for ax in axes:
        _clean(ax); ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    def node(ax, x, y, color):
        ax.plot(x, y, "o", color=color, ms=8.5, mec="none", zorder=3)

    def backbone(ax, label=False):
        for i in (0, 2, 4):
            ax.plot([bx[i], bx[i + 1]], [y0, y0], color=DARK, lw=1.3, zorder=1)
        for x, s in zip(bx, spine_seq):
            node(ax, x, y0, DARK)
            if label:
                ax.text(x, y0 - 0.19, s, ha="center", fontsize=5, color=DARK)

    def shared(ax, xl, xr, seq):
        xm = (xl + xr) / 2
        ax.plot([xl, xm, xr], [y0, y0, y0], color=GREEN, lw=1.7, zorder=2)
        node(ax, xm, y0, GREEN)
        ax.text(xm, y0 + 0.15, seq, ha="center", fontsize=5.5, color=GREEN, style="italic")

    def allele(ax, xl, xr, yapex, seq, color, above):
        xm = (xl + xr) / 2
        ax.plot([xl, xm, xr], [y0, yapex, y0], color=color, lw=1.4, zorder=1)
        node(ax, xm, yapex, color)
        ax.text(xm, yapex + (0.16 if above else -0.16), seq, ha="center", fontsize=5.5, color=color, style="italic")

    a = axes[0]; backbone(a, label=True)
    shared(a, bx[1], bx[2], "AC");  allele(a, bx[1], bx[2], up, "ATC",  BLUE, True)
    shared(a, bx[3], bx[4], "T");   allele(a, bx[3], bx[4], up, "TAGG", BLUE, True)
    a.text(0.0, 0.98, r"(a) input graph $A$", fontsize=6.5, va="top", color="#222")
    b = axes[1]; backbone(b)
    shared(b, bx[1], bx[2], "AC");  allele(b, bx[1], bx[2], dn, "A", RED, False)
    shared(b, bx[3], bx[4], "T")
    b.text(0.0, 0.98, r"(b) input graph $B$", fontsize=6.5, va="top", color="#222")
    c = axes[2]; backbone(c)
    shared(c, bx[1], bx[2], "AC")
    allele(c, bx[1], bx[2], up, "ATC", BLUE, True)
    allele(c, bx[1], bx[2], dn, "A",   RED,  False)
    shared(c, bx[3], bx[4], "T")
    allele(c, bx[3], bx[4], up, "TAGG", BLUE, True)
    c.text(0.0, 0.98, r"(c) delta $\Delta_{A\to B}$, $\Delta_{B\to A}$", fontsize=6.5, va="top", color="#222")
    fig.text(0.5, 0.375, r"$\Downarrow$ diff", ha="center", va="center", fontsize=8)
    for col, lab in [(DARK, "backbone"), (GREEN, "shared allele (in both)"),
                     (BLUE, r"$\Delta_{A\to B}$ ($A$-only)"), (RED, r"$\Delta_{B\to A}$ ($B$-only)")]:
        c.plot([], [], color=col, lw=1.6, label=lab)
    c.legend(loc="upper center", fontsize=5.6, ncol=2, frameon=False,
             bbox_to_anchor=(0.5, -0.02), handlelength=1.4, columnspacing=1.2)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.99, bottom=0.09, hspace=0.22)
    fig.savefig(path); plt.close(fig)


def fig_veb(path):
    fig, axes = plt.subplots(3, 1, figsize=(3.4, 3.6), gridspec_kw={"height_ratios": [1.0, 0.5, 0.5]})
    for ax in axes: _clean(ax)
    a0 = axes[0]; a0.set_xlim(0, 1); a0.set_ylim(0, 1)
    tree = {"R": (0.50, 0.90), "S1": (0.25, 0.65), "S2": (0.75, 0.65),
            "a": (0.13, 0.40), "b": (0.37, 0.40), "c": (0.63, 0.40), "S3": (0.87, 0.40),
            "d": (0.78, 0.16), "e": (0.94, 0.16)}
    edges = [("R","S1"),("R","S2"),("S1","a"),("S1","b"),("S2","c"),("S2","S3"),("S3","d"),("S3","e")]
    snarl_nodes = {"R", "S1", "S2", "S3"}
    st_colors = {"S1": "#cde4f7", "S3": "#fde0d0"}
    a0.add_patch(FancyBboxPatch((0.06, 0.30), 0.38, 0.48, boxstyle="round,pad=0.03", fc=st_colors["S1"], ec="none", alpha=0.4, zorder=0))
    a0.add_patch(FancyBboxPatch((0.70, 0.06), 0.30, 0.48, boxstyle="round,pad=0.03", fc=st_colors["S3"], ec="none", alpha=0.4, zorder=0))
    for p, c in edges:
        px, py = tree[p]; cx, cy = tree[c]
        a0.plot([px, cx], [py, cy], color="#666", lw=1.0, zorder=1)
    for name, (x, y) in tree.items():
        fc = "#eef3fb" if name in snarl_nodes else "#f5f5f5"
        ec = BLUE if name in snarl_nodes else GREY
        a0.add_patch(Circle((x, y), 0.05, fc=fc, ec=ec, lw=1.2, zorder=2))
        a0.text(x, y, name, ha="center", va="center", fontsize=6.5,
                fontweight="bold" if name in snarl_nodes else "normal", zorder=3)
    a0.set_title("(a) nested snarl tree", fontsize=8, pad=2)

    def cells(ax, order, highlights, title):
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        n = len(order); cw = 0.85 / n; x0 = (1 - n * cw) / 2
        for i, lab in enumerate(order):
            fc = highlights.get(lab, "#f5f5f5")
            ec_c = BLUE if lab in snarl_nodes else GREY
            ax.add_patch(Rectangle((x0 + i * cw, 0.25), cw * 0.92, 0.50, fc=fc, ec=ec_c, lw=1))
            ax.text(x0 + i * cw + cw * 0.46, 0.50, lab, ha="center", va="center", fontsize=6,
                    fontweight="bold" if lab in snarl_nodes else "normal")
        ax.set_title(title, fontsize=7, pad=2)

    bfs_order = ["R", "S1", "S2", "a", "b", "c", "S3", "d", "e"]
    veb_order = ["R", "S1", "a", "b", "S2", "c", "S3", "d", "e"]
    hi = {**{n: st_colors["S1"] for n in ["S1", "a", "b"]}, **{n: st_colors["S3"] for n in ["S3", "d", "e"]},
          "R": "#eef3fb", "S2": "#eef3fb", "c": "#f5f5f5"}
    cells(axes[1], bfs_order, hi, "(b) BFS layout — subtrees scattered")
    cells(axes[2], veb_order, hi, "(c) vEB layout — subtrees contiguous")
    a2 = axes[2]
    a2.annotate("", xy=(0.205, 0.15), xytext=(0.41, 0.15), arrowprops=dict(arrowstyle="|-|", color=BLUE, lw=0.8))
    a2.text(0.305, 0.03, "S1 subtree", ha="center", fontsize=5.5, color=BLUE)
    a2.annotate("", xy=(0.68, 0.15), xytext=(0.89, 0.15), arrowprops=dict(arrowstyle="|-|", color=RED, lw=0.8))
    a2.text(0.79, 0.03, "S3 subtree", ha="center", fontsize=5.5, color=RED)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


# ── driver ────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=R("latex", "figs"))
    ap.add_argument("--matching-csv",   default=R("results", "matching_cost_vs_k_cpp.csv"))
    ap.add_argument("--gini-csv",       default=R("results", "work_distribution.csv"))
    ap.add_argument("--alleles-csv",    default=R("results", "allele_count_hist.csv"))
    ap.add_argument("--refdiff-tsv",    default=R("results", "refdiff_scaling.tsv"))
    ap.add_argument("--throughput-csv", default=R("demo", "sample_outputs", "throughput_chr22.tsv"))
    ap.add_argument("--scaling-csv",    default=R("results", "ibex", "parallel_scaling_128.tsv"))
    ap.add_argument("--submitted-scaling-csv", default=R("results", "parallel_scaling.tsv"))
    ap.add_argument("--ablation-csv",   default=R("results", "accl", "ablation_results.tsv"))
    ap.add_argument("--biology-csv",    default=R("demo", "sample_outputs", "biology_results.tsv"))
    ap.add_argument("--fst-sweep-csv",  default=R("demo", "sample_outputs", "fst_sweep_chr22.tsv"))
    ap.add_argument("--pairing-csv",    default=R("results", "accl", "exp1b_pairing_error.csv"))
    ap.add_argument("--only", default=None, help="comma-separated figure names")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    figures = {
        "fig_paradigm":       lambda p: fig_paradigm(p),
        "fig_snarltree":      lambda p: fig_snarltree(p),
        "fig_pipeline":       lambda p: fig_pipeline(p),
        "fig_delta":          lambda p: fig_delta(p),
        "fig_veb":            lambda p: fig_veb(p),
        "fig_matching":       lambda p: fig_matching(p, a.matching_csv),
        "fig_workload":       lambda p: fig_workload(p, a.gini_csv, a.alleles_csv),
        "fig_refdiffscaling": lambda p: fig_refdiffscaling(p, a.refdiff_tsv),
        "fig_throughput":     lambda p: fig_throughput(p, a.throughput_csv),
        "fig_scaling":        lambda p: fig_scaling(p, a.scaling_csv),
        "fig_ablation":       lambda p: fig_ablation(p, a.ablation_csv),
        "fig_biology":        lambda p: fig_biology(p, a.biology_csv),
        "fig_fst_sweep":      lambda p: fig_fst_sweep(p, a.fst_sweep_csv),
        "fig_pairing_error":  lambda p: fig_pairing_error(p, a.pairing_csv),
        # combined figures used in the camera-ready layout
        "fig_workload_matching": lambda p: fig_workload_matching(p, a.gini_csv, a.matching_csv),
        "fig_perf":           lambda p: fig_perf(p, a.scaling_csv, a.throughput_csv),
    }
    only = set(a.only.split(",")) if a.only else None
    n = 0
    for name, fn in figures.items():
        if only and name not in only:
            continue
        p = os.path.join(a.outdir, name + ".pdf")
        fn(p); print("wrote", p); n += 1
    print(f"\n{n} figures -> {a.outdir}/")


if __name__ == "__main__":
    main()
