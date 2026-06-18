#!/usr/bin/env python3
"""
plot_exp6_biology.py — Biological validation figure (Figure 6 in paper).

Generates:
  - gnomAD-SV enrichment: observed vs expected overlap, enrichment ratio
  - ClinVar cross-population conservation: shared alleles vs pathogenic SVs
  - Population-specific delta alleles vs high-FST gnomAD SVs
  - Summary bar charts

Usage:
    python experiments/plot/plot_exp6_biology.py \\
        --gnomad-all  results/phase5/exp6/g2g_delta_vs_gnomad_all.bed \\
        --gnomad-fst  results/phase5/exp6/g2g_delta_vs_gnomad_highfst.bed \\
        --clinvar     results/phase5/exp6/g2g_shared_vs_clinvar_patho.bed \\
        --out         results/phase5/figures/
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


# ── BED field names from bedtools intersect output ──────────────────
BED_COLS_A = ["chrom_a", "start_a", "end_a", "name_a"]
BED_COLS_B = ["chrom_b", "start_b", "end_b", "name_b"]
INTERSECT_COLS = BED_COLS_A + BED_COLS_B


def load_intersect_bed(path: str) -> pd.DataFrame:
    """Load bedtools intersect output (8-column or more)."""
    if path is None or not os.path.exists(path):
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t", header=None,
                         names=INTERSECT_COLS[:8], usecols=range(8),
                         comment="#")
        return df
    except Exception as e:
        print(f"WARNING: {path}: {e}", file=sys.stderr)
        return pd.DataFrame()


def load_summary_tsv(path: str) -> pd.DataFrame:
    """Load a tab-separated summary file if it exists."""
    if path is None or not os.path.exists(path):
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", comment="#")
    except Exception as e:
        print(f"WARNING: {path}: {e}", file=sys.stderr)
        return pd.DataFrame()


def count_unique_alleles(df: pd.DataFrame, col: str = "name_a") -> int:
    return df[col].nunique() if not df.empty and col in df.columns else 0


def plot_biology(args, out_dir: str) -> None:
    df_gnomad_all  = load_intersect_bed(args.gnomad_all)
    df_gnomad_fst  = load_intersect_bed(args.gnomad_fst)
    df_clinvar     = load_intersect_bed(args.clinvar)
    df_summary     = load_summary_tsv(args.summary) if hasattr(args, "summary") and args.summary else pd.DataFrame()

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # ── Panel A: gnomAD-SV enrichment ──────────────────────────
    ax = axes[0, 0]
    _plot_gnomad_enrichment(ax, df_gnomad_all, df_gnomad_fst, args)

    # ── Panel B: ClinVar conservation ──────────────────────────
    ax = axes[0, 1]
    _plot_clinvar_overlap(ax, df_clinvar, args)

    # ── Panel C: Population-specific delta alleles breakdown ───
    ax = axes[1, 0]
    _plot_delta_breakdown(ax, df_gnomad_all, df_gnomad_fst)

    # ── Panel D: Summary enrichment bar chart ──────────────────
    ax = axes[1, 1]
    _plot_summary_bars(ax, df_gnomad_all, df_gnomad_fst, df_clinvar, args)

    fig.suptitle("G2G Biological Validation — Experiment 6\n(AFR n=20 vs EAS n=20 HPRC pangenomes)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()

    for ext in ["pdf", "png"]:
        path = os.path.join(out_dir, f"fig6_biology.{ext}")
        dpi = 150 if ext == "pdf" else 200
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        print(f"[plot] Saved: {path}")
    plt.close(fig)


def _plot_gnomad_enrichment(ax, df_all: pd.DataFrame, df_fst: pd.DataFrame, args) -> None:
    """Observed vs expected gnomAD-SV overlap with enrichment fold."""
    categories = ["All gnomAD SVs", "High-FST SVs\n(AFR-EAS FST > 0.10)"]
    observed = [
        count_unique_alleles(df_all),
        count_unique_alleles(df_fst),
    ]
    expected = [
        getattr(args, "expected_all", None),
        getattr(args, "expected_fst", None),
    ]

    x = np.arange(len(categories))
    width = 0.35
    colors_obs = ["#2166ac", "#d73027"]
    colors_exp = ["#92c5de", "#f4a582"]

    if all(v is not None for v in expected):
        bars_obs = ax.bar(x - width/2, observed, width, color=colors_obs, label="Observed")
        bars_exp = ax.bar(x + width/2, expected, width, color=colors_exp, label="Expected (shuffled)")
        for i, (obs, exp) in enumerate(zip(observed, expected)):
            if exp and exp > 0:
                fold = obs / exp
                ymax = max(obs, exp)
                ax.text(i, ymax * 1.05, f"{fold:.2f}×",
                        ha="center", fontsize=10, fontweight="bold",
                        color="#d73027" if fold > 1 else "#4575b4")
        ax.legend(fontsize=9)
    else:
        bars = ax.bar(x, observed, 0.5, color=colors_obs)
        for bar, val in zip(bars, observed):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    str(val), ha="center", va="bottom", fontsize=10)

    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=10)
    ax.set_ylabel("G2G delta alleles overlapping gnomAD SVs")
    ax.set_title("(A) gnomAD-SV Enrichment\n(population-divergent delta alleles)")
    ax.grid(True, alpha=0.3, axis="y")

    if df_all.empty and df_fst.empty:
        ax.text(0.5, 0.5, "No intersection data\n(run Exp6 first)",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=11, color="gray",
                bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))


def _plot_clinvar_overlap(ax, df_clinvar: pd.DataFrame, args) -> None:
    """Shared alleles overlapping ClinVar pathogenic SVs."""
    n_shared_clinvar = count_unique_alleles(df_clinvar)
    n_total_shared   = getattr(args, "n_shared_total", None)

    if n_total_shared and n_total_shared > 0:
        pct = 100.0 * n_shared_clinvar / n_total_shared
        wedge_sizes  = [n_shared_clinvar, n_total_shared - n_shared_clinvar]
        wedge_labels = [f"ClinVar\npathogenic\nn={n_shared_clinvar}", "Other shared\nalleles"]
        wedge_colors = ["#d73027", "#92c5de"]
        wedges, texts, autotexts = ax.pie(
            wedge_sizes, labels=wedge_labels, colors=wedge_colors,
            autopct="%1.1f%%", startangle=90,
            textprops={"fontsize": 9},
            wedgeprops={"edgecolor": "white", "linewidth": 1.5},
        )
        for at in autotexts:
            at.set_fontsize(9)
        ax.set_title(f"(B) Shared Alleles vs. ClinVar\n"
                     f"(n={n_total_shared:,} shared; {pct:.1f}% overlap pathogenic SVs)")
    else:
        # Fallback bar plot when total not provided
        vals   = [n_shared_clinvar, getattr(args, "n_delta_total", 0)]
        labels = ["Shared alleles\n∩ ClinVar pathogenic", "All G2G\nshared alleles"]
        bars = ax.bar(range(len(vals)), vals, color=["#d73027", "#92c5de"], width=0.5)
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels(labels, fontsize=10)
        ax.set_ylabel("Count")
        ax.set_title("(B) Shared Alleles vs. ClinVar Pathogenic SVs")
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    str(v), ha="center", va="bottom", fontsize=10)

    if df_clinvar.empty:
        ax.clear()
        ax.text(0.5, 0.5, "No ClinVar intersection data\n(run Exp6 first)",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=11, color="gray",
                bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
        ax.set_title("(B) Shared Alleles vs. ClinVar Pathogenic SVs")

    ax.grid(False)


def _plot_delta_breakdown(ax, df_all: pd.DataFrame, df_fst: pd.DataFrame) -> None:
    """Population-specific delta alleles by SV type (if available)."""
    if df_all.empty:
        ax.text(0.5, 0.5, "No data\n(run Exp6 first)",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=11, color="gray",
                bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
        ax.set_title("(C) Delta Allele Breakdown by SV Type")
        return

    # Classify by allele name suffix if present (DEL/INS/INV/DUP)
    sv_types = ["DEL", "INS", "INV", "DUP", "OTHER"]
    colors   = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#999999"]

    def classify(df):
        counts = {t: 0 for t in sv_types}
        if df.empty or "name_a" not in df.columns:
            return counts
        names = df["name_a"].dropna().unique()
        for name in names:
            name_u = str(name).upper()
            matched = False
            for sv in ["DEL", "INS", "INV", "DUP"]:
                if sv in name_u:
                    counts[sv] += 1
                    matched = True
                    break
            if not matched:
                counts["OTHER"] += 1
        return counts

    counts_all = classify(df_all)
    counts_fst = classify(df_fst)

    x = np.arange(len(sv_types))
    width = 0.35
    ax.bar(x - width/2, [counts_all[t] for t in sv_types], width,
           color=colors, alpha=0.9, label="All gnomAD overlap")
    ax.bar(x + width/2, [counts_fst[t] for t in sv_types], width,
           color=colors, alpha=0.5, label="High-FST overlap",
           edgecolor="black", linewidth=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(sv_types, fontsize=10)
    ax.set_ylabel("Delta allele count")
    ax.set_title("(C) Delta Alleles by SV Type\n(gnomAD-overlapping)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")


def _plot_summary_bars(ax, df_all: pd.DataFrame, df_fst: pd.DataFrame,
                       df_clinvar: pd.DataFrame, args) -> None:
    """Summary horizontal bar chart of key biological validation metrics."""
    metrics = []
    values  = []
    colors  = []

    n_all   = count_unique_alleles(df_all)
    n_fst   = count_unique_alleles(df_fst)
    n_clv   = count_unique_alleles(df_clinvar)

    metrics.append("Delta alleles\n∩ gnomAD (all SVs)")
    values.append(n_all)
    colors.append("#2166ac")

    metrics.append("Delta alleles\n∩ gnomAD high-FST")
    values.append(n_fst)
    colors.append("#d73027")

    metrics.append("Shared alleles\n∩ ClinVar pathogenic")
    values.append(n_clv)
    colors.append("#4dac26")

    # Enrichment fold if expected provided
    exp_all = getattr(args, "expected_all", None)
    exp_fst = getattr(args, "expected_fst", None)
    if exp_all and exp_all > 0 and n_all > 0:
        metrics.append("Enrichment fold\n(all gnomAD, ×)")
        values.append(round(n_all / exp_all, 2))
        colors.append("#984ea3")
    if exp_fst and exp_fst > 0 and n_fst > 0:
        metrics.append("Enrichment fold\n(high-FST, ×)")
        values.append(round(n_fst / exp_fst, 2))
        colors.append("#ff7f00")

    y = np.arange(len(metrics))
    bars = ax.barh(y, values, color=colors, alpha=0.85, edgecolor="white")
    for bar, val in zip(bars, values):
        ax.text(bar.get_width() + max(values) * 0.01, bar.get_y() + bar.get_height()/2,
                f"{val:,}", va="center", fontsize=9)

    ax.set_yticks(y)
    ax.set_yticklabels(metrics, fontsize=9)
    ax.set_xlabel("Count / Fold enrichment")
    ax.set_title("(D) Biological Validation Summary")
    ax.grid(True, alpha=0.3, axis="x")
    ax.set_xlim(0, max(values) * 1.15 if values else 10)

    if all(v == 0 for v in values):
        ax.text(0.5, 0.5, "No data\n(run Exp6 first)",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=11, color="gray",
                bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--gnomad-all",  default=None,
                        help="bedtools intersect output: delta alleles vs gnomAD all SVs")
    parser.add_argument("--gnomad-fst",  default=None,
                        help="bedtools intersect output: delta alleles vs gnomAD high-FST SVs")
    parser.add_argument("--clinvar",     default=None,
                        help="bedtools intersect output: shared alleles vs ClinVar pathogenic SVs")
    parser.add_argument("--summary",     default=None,
                        help="Optional TSV summary produced by Exp6 scripts")
    parser.add_argument("--expected-all", type=float, default=None,
                        help="Expected overlap count with all gnomAD SVs (from shuffled background)")
    parser.add_argument("--expected-fst", type=float, default=None,
                        help="Expected overlap count with high-FST gnomAD SVs (from shuffled background)")
    parser.add_argument("--n-shared-total", type=int, default=None,
                        help="Total shared alleles from G2G output (for ClinVar pie chart)")
    parser.add_argument("--out", default="results/phase5/figures/",
                        help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    plot_biology(args, args.out)


if __name__ == "__main__":
    main()
