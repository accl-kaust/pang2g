"""
profile_workload.py - Phase 2: Workload Characterization
=========================================================
Measures the four quantities called for in the project outline:

  1. Wall-clock time breakdown (snarl decomposition vs bipartite matching
     vs anchor detection) across graph sizes n ∈ {10, 50, 100, 500, 1000}.

  2. Memory usage peak at each graph size (tracemalloc).

  3. Work distribution: record k (allele count) for every snarl pair.
     Compute the Gini coefficient - high Gini = high imbalance =
     strong motivation for work-stealing over static partitioning.

  4. Brute-force matching cost vs k: measure time for
     bipartite_allele_match at k = 2, 3, 4, 5, 6, 7, 8 to find
     the crossover point where brute-force becomes infeasible.

Outputs CSV files and a summary table to stdout.
Run with:  python profiling/profile_workload.py
"""

from __future__ import annotations
import sys, os, time, tracemalloc, cProfile, pstats, io, csv



from synth_generator import generate_graph_pair
from g2g.decompose   import SnarlDecomposer
from g2g.align       import GraphAligner, bipartite_allele_match
from g2g.sequence    import seq_similarity
import synth_generator as _sg
_orig_zipf = _sg._zipf_k
_sg._zipf_k = lambda rng, k_min=2, k_max=6, alpha=2.2: _orig_zipf(rng, k_min=2, k_max=6, alpha=2.2)

# -- Config ----------------------------------------------------
N_SIZES   = [10, 50, 100, 200, 500]
K_VALUES  = [2, 3, 4, 5, 6, 7, 8]
SEQ_LEN   = 12      # representative inner-node sequence length for k-sweep
REPEATS   = 3       # repeats per size for stable timing
OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)


# -------------------------------------------------------------
# Helpers
# -------------------------------------------------------------

def gini(values: list[float]) -> float:
    """Gini coefficient of a list of non-negative values."""
    if not values or sum(values) == 0:
        return 0.0
    v = sorted(values)
    n = len(v)
    cumsum = 0.0
    for i, x in enumerate(v):
        cumsum += (2 * (i + 1) - n - 1) * x
    return cumsum / (n * sum(v))


def time_phase(fn, *args, **kwargs):
    """Return (elapsed_seconds, return_value)."""
    t0 = time.perf_counter()
    result = fn(*args, **kwargs)
    return time.perf_counter() - t0, result


# -------------------------------------------------------------
# Experiment 1 & 2: Timing + Memory by graph size
# -------------------------------------------------------------

def run_timing_and_memory() -> list[dict]:
    print("\n" + "=" * 64)
    print("  EXPERIMENT 1+2: Timing & Memory vs Graph Size")
    print("=" * 64)
    print(f"  {'n':>6}  {'decomp(s)':>10}  {'align(s)':>10}  "
          f"{'total(s)':>10}  {'mem_peak(MB)':>13}  {'n_snarls':>9}")
    print("  " + "-" * 62)

    rows = []
    for n in N_SIZES:
        t_decomp_all, t_align_all, mem_all = [], [], []
        for rep in range(REPEATS):
            ga, gb = generate_graph_pair(n_snarls=n, seed=rep)

            # Memory tracking
            tracemalloc.start()

            t_d, (sa, sb) = time_phase(
                lambda: (SnarlDecomposer(ga).decompose(),
                         SnarlDecomposer(gb).decompose())
            )

            aligner = GraphAligner()
            t_a, result = time_phase(aligner.align, ga, gb)

            _, mem_peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            t_decomp_all.append(t_d)
            t_align_all.append(t_a)
            mem_all.append(mem_peak / 1024 / 1024)

        t_d   = sum(t_decomp_all) / REPEATS
        t_a   = sum(t_align_all)  / REPEATS
        mem   = sum(mem_all)      / REPEATS
        total = t_d + t_a

        # Snarl count from last run
        sa = SnarlDecomposer(ga).decompose()
        sb = SnarlDecomposer(gb).decompose()
        n_snarls = len(sa) + len(sb)

        row = dict(n=n, t_decomp=t_d, t_align=t_a, t_total=total,
                   mem_mb=mem, n_snarls=n_snarls)
        rows.append(row)

        print(f"  {n:>6}  {t_d:>10.4f}  {t_a:>10.4f}  "
              f"{total:>10.4f}  {mem:>13.2f}  {n_snarls:>9}")

    return rows


# -------------------------------------------------------------
# Experiment 3: Work distribution & Gini coefficient
# -------------------------------------------------------------

def run_work_distribution() -> list[dict]:
    print("\n" + "=" * 64)
    print("  EXPERIMENT 3: Work Distribution & Gini Coefficient")
    print("=" * 64)

    rows = []
    for n in N_SIZES:
        ga, gb = generate_graph_pair(n_snarls=n, seed=42)
        snarls_a = SnarlDecomposer(ga).decompose()
        snarls_b = SnarlDecomposer(gb).decompose()

        # Work per snarl ≈ k³ (Hungarian algorithm cost)
        all_k: list[int]     = []
        all_work: list[float] = []

        for s in snarls_a + snarls_b:
            k = len(s.allele_paths)
            all_k.append(k)
            all_work.append(k ** 3)

        gini_k    = gini([float(x) for x in all_k])
        gini_work = gini(all_work)
        k_max     = max(all_k)
        k_mean    = sum(all_k) / len(all_k)

        # Distribution buckets
        hist = {}
        for k in all_k:
            hist[k] = hist.get(k, 0) + 1

        # % of snarls that are SNP-bubbles (k=2)
        pct_snp = 100.0 * hist.get(2, 0) / len(all_k)

        row = dict(n=n, total_snarls=len(all_k), k_mean=k_mean, k_max=k_max,
                   gini_k=gini_k, gini_work=gini_work, pct_snp=pct_snp,
                   hist=dict(sorted(hist.items())))
        rows.append(row)

        print(f"\n  n={n}: {len(all_k)} snarls  "
              f"k_mean={k_mean:.2f}  k_max={k_max}  "
              f"Gini(k)={gini_k:.3f}  Gini(k³)={gini_work:.3f}  "
              f"SNP-bubbles={pct_snp:.1f}%")
        print(f"  Allele-count histogram: {dict(sorted(hist.items()))}")

    return rows


# -------------------------------------------------------------
# Experiment 4: Brute-force matching cost vs k
# -------------------------------------------------------------

def _make_seqs(k: int, length: int, seed: int = 0) -> list[str]:
    import random
    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(length)) for _ in range(k)]


def run_matching_cost_vs_k() -> list[dict]:
    print("\n" + "=" * 64)
    print("  EXPERIMENT 4: Bipartite Matching Cost vs k")
    print("=" * 64)
    print(f"  {'k':>4}  {'time_ms':>10}  {'n_pairs_tried':>15}  {'speedup_vs_k2':>15}")
    print("  " + "-" * 46)

    rows = []
    base_time = None

    for k in K_VALUES:
        seqs_a = _make_seqs(k, SEQ_LEN, seed=1)
        seqs_b = _make_seqs(k, SEQ_LEN, seed=2)

        # Time over multiple calls for stability
        reps = max(1, 200 // (k ** 3))
        t0 = time.perf_counter()
        for _ in range(reps):
            bipartite_allele_match(seqs_a, seqs_b, gap_cost=0.5, min_match_sim=0.0)
        elapsed = (time.perf_counter() - t0) / reps * 1000  # ms

        if base_time is None:
            base_time = elapsed if elapsed > 0 else 1e-9

        # Theoretical pairs tried: sum over m of C(k,m)*P(k,m)
        from math import comb, perm
        pairs = sum(comb(k, m) * perm(k, m) for m in range(k + 1))

        speedup = elapsed / base_time
        row = dict(k=k, time_ms=elapsed, pairs_tried=pairs, speedup=speedup)
        rows.append(row)

        print(f"  {k:>4}  {elapsed:>10.3f}  {pairs:>15,d}  {speedup:>15.1f}x")

    return rows


# -------------------------------------------------------------
# cProfile hotspot dump for n=500
# -------------------------------------------------------------

def run_cprofile_hotspots():
    print("\n" + "=" * 64)
    print("  cPROFILE: Top 15 hotspots at n=500")
    print("=" * 64)

    ga, gb = generate_graph_pair(n_snarls=500, seed=42)
    aligner = GraphAligner()

    pr = cProfile.Profile()
    pr.enable()
    aligner.align(ga, gb)
    pr.disable()

    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
    ps.print_stats(15)
    output = s.getvalue()

    # Filter to only lines with g2g or profiling in the path
    lines = output.split("\n")
    header_done = False
    for line in lines:
        if "ncalls" in line or "cumtime" in line:
            header_done = True
        if header_done:
            print("  " + line)


# -------------------------------------------------------------
# Save CSV results
# -------------------------------------------------------------

def save_csv(filename: str, rows: list[dict], exclude_keys: list[str] = None):
    exclude_keys = exclude_keys or []
    path = os.path.join(OUTPUT_DIR, filename)
    keys = [k for k in rows[0].keys() if k not in exclude_keys]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in rows:
            w.writerow({k: row[k] for k in keys})
    print(f"\n  → Saved {path}")


# -------------------------------------------------------------
# Main
# -------------------------------------------------------------

if __name__ == "__main__":
    timing_rows  = run_timing_and_memory()
    work_rows    = run_work_distribution()
    matching_rows = run_matching_cost_vs_k()
    run_cprofile_hotspots()

    save_csv("timing_memory.csv",   timing_rows)
    save_csv("work_distribution.csv", work_rows, exclude_keys=["hist"])
    save_csv("matching_cost_vs_k.csv", matching_rows)

    print("\n" + "=" * 64)
    print("  Phase 2 profiling complete.")
    print("  Results written to:", OUTPUT_DIR)
    print("=" * 64)
