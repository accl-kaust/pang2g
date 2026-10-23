#!/usr/bin/env python3
"""
hprc_enrich.py - population-differentiation enrichment for a G2G delta graph.

The Q4 claim: the alleles G2G reports as present in one population pangenome but
not the other (the delta graph) should be enriched for *population-differentiated*
structural variants - i.e. SVs with a high F_ST between the two populations  - 
relative to chance. This script tests that with a permutation null.

Inputs
  --delta     G2G delta VCF (from `g2g --ref GRCh38 --vcf …`); the A-only/B-only
              records are the population-specific alleles.
  --gnomad-sv gnomAD-SV VCF (or any VCF whose INFO carries per-population allele
              frequencies). F_ST is computed per record from two AF fields.
  --sizes     chrom sizes (UCSC ".genome" or ".fai": "chrom<TAB>length …") - the
              universe the permutation shuffles within.
  --clinvar   (optional) ClinVar VCF; reports delta overlap with pathogenic sites.

Method
  observed = # delta alleles within --window bp of a high-F_ST gnomAD-SV.
  null     = the same count after shuffling every delta to a uniformly random
             position on its chromosome, repeated --perms times.
  fold     = observed / mean(null);  p = (#{null >= observed} + 1)/(perms + 1).

Self-test:  python3 hprc_enrich.py --self-test   (no external data needed).
Dependencies: Python standard library; numpy is used for the permutation loop
when importable (same statistic, ~100x faster on chromosome-scale deltas).
"""
import argparse, bisect, csv, gzip, math, os, random, sys


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_sizes(path):
    sizes = {}
    with _open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2 and p[1].isdigit():
                sizes[p[0]] = int(p[1])
    return sizes


def _info(field):
    d = {}
    for kv in field.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1); d[k] = v
    return d


def read_delta(path, min_len=0):
    """G2G delta VCF -> [(chrom, pos, group)]; group 'A'/'B' from INFO GRAPH=.
    With min_len>0 keep only alleles whose INFO LEN >= min_len (SV-sized deltas),
    the like-for-like set when the targets are structural variants."""
    out = []
    with _open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 8:
                continue
            info = _info(c[7])
            if min_len > 0:
                try:
                    if int(info.get("LEN", "0")) < min_len:
                        continue
                except ValueError:
                    continue
            out.append((c[0], int(c[1]), info.get("GRAPH", "?")))
    return out


def hudson_fst(p1, p2):
    """Hudson's F_ST estimator for one biallelic site between two populations."""
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return ((p1 - p2) ** 2) / den if den > 0 else 0.0


def read_gnomad_records(path, afr_field, eas_field):
    """Parse a gnomAD-SV VCF once -> ([(chrom,pos,F_ST)], n_total, n_usable)."""
    recs = []
    n_total = n_usable = 0
    with _open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 8:
                continue
            n_total += 1
            info = _info(c[7])
            try:
                p1 = float(info[afr_field]); p2 = float(info[eas_field])
            except (KeyError, ValueError):
                continue
            n_usable += 1
            recs.append((c[0], int(c[1]), hudson_fst(p1, p2)))
    return recs, n_total, n_usable


def targets_at(recs, fst_thresh):
    """{chrom: sorted [pos]} for records with F_ST >= fst_thresh (+ count)."""
    hi = {}
    for chrom, pos, fst in recs:
        if fst >= fst_thresh:
            hi.setdefault(chrom, []).append(pos)
    for chrom in hi:
        hi[chrom].sort()
    return hi, sum(len(v) for v in hi.values())


def read_gnomad_fst(path, afr_field, eas_field, fst_thresh):
    """-> ({chrom: sorted [pos]} for high-F_ST SVs, n_total, n_high, n_usable)."""
    recs, n_total, n_usable = read_gnomad_records(path, afr_field, eas_field)
    hi, n_high = targets_at(recs, fst_thresh)
    return hi, n_total, n_high, n_usable


def read_points(path, want_pathogenic=False):
    pts = {}
    with _open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 8:
                continue
            if want_pathogenic and "athogenic" not in c[7]:   # CLNSIG ~ Pathogenic
                continue
            pts.setdefault(c[0], []).append(int(c[1]))
    for chrom in pts:
        pts[chrom].sort()
    return pts


def _near(pos, sorted_pos, window):
    """True iff some target lies within `window` bp of pos."""
    if not sorted_pos:
        return False
    i = bisect.bisect_left(sorted_pos, pos)
    for j in (i - 1, i):
        if 0 <= j < len(sorted_pos) and abs(sorted_pos[j] - pos) <= window:
            return True
    return False


def count_overlaps(deltas, targets, window):
    return sum(_near(pos, targets.get(chrom, []), window) for chrom, pos, _ in deltas)


def permutation_test(deltas, targets, sizes, window, perms, seed=0):
    """observed vs. the null of uniformly re-placed deltas (same count per
    chromosome).  Uses numpy when available (vectorised searchsorted, ~100x
    faster on 10^5-10^6 deltas); the pure-stdlib loop below is the reference
    implementation and gives the same statistic."""
    observed = count_overlaps(deltas, targets, window)
    # group delta positions by chromosome for fast reshuffling
    by_chrom = {}
    for chrom, pos, _ in deltas:
        by_chrom.setdefault(chrom, 0)
        by_chrom[chrom] += 1
    null = []
    try:
        import numpy as np
    except ImportError:  # pragma: no cover - stdlib fallback
        np = None
    if np is not None:
        rng = np.random.default_rng(seed)
        prepared = []
        for chrom, n in by_chrom.items():
            tgt = targets.get(chrom, [])
            L = sizes.get(chrom, 0)
            if tgt and L > 0:
                prepared.append((np.asarray(tgt, dtype=np.int64), int(L), int(n)))
        for _ in range(perms):
            cnt = 0
            for T, L, n in prepared:
                pos = rng.integers(1, L + 1, size=n, dtype=np.int64)
                i = np.searchsorted(T, pos)                        # == bisect_left
                hit = np.zeros(n, dtype=bool)
                ok = i < len(T)
                hit[ok] = np.abs(T[i[ok]] - pos[ok]) <= window
                ok = i > 0
                hit[ok] |= np.abs(T[i[ok] - 1] - pos[ok]) <= window
                cnt += int(hit.sum())
            null.append(cnt)
    else:
        rng = random.Random(seed)
        for _ in range(perms):
            cnt = 0
            for chrom, n in by_chrom.items():
                tgt = targets.get(chrom, [])
                L = sizes.get(chrom, 0)
                if not tgt or L <= 0:
                    continue
                for _ in range(n):
                    if _near(rng.randint(1, L), tgt, window):
                        cnt += 1
            null.append(cnt)
    null_mean = sum(null) / perms if perms else float("nan")
    fold = observed / null_mean if null_mean > 0 else float("inf")
    pval = (sum(1 for x in null if x >= observed) + 1) / (perms + 1)
    return observed, null_mean, fold, pval


def analyse(delta, gnomad, sizes, args, label, out_rows):
    obs, exp, fold, p = permutation_test(delta, gnomad, sizes, args.window, args.perms, args.seed)
    print(f"  {label:<34} obs={obs:<6d} exp={exp:7.1f}  fold={fold:5.2f}x  p={p:.4g}")
    out_rows.append((label, obs, round(exp, 2), round(fold, 3), p))


# -- synthetic self-test -------------------------------------------------------
def self_test():
    print("[self-test] planted-enrichment + negative control")
    rng = random.Random(1)
    sizes = {f"chr{i}": 5_000_000 for i in range(1, 6)}
    # high-F_ST "SV" sites
    targets = {c: sorted(rng.randint(1, L) for _ in range(300)) for c, L in sizes.items()}
    window = 1000
    # positive: deltas sit ON target sites (+- small jitter) -> should enrich
    pos_delta = []
    for c, ts in targets.items():
        for t in rng.sample(ts, 120):
            pos_delta.append((c, t + rng.randint(-200, 200), "A"))
    # negative: deltas placed uniformly at random -> should NOT enrich
    neg_delta = [(c, rng.randint(1, sizes[c]), "A")
                 for c in sizes for _ in range(120)]
    rows = []
    analyse(pos_delta, targets, sizes, argparse.Namespace(window=window, perms=1000, seed=2),
            "planted (expect fold>>1, p<0.01)", rows)
    analyse(neg_delta, targets, sizes, argparse.Namespace(window=window, perms=1000, seed=3),
            "random control (expect fold~1)", rows)
    ok = rows[0][3] > 3 and rows[0][4] < 0.01 and 0.5 < rows[1][3] < 2.0 and rows[1][4] > 0.05
    print("[self-test]", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--delta"); ap.add_argument("--gnomad-sv"); ap.add_argument("--sizes")
    ap.add_argument("--clinvar", default=None)
    ap.add_argument("--afr-field", default="AFR_AF"); ap.add_argument("--eas-field", default="EAS_AF")
    ap.add_argument("--fst-threshold", type=float, default=0.2)
    ap.add_argument("--fst-sweep", default=None,
                    help="comma-separated F_ST thresholds (e.g. 0.05,0.1,0.2,0.3,0.5): "
                         "report enrichment at each to show a dose-response, reusing "
                         "the same delta (no graph rebuild). Writes a sweep TSV to --out.")
    ap.add_argument("--min-len", type=int, default=0,
                    help="keep only delta alleles with INFO LEN >= this (e.g. 50 for "
                         "SV-sized deltas - the like-for-like set vs structural variants)")
    ap.add_argument("--window", type=int, default=1000)
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None, help="write a results TSV (for fig_biology)")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        sys.exit(self_test())
    if not (a.delta and a.gnomad_sv and a.sizes):
        ap.error("--delta, --gnomad-sv and --sizes are required (or use --self-test)")

    sizes = read_sizes(a.sizes)
    deltas = read_delta(a.delta, a.min_len)
    if a.min_len:
        print(f"[filter] kept {len(deltas)} delta alleles with LEN>={a.min_len}")
    if not deltas:
        print("no delta alleles to test (check --min-len / the delta VCF)"); return

    # -- F_ST dose-response sweep -----------------------------------------------
    if a.fst_sweep:
        thresholds = sorted(float(t) for t in a.fst_sweep.split(","))
        recs, n_tot, n_use = read_gnomad_records(a.gnomad_sv, a.afr_field, a.eas_field)
        print(f"delta alleles: {len(deltas)}   gnomAD-SV: {n_tot} ({n_use} with AF) "
              f" -  F_ST dose-response over {thresholds}")
        rows = []
        for t in thresholds:
            hi, n_hi = targets_at(recs, t)
            obs, exp, fold, p = permutation_test(deltas, hi, sizes, a.window, a.perms, a.seed)
            print(f"  F_ST≥{t:<5} n_high={n_hi:<6d} obs={obs:<6d} exp={exp:7.1f} "
                  f"fold={fold:5.2f}x  p={p:.4g}")
            rows.append((t, n_hi, obs, round(exp, 2), round(fold, 3), p))
        if a.out:
            with open(a.out, "w", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["fst_threshold", "n_high", "observed", "expected", "fold", "p"])
                w.writerows(rows)
            print(f"wrote {a.out}")
        return

    hi, n_tot, n_hi, n_use = read_gnomad_fst(a.gnomad_sv, a.afr_field, a.eas_field, a.fst_threshold)
    print(f"delta alleles: {len(deltas)}   gnomAD-SV: {n_tot} ({n_use} with AF, "
          f"{n_hi} high-F_ST≥{a.fst_threshold})")
    rows = []
    analyse(deltas, hi, sizes, a, f"delta vs high-F_ST gnomAD-SV (≥{a.fst_threshold})", rows)
    if a.clinvar:
        cv = read_points(a.clinvar, want_pathogenic=True)
        analyse(deltas, cv, sizes, a, "delta vs ClinVar pathogenic", rows)
    if a.out:
        with open(a.out, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["category", "observed", "expected", "fold", "p"])
            w.writerows(rows)
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
