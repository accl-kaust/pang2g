#!/usr/bin/env python3
"""
pairing_error_sweep.py - Exp1b: snarl-pairing error as a function of graph size.

For each configuration (n_snarls, boundary-node length, private-snarl fraction,
anchor-linker spacing) and seed, this script

  1. generates a planted graph pair with gen_synthetic_gfa.py (planted mode),
  2. runs the snarl-pairing aligner (`g2g A B`, NOT the --ref diff),
  3. scores the reported delta against the planted truth, for BOTH graphs,
     by exact allele sequence (multiset match),
  4. records precision / recall / F1, the number of snarl pairs and wall time.

NOTE: the existing experiments/exp1_correctness/run.sh sweep parses the wrong
VCF column (it reads INFO 'GRAPH=A' as the allele sequence), so its recall
numbers are meaningless; use this script instead.

Usage (from the repo root, after building build/release/g2g):
    python3 experiments/exp1_correctness/pairing_error_sweep.py \
        --g2g build/release/g2g --out results/exp1b_pairing_error.csv

    # a quick smoke run
    python3 experiments/exp1_correctness/pairing_error_sweep.py --quick

Knobs (all optional): --sizes, --anchor-lens, --private-fracs, --linker-every,
--seeds, --pairing {default,boundary,interval,both}, --threads.
The --anchor-len / --private-frac / --linker-every generator options are
planned extensions of gen_synthetic_gfa.py (planted mode); until they land,
run with --anchor-lens 8 --private-fracs 0 --linker-every 0 and the generator
flags are simply not passed.
"""
import argparse, csv, os, re, subprocess, sys, tempfile, time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
GEN = os.path.join(ROOT, "experiments", "exp3_throughput", "gen_synthetic_gfa.py")


def load_truth(path):
    with open(path) as f:
        return Counter(l.strip() for l in f if l.strip())


def load_delta(vcf):
    """-> (Counter of A-only ALT sequences, Counter of B-only ALT sequences)."""
    a, b = Counter(), Counter()
    with open(vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            alt = p[4]
            if alt == ".":
                alt = ""
            if "GRAPH=A" in p[7]:
                a[alt] += 1
            elif "GRAPH=B" in p[7]:
                b[alt] += 1
    return a, b


def pr(reported: Counter, truth: Counter):
    tp = sum((reported & truth).values())
    fp = sum(reported.values()) - tp
    fn = sum(truth.values()) - tp
    prec = tp / (tp + fp) if tp + fp else 1.0
    rec = tp / (tp + fn) if tp + fn else 1.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return tp, fp, fn, prec, rec, f1


def run_one(g2g, n, k_shared, k_only, anchor_len, private_frac, linker_every,
            seed, pairing, threads, workdir):
    ga = os.path.join(workdir, "ga.gfa"); gb = os.path.join(workdir, "gb.gfa")
    ta = os.path.join(workdir, "truth_a.txt"); tb = os.path.join(workdir, "truth_b.txt")
    vcf = os.path.join(workdir, "out.vcf")
    if os.path.exists(vcf):
        os.remove(vcf)
    gen_cmd = [sys.executable, GEN, "--mode", "planted", "--n-snarls", str(n),
               "--k-shared", str(k_shared), "--k-only", str(k_only),
               "--seed", str(seed),
               "--out-a", ga, "--out-b", gb, "--truth-a", ta, "--truth-b", tb]
    # Only pass the pairing-stress knobs when they differ from the generator
    # defaults, so the script also works with the original generator.
    if anchor_len != 8:
        gen_cmd += ["--anchor-len", str(anchor_len)]
    if private_frac > 0:
        gen_cmd += ["--private-frac", str(private_frac)]
    if linker_every > 0:
        gen_cmd += ["--linker-every", str(linker_every)]
    gen = subprocess.run(gen_cmd, capture_output=True, text=True)
    if gen.returncode != 0:
        sys.exit(f"generator failed:\n{gen.stderr}")
    m = re.search(r"slots=(\{.*\})", gen.stderr)
    slots = eval(m.group(1)) if m else {}

    cmd = [g2g, ga, gb, "--threads", str(threads), "--gaf", os.devnull,
           "--vcf", vcf, "--verbose-timing"]
    if pairing:
        cmd += ["--pairing", pairing]
    t0 = time.perf_counter()
    run = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if run.returncode not in (0, 2):
        sys.exit(f"g2g failed (rc={run.returncode}):\n{run.stderr}")
    err = run.stderr
    pairs = int(re.search(r"snarl_pairs=(\d+)", err).group(1)) if "snarl_pairs=" in err else 0
    tw = re.search(r" wall=([0-9.]+)", err)
    g2g_wall = float(tw.group(1)) if tw else wall

    rep_a, rep_b = load_delta(vcf) if os.path.exists(vcf) else (Counter(), Counter())
    tru_a, tru_b = load_truth(ta), load_truth(tb)
    ra = pr(rep_a, tru_a); rb = pr(rep_b, tru_b)
    return {
        "n_snarls": n, "k_shared": k_shared, "k_only": k_only,
        "anchor_len": anchor_len, "private_frac": private_frac,
        "linker_every": linker_every, "seed": seed, "pairing": pairing or "default",
        "threads": threads,
        "shared_slots": slots.get("shared_slots", ""),
        "a_private_slots": slots.get("a_private_slots", ""),
        "b_private_slots": slots.get("b_private_slots", ""),
        "linkers": slots.get("linkers", ""),
        "snarl_pairs": pairs,
        "tp_a": ra[0], "fp_a": ra[1], "fn_a": ra[2],
        "precision_a": f"{ra[3]:.6f}", "recall_a": f"{ra[4]:.6f}", "f1_a": f"{ra[5]:.6f}",
        "tp_b": rb[0], "fp_b": rb[1], "fn_b": rb[2],
        "precision_b": f"{rb[3]:.6f}", "recall_b": f"{rb[4]:.6f}", "f1_b": f"{rb[5]:.6f}",
        "wall_s": f"{g2g_wall:.4f}",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--g2g", default=os.path.join(ROOT, "build", "release", "g2g"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "exp1b_pairing_error.csv"))
    ap.add_argument("--sizes", default="5,10,20,50,100,1000,10000,100000")
    ap.add_argument("--anchor-lens", default="8")
    ap.add_argument("--private-fracs", default="0")
    ap.add_argument("--linker-every", default="0")
    ap.add_argument("--seeds", default="42,123,999")
    ap.add_argument("--k-shared", type=int, default=2)
    ap.add_argument("--k-only", type=int, default=1)
    ap.add_argument("--pairing", default="default",
                    help="'default' (no flag), 'boundary', 'interval' or 'both'")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--quick", action="store_true", help="tiny smoke sweep")
    a = ap.parse_args()

    if a.quick:
        a.sizes, a.seeds = "5,50,1000", "42"

    sizes = [int(x) for x in a.sizes.split(",")]
    alens = [int(x) for x in a.anchor_lens.split(",")]
    pfs = [float(x) for x in a.private_fracs.split(",")]
    links = [int(x) for x in a.linker_every.split(",")]
    seeds = [int(x) for x in a.seeds.split(",")]
    pairings = {"default": [None], "boundary": ["boundary"], "interval": ["interval"],
                "both": ["boundary", "interval"]}[a.pairing]

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    rows = []
    with tempfile.TemporaryDirectory(prefix="g2g_exp1b_") as wd:
        for pairing in pairings:
            for n in sizes:
                for al in alens:
                    for pf in pfs:
                        for le in links:
                            for sd in seeds:
                                r = run_one(a.g2g, n, a.k_shared, a.k_only, al, pf, le,
                                            sd, pairing, a.threads, wd)
                                rows.append(r)
                                print(f"[exp1b] pairing={r['pairing']:<8} n={n:<7} L={al} "
                                      f"priv={pf:<4} link={le} seed={sd:<4} pairs={r['snarl_pairs']:<7} "
                                      f"P/R(A)={float(r['precision_a']):.3f}/{float(r['recall_a']):.3f} "
                                      f"P/R(B)={float(r['precision_b']):.3f}/{float(r['recall_b']):.3f} "
                                      f"wall={r['wall_s']}s", flush=True)
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"[exp1b] wrote {len(rows)} rows -> {a.out}")


if __name__ == "__main__":
    main()
