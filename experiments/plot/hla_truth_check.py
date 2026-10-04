#!/usr/bin/env python3
"""
hla_truth_check.py — score a `g2g --ref` delta VCF against the variant lists
the two population graphs were BUILT from (demo/04_real_hla_demo.sh).

Ground truth: an atomic variant (0-based pos, ref, alt) is A-only if it is in
popA.vcf and not in popB.vcf (and symmetrically).  VCF records are atomised
the same way PanG2G does: the shared leading anchor base is stripped and
multi-allelic ALTs are split.

Usage:
    python3 experiments/plot/hla_truth_check.py popA.vcf popB.vcf delta.vcf
"""
import sys


def atomise(pos1, ref, alt):
    """VCF (1-based POS, REF, ALT) -> (0-based pos, ref, alt) without the
    common prefix/suffix, matching PanG2G's atomic-variant keys."""
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    ref2, alt2, pos0 = ref[i:], alt[i:], pos1 - 1 + i
    j = 0
    while j < len(ref2) and j < len(alt2) and ref2[-1 - j] == alt2[-1 - j]:
        j += 1
    if j:
        ref2, alt2 = ref2[:-j], alt2[:-j]
    return (pos0, ref2, alt2)


def read_vcf_variants(path):
    out = set()
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        for alt in c[4].split(","):
            if alt == "." or alt.startswith("<"):
                continue
            out.add(atomise(int(c[1]), c[3], alt))
    return out


def read_delta(path):
    a, b = set(), set()
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 8:
            continue
        ref = "" if c[3] == "." else c[3]
        alt = "" if c[4] == "." else c[4]
        key = (int(c[1]), ref, alt)
        (a if "GRAPH=A" in c[7] else b).add(key)
    return a, b


def pr(rep, tru):
    tp = len(rep & tru); fp = len(rep - tru); fn = len(tru - rep)
    p = tp / (tp + fp) if tp + fp else 1.0
    r = tp / (tp + fn) if tp + fn else 1.0
    return tp, fp, fn, p, r


def main():
    va, vb, delta = sys.argv[1:4]
    A, B = read_vcf_variants(va), read_vcf_variants(vb)
    tru_a, tru_b = A - B, B - A
    rep_a, rep_b = read_delta(delta)
    print(f"truth: |A|={len(A)} |B|={len(B)} shared={len(A & B)} A-only={len(tru_a)} B-only={len(tru_b)}")
    for lab, rep, tru in (("A-only", rep_a, tru_a), ("B-only", rep_b, tru_b)):
        tp, fp, fn, p, r = pr(rep, tru)
        print(f"{lab}: reported={len(rep)} tp={tp} fp={fp} fn={fn} precision={p:.4f} recall={r:.4f}")
        if fp:
            ex = sorted(rep - tru)[:5]
            print(f"   e.g. false positives: {ex}")
        if fn:
            ex = sorted(tru - rep)[:5]
            print(f"   e.g. missed: {ex}")


if __name__ == "__main__":
    main()
