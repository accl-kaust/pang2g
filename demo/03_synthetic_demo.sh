#!/usr/bin/env bash
# Demo 1 — controlled synthetic pair with KNOWN ground truth.
#
# Two graphs share a reference path but plant a known set of variants and chop
# the reference into DIFFERENT nodes. `g2g --ref` must recover exactly the
# planted shared / A-only / B-only variants — proving correctness and
# robustness to node segmentation. Also runs a self-diff (graph vs itself),
# which must report zero differences.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
need_g2g
OUT="$WORK/synthetic"; mkdir -p "$OUT"

log "1. Generate a controlled pair (5 shared, 3 A-only, 2 B-only; chop A=25 B=17)"
python3 "$DEMO_DIR/make_synthetic_pair.py" \
    --n-shared 5 --n-a 3 --n-b 2 --chop-a 25 --chop-b 17 --outdir "$OUT"
note "graph A: $(grep -c '^S' "$OUT/demoA.gfa") segments, $(grep -c '^P' "$OUT/demoA.gfa") path"
note "graph B: $(grep -c '^S' "$OUT/demoB.gfa") segments, $(grep -c '^P' "$OUT/demoB.gfa") path"

log "2. Diff:  g2g demoA.gfa demoB.gfa --ref GRCh38 --vcf diff.vcf"
"$G2G" "$OUT/demoA.gfa" "$OUT/demoB.gfa" --ref GRCh38 \
    --gaf "$OUT/diff.gaf" --vcf "$OUT/diff.vcf"

log "3. Diff output (variant-level, real coordinates)"
grep -v '^##' "$OUT/diff.vcf"

log "4. Check the diff against the planted ground truth"
python3 - "$OUT/truth.tsv" "$OUT/diff.vcf" <<'PY'
import sys
truth = {}
for ln in open(sys.argv[1]):
    if ln.startswith("pos"): continue
    pos, ref, alt, cls = ln.split()
    truth[(int(pos), ref, alt)] = cls
got = {"shared": set(), "A": set(), "B": set()}
for ln in open(sys.argv[2]):
    if ln.startswith("#"): continue
    f = ln.split("\t")
    key = (int(f[1]), f[3], f[4])
    got["A" if "GRAPH=A" in ln else "B"].add(key)
# shared = planted shared variants that are NOT reported as unique
exp = {"shared": set(), "A": set(), "B": set()}
for k, c in truth.items():
    exp[c].add(k)
ok = (got["A"] == exp["A"] and got["B"] == exp["B"])
print(f"  planted : shared={len(exp['shared'])} A_only={len(exp['A'])} B_only={len(exp['B'])}")
print(f"  reported: A_only={len(got['A'])} B_only={len(got['B'])}")
print("  RESULT  :", "PASS — diff matches ground truth exactly" if ok else "FAIL")
sys.exit(0 if ok else 1)
PY

log "5. Self-diff (graph vs itself) — must be zero differences"
"$G2G" "$OUT/demoA.gfa" "$OUT/demoA.gfa" --ref GRCh38 2>&1 \
    | grep -o "a_only=[0-9]*  b_only=[0-9]*" | sed 's/^/  /'

echo; echo "OK — synthetic demo complete (outputs in $OUT)"
