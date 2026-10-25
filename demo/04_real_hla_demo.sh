#!/usr/bin/env bash
# Demo 2 - REAL biological data: HLA-DRB1 pangenomes.
#
# Builds two population pangenome graphs from real HLA-DRB1 haplotypes the way
# real ones are made (align to GRCh38 -> call variants -> vg construct), then
# diffs them with `g2g --ref GRCh38`. HLA-DRB1 is a standard, highly polymorphic
# pangenome benchmark.
#
# Requires the tools from 01_setup_tools.sh (vg, minimap2, paftools, pysam).
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
need_g2g
for t in vg minimap2 k8; do
    have "$t" || { echo "ERROR: '$t' not found - run demo/01_setup_tools.sh" >&2; exit 1; }
done
PAFTOOLS="$TOOLS/paftools.js"
OUT="$WORK/hla"; mkdir -p "$OUT"

log "1. Download real HLA-DRB1 haplotypes (12 GRCh38 sequences, pggb test data)"
if [ ! -s "$OUT/DRB1.fa" ]; then
    curl -fsSL --retry 3 -o "$OUT/DRB1.fa.gz" \
        "https://raw.githubusercontent.com/pangenome/pggb/master/data/HLA/DRB1-3123.fa.gz"
    gunzip -kf "$OUT/DRB1.fa.gz"
fi
note "$(grep -c '^>' "$OUT/DRB1.fa") haplotypes"

log "2. Designate GRCh38 (haplotype 0) as the shared reference; split the rest"
python3 - "$OUT" <<'PY'
import sys
out = sys.argv[1]
recs, name, seq = [], None, []
for L in open(f"{out}/DRB1.fa"):
    if L.startswith(">"):
        if name: recs.append((name, "".join(seq)))
        name, seq = L[1:].split()[0], []
    else: seq.append(L.strip())
if name: recs.append((name, "".join(seq)))
open(f"{out}/ref.fa", "w").write(f">GRCh38\n{recs[0][1]}\n")
for i, (_, s) in enumerate(recs[1:]):
    open(f"{out}/hap{i}.fa", "w").write(f">hap{i}\n{s}\n")
print(f"    reference GRCh38 ({len(recs[0][1])} bp) + {len(recs)-1} haplotypes")
PY
python3 -c "import pysam; pysam.faidx('$OUT/ref.fa')"

log "3. Call variants per haplotype (minimap2 align -> paftools call)"
for h in 0 1 2 3 4 5; do
    minimap2 -cx asm20 --cs "$OUT/ref.fa" "$OUT/hap$h.fa" 2>/dev/null \
        | sort -k6,6 -k8,8n > "$OUT/hap$h.paf"
    k8 "$PAFTOOLS" call -f "$OUT/ref.fa" -s "hap$h" -L 50 "$OUT/hap$h.paf" \
        2>/dev/null > "$OUT/hap$h.vcf"
    note "hap$h: $(grep -vc '^#' "$OUT/hap$h.vcf") variants"
done

log "4. Merge into two population VCFs (A = haps 0,2,4 ; B = haps 1,3,5)"
python3 - "$OUT" <<'PY'
import sys, collections, pysam
out = sys.argv[1]
def merge(haps, dst):
    sites = collections.defaultdict(set)
    for h in haps:
        for L in open(f"{out}/hap{h}.vcf"):
            if L.startswith("#"): continue
            f = L.rstrip("\n").split("\t")
            if len(f) < 5 or "," in f[3] or "," in f[4]: continue
            sites[(int(f[1]), f[3])].add(f[4])
    with open(f"{out}/{dst}.vcf", "w") as o:
        o.write("##fileformat=VCFv4.2\n##contig=<ID=GRCh38>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for pos, ref in sorted(sites):
            o.write(f"GRCh38\t{pos}\t.\t{ref}\t{','.join(sorted(sites[(pos,ref)]))}\t60\tPASS\t.\n")
    pysam.tabix_compress(f"{out}/{dst}.vcf", f"{out}/{dst}.vcf.gz", force=True)
    pysam.tabix_index(f"{out}/{dst}.vcf.gz", preset="vcf", force=True)
    print(f"    {dst}: {len(sites)} variant sites")
merge([0,2,4], "popA"); merge([1,3,5], "popB")
PY

log "5. Build acyclic, base-level pangenome graphs (vg construct)"
for g in A B; do
    vg construct -r "$OUT/ref.fa" -v "$OUT/pop$g.vcf.gz" 2>/dev/null > "$OUT/cg_$g.vg"
    vg convert -f "$OUT/cg_$g.vg" 2>/dev/null > "$OUT/cg_$g.gfa"
    note "cg_$g: $(vg stats -z "$OUT/cg_$g.vg" 2>/dev/null | tr '\n' ' ') $(vg stats -A "$OUT/cg_$g.vg" 2>/dev/null | head -1)"
done

log "6. Diff:  g2g cg_A.gfa cg_B.gfa --ref GRCh38 --vcf hla_diff.vcf"
"$G2G" "$OUT/cg_A.gfa" "$OUT/cg_B.gfa" --ref GRCh38 \
    --gaf "$OUT/hla_diff.gaf" --vcf "$OUT/hla_diff.vcf"

log "7. Sample of the variant-level diff (real chr coordinates)"
grep -v '^##' "$OUT/hla_diff.vcf" | sed -n '1,8p'   # not head: SIGPIPE under pipefail

log "8. Self-diff (cg_A vs cg_A) - must be zero differences"
"$G2G" "$OUT/cg_A.gfa" "$OUT/cg_A.gfa" --ref GRCh38 2>&1 \
    | grep -o "a_only=[0-9]*  b_only=[0-9]*" | sed 's/^/  /'

echo; echo "OK - real HLA demo complete (outputs in $OUT)"
