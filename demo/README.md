# G2G demo — diffing two pangenome graphs

This directory is a **self-contained, reproducible demonstration** of G2G's
graph-to-graph diff: give it two pangenome graphs that share a reference, and it
reports the variants that differ between them, with real reference coordinates.

```
g2g  graph_A.gfa  graph_B.gfa  --ref GRCh38  --vcf diff.vcf
```

It outputs a VCF-style diff — one record per allele present in only one of the
two graphs:

```
#CHROM  POS   ID    REF  ALT    ...  INFO
GRCh38  46    A_0   C    G      ...  GRAPH=A      <- unique to graph A
GRCh38  1516  A_2   GAA  .      ...  GRAPH=A      <- a deletion, unique to A
GRCh38  249   B_0   G    C      ...  GRAPH=B      <- unique to graph B
```

Variants are matched between the two graphs by **reference coordinate**, after
decomposing each bubble into atomic SNPs/indels. This is robust to the two
graphs chopping the shared reference into different nodes — the situation that
breaks naive bubble-boundary matching on real graphs.

---

## Quick start

```bash
# 1. everything (build + synthetic + scaling sweep + fetch tools + real HLA demo)
bash demo/run_all.sh

# synthetic + scaling only (no downloads, no extra tools):
SKIP_REAL=1 bash demo/run_all.sh

# skip just the scaling sweep (it adds a few seconds):
SKIP_SCALING=1 bash demo/run_all.sh
```

Logs from a real run are committed under [`logs/`](logs/); example outputs under
[`sample_outputs/`](sample_outputs/).

## Prerequisites

| For | You need |
|-----|----------|
| Building G2G | a C++20 compiler (GCC ≥ 9 or Clang ≥ 10), `cmake`, `ninja` |
| Synthetic demo | the above + `python3` |
| Real HLA demo | the above + outbound HTTPS (the setup script fetches `vg`, `minimap2`, `paftools`, and `pysam`) |

The setup downloads pinned, static binaries into `demo/tools/` — nothing is
installed system-wide. `tools/` and `work/` are git-ignored.

> **Older clusters (GCC 9):** GCC 9's libstdc++ has no `<span>` header. G2G
> ships a header-only `g2g::span` shim (`include/g2g/span.hpp`) that transparently
> falls back to a built-in view when `std::span` is unavailable, so the stock
> system compiler on e.g. Ubuntu 20.04 (GCC 9.4) builds the demo as-is — no
> `module load` or conda compiler required.

---

## Demo 1 — controlled synthetic pair (known ground truth)

`03_synthetic_demo.sh` → [`logs/03_synthetic_demo.log`](logs/03_synthetic_demo.log)

[`make_synthetic_pair.py`](make_synthetic_pair.py) builds two graphs that embed
the **same reference path** but (a) plant a known set of variants and (b) chop
the reference into **different node sizes** (A = 25 bp, B = 17 bp). `g2g --ref`
must recover exactly the planted variants regardless of that segmentation.

```
planted : shared=5 A_only=3 B_only=2
reported: A_only=3 B_only=2
RESULT  : PASS — diff matches ground truth exactly
```

It also runs a **self-diff** (a graph against itself), which must report
**zero** differences — a basic correctness invariant.

This demo needs no network and no third-party tools — just `g2g` + `python3`.

## Demo 2 — real biological data (HLA-DRB1)

`04_real_hla_demo.sh` → [`logs/04_real_hla_demo.log`](logs/04_real_hla_demo.log)

Builds two population pangenomes from **real HLA-DRB1 haplotypes** (12 GRCh38
sequences, a standard, highly polymorphic pangenome benchmark) the way real ones
are made, then diffs them:

1. download the haplotypes (pggb test data);
2. designate GRCh38 (haplotype 0) as the shared reference;
3. align the other haplotypes to it (`minimap2`) and call variants (`paftools`);
4. split into two populations (A = haps 0,2,4; B = haps 1,3,5) and build a
   base-level, **acyclic** graph for each (`vg construct`);
5. diff with `g2g --ref GRCh38`.

```
cg_A: 4500 nodes (acyclic)   cg_B: 3742 nodes (acyclic)
g2g --ref GRCh38  ->  snarl_pairs=863  a_only=266  b_only=93   (~0.1 s)
self-diff (cg_A vs cg_A)  ->  a_only=0  b_only=0
```

The output is a variant-level VCF with real `GRCh38:pos` coordinates, including
indels (see [`sample_outputs/hla_diff_sample.vcf`](sample_outputs/hla_diff_sample.vcf)).

## Demo 3 — scaling beyond the gene region

`05_scaling_demo.sh` → [`logs/05_scaling_demo.log`](logs/05_scaling_demo.log)

The same controlled diff, driven across graphs growing from ~100 kb to **8 Mb**
(~700× the HLA gene region) at a dense ~1 variant / 200 bp. At every size g2g
must still recover the planted A-only / B-only set **exactly** — so this measures
speed and memory *without giving up correctness*. Wall time and peak RSS are
g2g's own (`--verbose-timing` → `[timing] wall=`, `[mem] peak_rss_mb=`), so no
external profiler is needed.

| ref bp | nodes (A+B) | snarls | a_only / b_only | exact? | wall | peak RSS |
|-------:|------------:|-------:|----------------:|:------:|-----:|---------:|
| 100 kb |        12 k |    400 | 50 / 50         |   ✓    | 0.02 s |   8 MB |
| 500 kb |        59 k |  2 000 | 250 / 250       |   ✓    | 0.10 s |  20 MB |
| 2 Mb   |       237 k |  8 000 | 1 000 / 1 000   |   ✓    | 0.49 s |  68 MB |
| 8 Mb   |       948 k | 32 000 | 4 000 / 4 000   |   ✓    | 2.27 s | 258 MB |

Time and memory both scale ~linearly (each 4× in size costs ~4.6× time):
~18 k variants/s, ~420 k nodes/s on a single core. Push further with
`SCALE_MAX_MB=32` or a custom `SCALE_SIZES="bp:variants ..."`.

> The built-in snarl finder used to be quadratic here (a whole-graph reachability
> search per bubble). It now uses an exact, topologically-windowed sink search —
> falling back to the full search only for bubbles larger than the window — which
> cut the 4 Mb point from ~52 s to ~1 s with byte-identical boundaries.

---

## Using G2G on your own graphs

```bash
g2g A.gfa B.gfa --ref <name> --vcf diff.vcf       # variant-level diff
```

- Both graphs must embed the reference as a path. GFA **P-lines and W-lines**
  are both supported; resolve a W-line reference by its sample or contig, e.g.
  `--ref GRCh38` or `--ref chr1` (no need for the full `GRCh38#0#chr1:...` name).
- A reference spanning several contigs is diffed per-contig.
- Other useful flags: `--snarls vg` (use the vg toolkit's snarl finder instead
  of the built-in one; needs `vg` on `PATH` or `$G2G_VG`), `--verbose-timing`.

## How the diff works (one paragraph)

Each graph is decomposed into snarls (bubbles). For every bubble anchored on the
reference path, each allele is aligned to the bubble's reference allele and the
constituent SNPs/indels are emitted, keyed by their **absolute reference
position** (`bubble_open_offset + offset_within_allele`). Because that position
depends only on the reference coordinate — not on node boundaries or on how many
adjacent variants a bubble happens to span — a merged bubble in one graph yields
the same atomic variants as the separate bubbles in the other, and corresponding
variants match. Variants present in only one graph become that graph's delta.

## Validation evidence

- **Self-consistency:** a graph diffed against itself yields 0 differences
  (both demos).
- **Exact on known truth:** the synthetic pair recovers the planted
  shared / A-only / B-only set exactly, despite different node segmentation.
- **Merged vs split bubbles:** unit test `ReferenceDiff.MergedBubbleMatchesSplitBubbles`
  (in `tests/test_refdiff.cpp`) proves a single merged bubble matches two
  separate bubbles encoding the same SNPs.
- **Real data:** runs end-to-end on real HLA-DRB1 graphs in ~0.1 s and emits a
  coordinate-stamped variant diff.
- **Scales near-linearly:** the diff stays exact from 100 kb to 8 Mb (~950k
  nodes, ~32k bubbles) with ~linear time and memory (Demo 3).

## Limitations (honest)

- **DAG-only.** Decomposition requires an acyclic graph; cyclic input errors
  cleanly (`Graph contains cycles — DAG required`). Reference-anchored graphs
  from minigraph-cactus / `vg construct` are acyclic; some other builders
  (e.g. `vg msga`) can produce cycles.
- **Indels are not left-normalised**, so indel calls in repeats may differ by a
  base from a tool like `bcftools norm`.
- Validated to multi-megabase scale (8 Mb synthetic, ~950k nodes, ~2 s
  single-core; Demo 3). Whole-chromosome (50–250 Mb) and many-haplotype
  HPRC-scale graphs are the next step.
- The diff itself is single-threaded; the parallel work-stealing runtime
  accelerates the snarl-correspondence aligner, not the `--ref` coordinate diff.

## Files

| File | Purpose |
|------|---------|
| `run_all.sh` | run the whole demo, teeing logs to `logs/` |
| `01_setup_tools.sh` | fetch `vg`, `minimap2`, `paftools`, `pysam` into `tools/` |
| `02_build_g2g.sh` | build the `g2g` binary |
| `03_synthetic_demo.sh` | controlled pair vs known ground truth + self-diff |
| `04_real_hla_demo.sh` | real HLA-DRB1 pangenome diff |
| `05_scaling_demo.sh` | scaling sweep (100 kb → 8 Mb) with timing + peak memory |
| `make_synthetic_pair.py` | generator for the controlled / scaling pairs |
| `lib.sh` | shared paths/helpers |
| `logs/` | captured logs from a real run |
| `sample_outputs/` | example diff/truth files |
