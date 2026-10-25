# The demos

This directory is a self-contained demonstration of G2G's graph-to-graph
diff. Give it two pangenome graphs that share a reference and it reports the
variants that differ between them, with reference coordinates:

```
g2g  graph_A.gfa  graph_B.gfa  --ref GRCh38  --vcf diff.vcf
```

The output is a VCF-style diff with one record per allele present in only
one of the two graphs. POS is the 0-based reference offset of the site, so
add 1 to compare with a standard 1-based VCF:

```
#CHROM  POS   ID    REF  ALT    ...  INFO
GRCh38  46    A_0   C    G      ...  GRAPH=A      <- unique to graph A
GRCh38  1516  A_2   GAA  .      ...  GRAPH=A      <- a deletion, unique to A
GRCh38  249   B_0   G    C      ...  GRAPH=B      <- unique to graph B
```

Variants are matched between the two graphs by reference coordinate, after
each bubble has been decomposed into atomic SNPs and indels. Two graphs that
chop the shared reference into different nodes therefore still match, which
is exactly the situation that breaks naive bubble-boundary matching on real
graphs.

## Quick start

```bash
# the test case (build + planted ground truth + synthetic biological validation), from the repo root:
bash run_test_case.sh            # add --unit for the GoogleTests, --real for the HLA demo

# everything (build + synthetic + scaling sweep + fetch tools + real HLA demo)
bash demo/run_all.sh

# synthetic + scaling only (no downloads, no extra tools):
SKIP_REAL=1 bash demo/run_all.sh

# skip just the scaling sweep (it adds a few seconds):
SKIP_SCALING=1 bash demo/run_all.sh
```

Logs from a real run are committed under [`logs/`](logs/), and example
outputs under [`sample_outputs/`](sample_outputs/).

## Prerequisites

Building G2G needs a C++20 compiler (GCC 9 or newer, Clang 10 or newer) and
`cmake`; `ninja` is used if present. The synthetic demo needs `python3` on
top of that. The real HLA demo also needs outbound HTTPS, because the setup
script fetches `vg`, `minimap2`, `paftools` and `pysam`. The setup downloads
pinned static binaries into `demo/tools/`; nothing is installed system-wide,
and `tools/` and `work/` are git-ignored.

GCC 9's libstdc++ has no
`<span>` header. G2G carries a small header-only `g2g::span` shim
(`include/g2g/span.hpp`) that falls back to a built-in view when `std::span`
is unavailable, so the stock system compiler on, say, Ubuntu 20.04 builds the
demo as it is.

## Demo 1: a controlled synthetic pair with known ground truth

`03_synthetic_demo.sh`, log in [`logs/03_synthetic_demo.log`](logs/03_synthetic_demo.log).

[`make_synthetic_pair.py`](make_synthetic_pair.py) builds two graphs that
embed the same reference path but plant a known set of variants and chop the
reference into different node sizes (25 bp in A, 17 bp in B). `g2g --ref`
has to recover exactly the planted variants regardless of that segmentation:

```
planted : shared=5 A_only=3 B_only=2
reported: A_only=3 B_only=2
RESULT  : PASS - diff matches ground truth exactly
```

It also runs a self-diff, a graph against itself, which has to report zero
differences. This demo needs no network and no third-party tools, only
`g2g` and `python3`.

## Demo 2: real biological data (HLA-DRB1)

`04_real_hla_demo.sh`, log in [`logs/04_real_hla_demo.log`](logs/04_real_hla_demo.log).

This one builds two population pangenomes from real HLA-DRB1 haplotypes (12
GRCh38 sequences, a standard and highly polymorphic pangenome benchmark) the
way real pangenomes are built, then diffs them:

1. download the haplotypes (the pggb test data);
2. designate GRCh38 (haplotype 0) as the shared reference;
3. align the other haplotypes to it with `minimap2` and call variants with `paftools`;
4. split them into two populations (A = haplotypes 0, 2, 4; B = haplotypes 1, 3, 5)
   and build a base-level, acyclic graph for each with `vg construct`;
5. diff with `g2g --ref GRCh38`.

```
cg_A: 4500 nodes (acyclic)   cg_B: 3742 nodes (acyclic)
g2g --ref GRCh38  ->  snarl_pairs=1024  a_only=349  b_only=12   (about 0.1 s)
self-diff (cg_A vs cg_A)  ->  a_only=0  b_only=0
```

The output is a variant-level VCF with real `GRCh38:pos` coordinates,
indels included; see
[`sample_outputs/hla_diff_sample.vcf`](sample_outputs/hla_diff_sample.vcf).

## Demo 3: scaling beyond the gene region

`05_scaling_demo.sh`, log in [`logs/05_scaling_demo.log`](logs/05_scaling_demo.log).

The same controlled diff, driven across graphs growing from about 100 kb to
8 Mb (roughly 700 times the HLA gene region) at a dense one variant per
200 bp. At every size g2g still has to recover the planted A-only and B-only
sets exactly, so this measures speed and memory without giving up
correctness. Wall time and peak RSS come from g2g itself (`--verbose-timing`
prints `[timing] wall=` and `[mem] peak_rss_mb=`), so no external profiler is
needed.

| ref bp | nodes (A+B) | snarls | a_only / b_only | exact? | wall | peak RSS |
|-------:|------------:|-------:|----------------:|:------:|-----:|---------:|
| 100 kb |        12 k |    400 | 50 / 50         |   yes    | 0.02 s |   8 MB |
| 500 kb |        59 k |  2 000 | 250 / 250       |   yes    | 0.10 s |  20 MB |
| 2 Mb   |       237 k |  8 000 | 1 000 / 1 000   |   yes    | 0.49 s |  68 MB |
| 8 Mb   |       948 k | 32 000 | 4 000 / 4 000   |   yes    | 2.27 s | 258 MB |

Time and memory both scale close to linearly (each 4x in size costs about
4.6x in time): about 18 k variants/s and 420 k nodes/s on a single core.
`SCALE_MAX_MB=32` adds a larger point, and `SCALE_SIZES="bp:variants ..."`
replaces the sweep.

The built-in snarl finder used to be quadratic here (a
whole-graph reachability search per bubble). It now uses an exact,
topologically windowed sink search, and only falls back to the full search
for bubbles larger than the window, which cut the 4 Mb point from about 52 s
to about 1 s with byte-identical boundaries.

## Demo 3b: pushing to the machine's limit

`06_scaling_limit.sh`, log in [`logs/06_scaling_limit.log`](logs/06_scaling_limit.log).

Where `05` runs a fixed sweep, `06` keeps doubling the reference size until
the next point would exceed the node's RAM or a wall-clock budget, then stops
cleanly and reports the ceiling. It projects each point's peak memory from
the previous one and refuses to attempt a point it expects to run out of
memory on, so an unattended cluster run is never killed half way through and
the allocation is not wasted; the diff is checked for exactness at every
size.

On a 512 GB node it climbs to an 8 Gbp graph pair, past the whole human
genome, before hitting the RAM ceiling, exact at every step:

| ref | nodes (A+B) | variants | exact? | wall | peak RSS |
|----:|------------:|---------:|:------:|-----:|---------:|
| 256 Mb | 30.3 M | 1.28 M | yes |  104 s |   7.7 GiB |
|   1 Gb |  121 M  | 5.12 M | yes |  415 s |    31 GiB |
|   2 Gb |  243 M  | 10.2 M | yes |  817 s |    60 GiB |
|   4 Gb |  485 M  | 20.5 M | yes | 1880 s |   121 GiB |
|   8 Gb |  971 M  | 41.0 M | yes | 4292 s |   242 GiB |

That is about 971 million nodes and 41 million variants diffed exactly on a
single core in about 72 minutes, roughly 2.6 times the human genome
(3.1 Gb); the next doubling, 16 Gb, is the only thing the node's RAM cannot
hold. Across the whole sweep, from a million to nearly a billion nodes, time
and memory stay linear (about 226 k nodes/s on one core), and a smaller node
stops proportionally sooner (a 16 GB machine tops out around 256 Mb
and 30 M nodes).

```bash
# adapts to the node's RAM; tee a log to keep as evidence
bash demo/06_scaling_limit.sh 2>&1 | tee demo/logs/06_scaling_limit.log
# knobs: BUDGET_SEC (default 9 h), MEM_CEIL_MB (default 80% RAM), START_MB, MAX_MB
```

Fitting a power law to the `05` and `06` result TSVs, from 8 Mb to 8 Gb (11
points, a million to nearly a billion nodes), wall time scales as N^1.04 and
peak memory as N^0.99, both within a hair of linear, and snarl decomposition
is exactly N^1.00 (the old quadratic is gone).

## Demo 4: biological validation (population-differentiation enrichment)

`09_hprc_biology.sh`, with [`hprc_enrich.py`](hprc_enrich.py) as the engine.

This tests the biological claim. When G2G diffs two population
pangenomes (here AFR against EAS), the alleles it reports in one but not the
other, the delta graph, should sit closer than chance to structural variants
that are highly differentiated between those populations (high Hudson F_ST
in gnomAD-SV), and should touch known pathogenic (ClinVar) sites. The test
is a permutation null: the observed delta-to-target proximity against the
same count after shuffling every delta to a random position on its
chromosome.

```
fold = observed / mean(null)      p = (#{null >= observed} + 1) / (perms + 1)
```

The enrichment engine does not depend on any data and tests itself: it
plants a known enrichment, which it has to recover (fold well above 1,
p < 0.01), and a random negative control, which has to stay at a fold
around 1:

```bash
python3 demo/hprc_enrich.py --self-test     # stdlib only, no data, no tools
```

The synthetic end-to-end mode runs the whole pipeline with the real `g2g`
binary and needs no HPRC download and no `vg` or `bcftools`. It plants a
population difference, runs `g2g --ref`, synthesises a gnomAD-SV and ClinVar
fixture around the recovered deltas, and checks both that the signal
enriches and that a shuffled-delta control does not:

```bash
bash demo/09_hprc_biology.sh --synthetic
```
```
g2g --ref chr1  ->  a_only=140  b_only=160   (recovers all 300 planted deltas)
delta vs high-F_ST gnomAD-SV   fold=21.3x  p=0.001   [PASS]
delta vs ClinVar pathogenic    fold=42.0x  p=0.001   [PASS]
shuffled control (null)        fold= 1.3x  p=0.23    [PASS]   (specificity)
```

It writes its summary to `work/biology/biology_results.tsv`.

The real mode (`bash demo/09_hprc_biology.sh` without a flag) builds the AFR
and EAS pangenomes with `vg construct` from a population-labelled phased
panel, diffs them per chromosome with `g2g --ref`, and enriches against a
real gnomAD-SV sites VCF and ClinVar. It needs `vg`, `bcftools`, `tabix` and
the data files, supplied through environment variables (`REF_FA`,
`PANEL_VCF`, `SAMPLES_A`, `SAMPLES_B`, `GNOMAD_SV_VCF`, `CLINVAR_VCF`,
`CHROM_SIZES`, `CHROMS`); the header of the script lists them, and
`09_fetch_biology_data.sh` stages the public inputs. The analysis stage is
byte for byte the same engine as in the synthetic mode. A real chromosome-22
run (AFR against EAS, 1000 Genomes phase 3) is logged in
`logs/09_real_chr22.log`.

## Demos 5 and 6: parallel scaling and ablation

`07_parallel_scaling.sh` sweeps the thread count on one graph pair and
reports DP-phase against end-to-end speed-up, with a one-thread baseline; it
flags the run as invalid if the parallel result diverges from the sequential
one. `08_ablation.sh` builds one binary per design choice with that choice
disabled and compares each against the optimised baseline on a control and a
treatment input. Both write a TSV under `work/`.

## Using G2G on your own graphs

```bash
g2g A.gfa B.gfa --ref <name> --vcf diff.vcf       # variant-level diff
```

Both graphs must embed the reference as a path. GFA P-lines and W-lines are
both supported; a W-line reference resolves by its sample or contig, for
example `--ref GRCh38` or `--ref chr1`, without the full
`GRCh38#0#chr1:...` name. A reference spanning several contigs is diffed per
contig. Other useful flags are `--snarls vg`, which uses the vg toolkit's
snarl finder instead of the built-in one (needs `vg` on the PATH or in
`$G2G_VG`), and `--verbose-timing`.

## How the diff works

Each graph is decomposed into snarls (bubbles). For every bubble anchored on
the reference path, each allele is aligned to the bubble's reference allele
and the constituent SNPs and indels are emitted, keyed by their absolute
reference position (`bubble_open_offset + offset_within_allele`). Because
that position depends only on the reference coordinate, not on node
boundaries or on how many adjacent variants a bubble happens to span, a
merged bubble in one graph yields the same atomic variants as the separate
bubbles in the other, and corresponding variants match. Variants present in
only one graph become that graph's delta.

## What backs this up

A graph diffed against itself yields no differences (both demos check this).
The synthetic pair recovers the planted shared, A-only and B-only sets
exactly, despite the different node segmentation. The unit test
`ReferenceDiff.MergedBubbleMatchesSplitBubbles` in `tests/test_refdiff.cpp`
shows that a single merged bubble matches two separate bubbles encoding the
same SNPs. The real HLA-DRB1 graphs run end to end in about 0.1 s and give a
coordinate-stamped variant diff. And the diff stays exact from 100 kb to
8 Mb with close to linear time and memory (demo 3), and on to 8 Gbp, about
971 million nodes and 41 million variants, on a 512 GB node in about 72
minutes on one core (demo 3b).

## Limitations

Decomposition requires an acyclic graph; cyclic input stops cleanly with
`Graph contains cycles - DAG required`. Reference-anchored graphs from
Minigraph-Cactus and `vg construct` are acyclic; some other builders, `vg
msga` for example, can produce cycles. Indels are not left-normalised, so
indel calls in repeats may differ by a base from a tool like `bcftools norm`.
The scaling past whole-genome size was measured on synthetic graphs (demo
3b); the biological pipeline (demo 4) is validated end to end on synthetic
data with the real binary and on real chromosome-22 panels, while running it
on many-haplotype, HPRC-scale AFR and EAS pangenomes needs `vg`, `bcftools`
and the gnomAD-SV and ClinVar panels. The `--ref` diff itself is
single-threaded; the parallel work-stealing runtime accelerates the
snarl-correspondence aligner, not the coordinate diff.

## Files

| File | Purpose |
|------|---------|
| `run_all.sh` | run the whole demo, teeing logs to `logs/` |
| `01_setup_tools.sh` | fetch `vg`, `minimap2`, `paftools`, `pysam` into `tools/` |
| `02_build_g2g.sh` | build the `g2g` binary |
| `03_synthetic_demo.sh` | controlled pair vs known ground truth + self-diff |
| `04_real_hla_demo.sh` | real HLA-DRB1 pangenome diff |
| `05_scaling_demo.sh` | scaling sweep (100 kb to 8 Mb) with timing + peak memory |
| `06_scaling_limit.sh` | push scaling to the node's RAM / time limit (demo 3b) |
| `07_parallel_scaling.sh` | DP-phase vs end-to-end parallel speed-up sweep |
| `08_ablation.sh` | each design choice disabled vs the optimised baseline |
| `09_hprc_biology.sh` | biological validation: delta vs gnomAD-SV / ClinVar (synthetic + real) |
| `09_fetch_biology_data.sh` | stage GRCh37 inputs for the real run (1000G, gnomAD-SV, ClinVar) |
| `hprc_enrich.py` | permutation-null enrichment engine (`--self-test`, `--fst-sweep`, `--min-len`) |
| `10_throughput_baseline.sh` | G2G diff vs linearize-and-map, sweeping haplotype count |
| `make_synthetic_pair.py` | generator for the controlled / scaling pairs |
| `lib.sh` | shared paths and helpers |
| `logs/` | captured logs from a real run |
| `sample_outputs/` | example diff/truth files |
