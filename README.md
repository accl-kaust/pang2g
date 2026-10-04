# PanG2G — graph-to-graph alignment of pangenome variation graphs

PanG2G aligns two pangenome variation graphs to each other at base-pair
resolution and reports the **delta graph**: every allele present in one graph
but absent from the other, with its exact sequence. It is the code behind

> A. Shakleya and S. A. Fahmy, *PanG2G: Scalable Base-Pair-Resolution
> Graph-to-Graph Alignment of Pangenome Graphs*, IEEE BIBM 2026.

Two modes:

| mode | command | pairs snarls by | best for |
|---|---|---|---|
| reference-anchored diff | `g2g A.gfa B.gfa --ref <path> --vcf diff.vcf` | absolute reference coordinate (exact under identical variant representation) | graphs that embed a common reference (vg construct, Minigraph-Cactus) |
| graph-to-graph aligner | `g2g A.gfa B.gfa [--pairing boundary\|interval] --vcf diff.vcf` | boundary-sequence key (default, linear) or anchor-interval dynamic program | graphs without a shared reference; parallel (OpenMP) |

Everything reported in the paper (numbers, tables, figures) can be regenerated
from this repository; see [Reproducing the paper](#reproducing-the-paper).

## Layout

```
include/g2g/, src/     C++20 engine (library + `g2g` CLI)
tests/                 GoogleTest suite (fetched by CMake; -DG2G_BUILD_TESTS=OFF to skip)
tools/                 bench_matching: matching cost vs. allele count (fig_matching)
python/g2g/            Python reference prototype (correctness oracle; quadratic greedy pairing)
python/                profile_workload.py, synth_generator.py (workload characterization)
demo/                  end-to-end demos + harnesses: synthetic truth, real HLA-DRB1 graphs,
                       scaling sweeps, parallel scaling, ablations, Q4 biology, throughput
experiments/           Slurm-ready experiment suite (exp1-exp6), generators and plot scripts
results/               measured result tables behind the figures
demo/sample_outputs/   result tables/figures of the real-data runs
make_figures.py        regenerates every figure in latex/figs/ from results/ + demo/sample_outputs/
latex/                 the paper
```

## Build

Requirements: a C++20 compiler (GCC >= 9 or Clang >= 10), CMake >= 3.18,
OpenMP (ships with GCC), Python 3 with numpy/matplotlib for figures.

```bash
cmake -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j"$(nproc)"
cd build/release && ctest            # 161 tests; vg/minimap2 tests skip if the tools are absent
```

`demo/02_build_g2g.sh` does the same without the test fetch.

## Run

```bash
# reference-anchored diff (both graphs must embed the path; P- or W-lines)
build/release/g2g A.gfa B.gfa --ref GRCh38 --vcf diff.vcf --gaf aln.gaf

# graph-to-graph aligner, 16 threads, anchor-interval pairing
build/release/g2g A.gfa B.gfa --threads 16 --pairing interval --vcf diff.vcf --gaf aln.gaf

build/release/g2g --help     # all options: --snarls vg, --verbose-timing, --dump-snarls, weights
```

Output: a VCF-style delta (`GRAPH=A` / `GRAPH=B` records with the differing
sequence; real coordinates in `--ref` mode) and a GAF-style alignment file
(anchors, snarl pairs, allele matches). Input is GFA 1.0/1.1 (S/L/P/W lines).

## Reproducing the paper

| paper item | script | needs |
|---|---|---|
| Q1 planted-graph correctness, test suite | `ctest` in `build/release`; `experiments/exp1_correctness/run.sh` | nothing |
| Q1 pairing error vs. size (fig_pairing_error) | `python3 experiments/exp1_correctness/pairing_error_sweep.py --pairing both --anchor-lens 8,2,1 --private-fracs 0,0.1 --linker-every 4` | nothing (~15 min) |
| Q1 real HLA-DRB1 graphs | `demo/01_setup_tools.sh` then `demo/04_real_hla_demo.sh`; score with `experiments/plot/hla_truth_check.py demo/work/hla/popA.vcf demo/work/hla/popB.vcf demo/work/hla/hla_diff.vcf` | vg, minimap2, paftools (downloaded by 01) |
| Finding 1 (fig_matching) | `build/release/bench_matching > results/matching_cost_vs_k_cpp.csv` | nothing |
| Finding 2 (fig_workload) | `cd python && python3 profile_workload.py` (Gini); allele histogram as in `results/allele_count_hist.csv` | networkx |
| Q2 parallel scaling (fig_scaling) | `demo/07_parallel_scaling.sh` | a many-core node |
| Q2 throughput vs. linearize-and-map (fig_throughput) | `demo/10_throughput_baseline.sh` | vg, bcftools, minimap2 + chr22 inputs (below) |
| reference-anchored scaling (fig_refdiffscaling) | `demo/05_scaling_demo.sh`, `demo/06_scaling_limit.sh` | a large-memory node for the 8-Gbp point |
| Q3 ablations (fig_ablation) | `demo/08_ablation.sh` | nothing |
| Q4 biology (fig_biology, F_ST sweep) | `demo/09_fetch_biology_data.sh` then `demo/09_hprc_biology.sh`; sweep: `demo/hprc_enrich.py --fst-sweep 0.05,0.1,0.2,0.3,0.5` | vg, bcftools, tabix + ~1.5 GB of downloads |
| all figures | `python3 make_figures.py` (reads `results/` and `demo/sample_outputs/`) | numpy, matplotlib |

Every harness writes the TSV/CSV that `make_figures.py` reads; pass the
matching `--*-csv` flag to plot a fresh run. The committed tables in
`results/` and `demo/sample_outputs/` are the measurements used in the paper
(see the header of `make_figures.py` for the mapping). `demo/README.md`
describes each demo and its expected output; `experiments/README.md` describes
the Slurm suite.

### Public inputs used

| data | version / URL |
|---|---|
| 1000 Genomes phase 3 | release 20130502 (`ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`, `integrated_call_samples_v3.20130502.ALL.panel`), https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ |
| GRCh37 reference | `human_g1k_v37.fasta` (1000 Genomes technical/reference), or chromosome 22 from Ensembl GRCh37 (`Homo_sapiens.GRCh37.dna.chromosome.22.fa.gz`, identical sequence) |
| gnomAD-SV | v2.1 sites VCF (`gnomad_v2.1_sv.sites.vcf.gz`, GRCh37), https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/ |
| ClinVar | GRCh37 VCF, weekly release of 28 September 2026 (`##fileDate=2026-09-28`), https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz |
| HLA-DRB1 haplotypes | PGGB test data `DRB1-3123.fa.gz` (12 GRCh38 sequences) |
| tools | vg 1.75-1.76, minimap2 2.28-2.31, bcftools/samtools 1.24 (conda-forge/bioconda) |

The HPRC assembly pipeline under `experiments/exp6_hprc/` (Minigraph-Cactus
population graphs) is provided as scripts only; the paper reports no results
on Minigraph-Cactus graphs.

## Clusters

The Slurm scripts in `experiments/` and `demo/` use generic directives
(`--partition`, `--mem`, `--time`); set `--partition` (e.g. `batch` on KAUST
Ibex, where CPU jobs are submitted from `ilogin`) and `--account` for your
site, and run anything that downloads on a login node, since compute nodes
usually have no outbound network. `experiments/accl/` holds the three job
scripts used for the camera-ready re-measurements on the authors' Slurm
cluster (AMD EPYC 7763 nodes, 48-CPU / 192 GB job limit): build + tests +
pairing sweep + parallel scaling + ablation, the single-core scaling sweep,
and the real chr22 throughput and Q4 runs.

## Citation and license

MIT license (see `LICENSE`). Cite the BIBM 2026 paper (`CITATION.cff`); the
camera-ready release is archived on Zenodo (DOI in the paper's Availability
paragraph once minted).
