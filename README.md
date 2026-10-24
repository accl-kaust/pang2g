# PanG2G: graph-to-graph alignment of pangenome variation graphs

PanG2G compares two pangenome variation graphs **directly**, graph against
graph, and reports the **delta graph**: every allele (SNP, indel, structural
variant) that is present in one graph but absent from the other, as a
VCF-style table with its sequence and, when both graphs embed a common
reference, its reference coordinate.  It does for two graphs what
`bcftools isec` does for two VCFs, without first flattening either graph to
a linear reference.

```bash
g2g popA.gfa popB.gfa --ref GRCh38 --vcf delta.vcf
```

It is the code behind

> A. Shakleya and S. A. Fahmy, *PanG2G: Scalable Base-Pair-Resolution
> Graph-to-Graph Alignment of Pangenome Graphs*, IEEE BIBM 2026.

Contents: [What it does](#1-what-it-does) ·
[Quick start](#2-quick-start-any-linux-machine) ·
[Running on a cluster](#3-running-on-a-cluster) ·
[The test case](#4-the-test-case) ·
[Your own graphs](#5-using-pang2g-on-your-own-graphs) ·
[Reproducing the paper](#6-reproducing-the-paper) ·
[Repository layout](#7-repository-layout) ·
[Troubleshooting](#8-troubleshooting) ·
[Citation](#9-citation-and-license)

---

## 1. What it does

A pangenome variation graph (from `vg construct`, Minigraph-Cactus, ...)
stores a set of haplotypes as a graph: shared sequence is a chain of nodes,
and every variant site is a **bubble** (a *snarl*) with one branch per
**allele**.  Given two such graphs, say one built from African and one from
East Asian haplotypes, PanG2G

1. decomposes each graph into its bubbles,
2. pairs up the bubbles that describe the same site in both graphs,
3. matches the alleles inside each paired bubble by sequence (an exact
   bipartite assignment, LAPJV), and
4. emits every allele that found no partner, with its sequence: the
   **delta graph**, i.e. the variation private to one population.

Two ways of pairing bubbles are available:

| mode | command | how bubbles are paired | when to use |
|---|---|---|---|
| **reference-anchored diff** | `g2g A.gfa B.gfa --ref GRCh38 --vcf delta.vcf` | by the absolute coordinate of the site on a reference path embedded in both graphs; each allele is first decomposed into atomic SNPs/indels, so it is robust to the two graphs cutting the reference into different nodes | the normal case: graphs built on the same reference (vg construct, Minigraph-Cactus, PGGB graphs that carry the reference as a path) |
| **graph-to-graph aligner** | `g2g A.gfa B.gfa --threads 16 --vcf delta.vcf` | by bubble boundary sequences (default) or by an anchor-guided order-preserving dynamic program (`--pairing interval`); work-stealing parallel runtime (OpenMP) | graphs with no shared reference path |

The output is a VCF-style file with one record per private allele; see
[Output format](#output-format) for the columns.  Section
[4](#4-the-test-case) shows a complete example.

## 2. Quick start (any Linux machine)

Requirements: GCC >= 9 (or Clang >= 10), CMake >= 3.18, `python3` (>= 3.6),
`bash`.  Nothing else is needed for the build or the test case; the Python
scripts use only the standard library (numpy, if present, speeds one of them
up).

```bash
git clone https://github.com/ahmadshakleya/G2G.git     # or unpack the archive you received
cd G2G
bash run_test_case.sh            # build + three checks with known answers, ~2 minutes
```

The run ends with

```
==== summary ====
  PASS  build g2g (build/release)                                     62s
  PASS  1. planted ground truth (demo/03_synthetic_demo.sh)            0s
  PASS  2. enrichment self-test (demo/hprc_enrich.py --self-test)      1s
  PASS  3. Q4 pipeline on synthetic data (demo/09_hprc_biology.sh --synthetic)     4s

ALL CHECKS PASSED
```

and leaves the binary at `build/release/g2g`.  Add `--unit` to also run the
161 GoogleTest unit tests (CMake downloads googletest, so this needs network
access), and `--real` to download `vg` + `minimap2` and diff two pangenomes
built from real HLA-DRB1 haplotypes.  `bash run_test_case.sh --help` lists
every option.

### Container image

`Dockerfile` builds an image with the `g2g` command, the demo scripts and the
test case; `container/g2g.def` is the same recipe for Apptainer/Singularity.
Both compile the engine for x86-64-v3, so one image runs on any x86-64 node.

```bash
docker build -t g2g .                                  # or: docker pull ghcr.io/ahmadshakleya/g2g:latest
docker run --rm g2g g2g --help
docker run --rm -v "$PWD":/data g2g g2g A.gfa B.gfa --ref GRCh38 --vcf delta.vcf
docker run --rm g2g /opt/g2g/container/selftest.sh    # the test case of section 4
```

On a cluster use Apptainer or Singularity.  Building from the definition
file needs root or fakeroot; converting the Docker image does not:

```bash
apptainer pull g2g.sif docker://ghcr.io/ahmadshakleya/g2g:latest   # or: apptainer build g2g.sif container/g2g.def
apptainer exec g2g.sif g2g A.gfa B.gfa --ref GRCh38 --vcf delta.vcf
apptainer exec g2g.sif /opt/g2g/container/selftest.sh
```

(`singularity` takes the same commands.)  The workflow in
`.github/workflows/container.yml` builds both images on every push, runs the
test case inside each, pushes the Docker image to ghcr.io and attaches
`g2g.sif` to the run as a downloadable artifact.  While the repository is
private, pulling from ghcr.io needs a GitHub token (`docker login ghcr.io`,
`apptainer registry login docker://ghcr.io`); the alternative is to download
the `g2g.sif` artifact and copy the file to the cluster.

## 3. Running on a cluster

PanG2G needs only a C++ compiler, CMake and Python, so it fits the usual
Slurm workflow without any site-specific setup.

**On the login node:**

```bash
git clone https://github.com/ahmadshakleya/G2G.git && cd G2G     # or unpack the archive
module load gcc cmake python     # whatever provides GCC >= 9, CMake >= 3.18 and python3 at your site
bash run_test_case.sh --unit     # ~3 minutes; the compile is light enough for a login node
```

**As a batch job** (`experiments/slurm/test_case.slurm`; the output file
lands in `logs/`):

```bash
G2G_MODULES="gcc cmake python" sbatch experiments/slurm/test_case.slurm   # set --partition / --account in the file
squeue --me
tail -f logs/test_case_<jobid>.out
```

The job loads the modules named in `G2G_MODULES` (leave it unset if your
cluster needs none), builds into `build/slurm` with `-DG2G_MARCH=x86-64-v3`,
and runs `run_test_case.sh`.  Points that matter on most clusters:

- **Which binary runs where.**  The default build uses `-march=native`, so a
  binary compiled on the login node may not run on a compute node of a
  different CPU generation.  Build with
  `G2G_MARCH=x86-64-v3 bash run_test_case.sh` (or
  `cmake -DG2G_MARCH=x86-64-v3 ...`) for a binary that runs on any AVX2 node,
  or build inside the job that will use it, as the paper's job scripts do.
- **No network on compute nodes.**  Anything that downloads
  (`demo/01_setup_tools.sh`, `demo/04_real_hla_demo.sh`,
  `demo/09_fetch_biology_data.sh`, the googletest fetch) runs on the login
  node.  To include the unit tests in the batch job, unpack googletest once
  on the login node and point `G2G_GTEST_DIR` at it:

  ```bash
  curl -fsSL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz | tar xz -C /path/to/scratch
  G2G_GTEST_DIR=/path/to/scratch/googletest-1.14.0 sbatch experiments/slurm/test_case.slurm
  ```

- **Disk.**  Keep large inputs and generated graphs on a scratch file system
  rather than under a home quota; the demo harnesses write under
  `demo/work/`, so clone the repository on scratch or symlink `demo/work`
  there.
- **Threads.**  The parallel aligner uses all CPUs Slurm gives the job
  (`--cpus-per-task`); pass `--threads N` to the binary to use fewer.
- **Figures** need numpy + matplotlib:
  `python3 -m venv ~/g2g-venv && . ~/g2g-venv/bin/activate && pip install numpy matplotlib`.
- **Larger runs.**  `experiments/accl/` and `experiments/ibex/` hold the
  whole-node job scripts used for the paper's measurements on the authors'
  clusters (parallel scaling to 128 threads, ablation, the 8-Gbp sweep);
  `experiments/exp*/slurm.sh` is the broader experiment suite, which expects
  the conda environment described in `experiments/README.md`.  All are
  submitted from the repository root; set `--partition` and `--account` for
  your site.

## 4. The test case

`run_test_case.sh` builds the engine and runs three checks whose correct
answer is known in advance.  None needs network access or third-party
tools.

**Check 1: planted ground truth** (`demo/03_synthetic_demo.sh`, < 1 s).
`demo/make_synthetic_pair.py` writes two GFA graphs that embed the *same*
400 bp reference path (`GRCh38`) but cut it into different node sizes
(25 bp in A, 17 bp in B), and plants 5 SNPs in both graphs, 3 only in A
and 2 only in B.  `g2g --ref GRCh38` must report exactly the 3 A-only and 2
B-only sites at the planted coordinates, and a graph diffed against itself
must report nothing:

```
[g2g] GA: 38 nodes  GB: 38 nodes
[g2g] reference-coordinate diff along path 'GRCh38'
[g2g] total_score=5  snarl_pairs=5  a_only=3  b_only=2

#CHROM  POS   ID    REF  ALT  QUAL  FILTER  INFO
GRCh38  27    A_0   A    T    .     PASS    GRAPH=A;SNARL=GRCh38:27;LEN=1
GRCh38  265   A_1   C    A    .     PASS    GRAPH=A;SNARL=GRCh38:265;LEN=1
GRCh38  333   A_2   G    C    .     PASS    GRAPH=A;SNARL=GRCh38:333;LEN=1
GRCh38  129   B_0   G    T    .     PASS    GRAPH=B;SNARL=GRCh38:129;LEN=1
GRCh38  163   B_1   G    T    .     PASS    GRAPH=B;SNARL=GRCh38:163;LEN=1

  planted : shared=5 A_only=3 B_only=2
  reported: A_only=3 B_only=2
  RESULT  : PASS - diff matches ground truth exactly

=== 5. Self-diff (graph vs itself) - must be zero differences ===
  a_only=0  b_only=0
```

**Check 2: enrichment statistic self-test** (`demo/hprc_enrich.py
--self-test`, ~1 s).  The paper's biological validation asks whether the
delta alleles lie closer than chance to population-differentiated structural
variants (high F_ST in gnomAD-SV) and to ClinVar pathogenic sites, using a
permutation null.  The self-test plants an enrichment and must recover it
(fold >> 1, p < 0.01), and must report none for a random control (fold ~ 1).

**Check 3: the Q4 pipeline end to end on synthetic data**
(`demo/09_hprc_biology.sh --synthetic`, ~5 s).  A 10 Mb pair with 300
planted population-specific sites is diffed with the real binary, a gnomAD-SV
and ClinVar fixture is synthesised around the recovered deltas, and the
pipeline must find the signal and must *not* find one for shuffled deltas:

```
  [PASS] delta vs high-F_ST (signal)        fold=20.84x  p=0.000999
  [PASS] delta vs ClinVar pathogenic        fold=41.14x  p=0.000999
  [PASS] shuffled control (null)            fold=1.30x   p=0.23
```

**Unit tests** (`--unit`): 161 GoogleTests over the GFA reader, snarl
decomposition (checked against `vg snarls` when `vg` is on `PATH`), allele
matching (LAPJV against brute force), both pairing modes, the parallel
runtime (results identical to the sequential aligner at every thread
count), the reference-anchored diff (merged bubble = split bubbles), and
planted-graph recall/precision.  The 9 tests that need `vg` or `minimap2`
are skipped when those tools are absent.

**Real data** (`--real`, needs network): `demo/04_real_hla_demo.sh`
downloads 12 real HLA-DRB1 haplotypes (PGGB test data), builds two
population graphs with `minimap2` + `vg construct` the way real pangenomes
are built, diffs them and self-diffs one of them (expect `a_only=349
b_only=12` and a self-diff of 0; see `demo/logs/04_real_hla_demo.log`).

Every step's full output is kept in `logs/test_case/`.  `demo/README.md`
walks through all the demos (scaling to 8 Gbp, parallel speed-up, ablation,
throughput against linearize-and-map) with their expected numbers.

## 5. Using PanG2G on your own graphs

### Input

Two GFA 1.0/1.1 files (`S`, `L`, `P` and `W` lines are read).  For the
reference-anchored diff both graphs must embed the reference as a path with
the same name: a `P` line named e.g. `GRCh38` or `chr22`, or a `W` line whose
PanSN name has that sample or contig (`--ref GRCh38` matches
`GRCh38#0#chr22`).  A reference that spans several contigs is diffed contig
by contig.  Graphs must be acyclic (true of reference-anchored graphs from
`vg construct` and Minigraph-Cactus; `vg msga` and some PGGB graphs contain
cycles).

Building two population graphs from a phased panel VCF, for one chromosome:

```bash
bcftools view -S afr_samples.txt -r 22 --force-samples -c1 -Oz panel.vcf.gz > afr.vcf.gz && tabix afr.vcf.gz
vg construct -r ref.fa -v afr.vcf.gz -R 22 | vg convert -f - > afr.gfa      # likewise eas.gfa
g2g afr.gfa eas.gfa --ref 22 --vcf delta.vcf --gaf delta.gaf
```

(`demo/09_hprc_biology.sh` and `demo/10_throughput_baseline.sh` do exactly
this for the paper's chromosome-22 runs; `demo/09_fetch_biology_data.sh`
downloads the public inputs.)

### Command line

```
g2g <graph_a.gfa> <graph_b.gfa> [options]
```

| option | meaning |
|---|---|
| `--ref <name>` | reference-anchored diff along the path `<name>` (recommended whenever both graphs embed a reference) |
| `--vcf <file>` | write the delta (private alleles) as a VCF-style file |
| `--gaf <file>` | write the alignment (snarl pairs and allele matches) in a GAF-like text format; default stdout |
| `--threads <n>` | worker threads for the graph-to-graph aligner (default: all available; also `G2G_NUM_THREADS`) |
| `--pairing boundary\|interval` | bubble pairing for the graph-to-graph aligner (no effect with `--ref`): boundary-sequence key (default) or anchor-interval dynamic program |
| `--snarls builtin\|vg` | bubble finder: the built-in superbubble scan (default) or `vg snarls` (`vg` on `PATH` or `$G2G_VG`) |
| `--sequential` | single-threaded greedy aligner (no effect with `--ref`) |
| `--verbose-timing` | per-phase timings, peak memory and scheduler statistics on stderr |
| `--dump-snarls <file>` | write every bubble of both graphs (for debugging) |
| `--alpha --beta --gamma --gap` | scoring weights of the graph-to-graph aligner |

Exit status: 0 success; 1 bad arguments, unreadable input, cyclic graph or
reference path missing from a graph; 2 no bubble could be paired (usually
wrong `--ref` name or unrelated inputs).  Progress and counts
(`snarl_pairs`, `a_only`, `b_only`) go to stderr.

### Output format

`--vcf` writes one record per private allele.  The header line
`##fileformat=G2GdeltaVCF-0.2` marks it as VCF-*style*: it has the standard
eight columns but, in `--ref` mode, **POS is 0-based** (the reference offset
of the site, as in BED), so add 1 before comparing with a 1-based VCF.

| column | `--ref` mode | graph-to-graph mode |
|---|---|---|
| CHROM | reference path (contig) name | snarl id (`source>sink` node ids) |
| POS | 0-based reference offset of the site | 0 |
| ID | `A_<i>` / `B_<i>` | same |
| REF | reference sequence of the atomic variant; `.` for an insertion | `.` |
| ALT | alternate sequence; `.` for a deletion | the sequence private to one graph |
| INFO | `GRAPH=A` or `GRAPH=B` (which input the allele is private to); `SNARL=<contig>:<pos>`; `LEN=` ALT length | `GRAPH=`, `SNARL=`, `LEN=` |

Alleles are decomposed into atomic SNPs and indels before comparison, so a
bubble that spans two adjacent SNPs in one graph matches two separate bubbles
in the other.  Indels are not left-normalised; in repeats a call may differ
by a base from `bcftools norm`.

### Limits

- Acyclic graphs only (a cyclic input stops with `Graph contains cycles -
  DAG required`).
- Memory is linear in graph size: about 30 GB per 1 Gbp of reference in
  `--ref` mode with one variant per 200 bp (measured: 8 Gbp / 971 M nodes
  in 242 GB).  Human chromosome 22 takes seconds.
- The `--ref` diff is single-threaded; the parallel runtime accelerates the
  graph-to-graph aligner.

## 6. Reproducing the paper

Everything reported in the paper (numbers, tables, figures) can be
regenerated from this repository.  The committed tables in `results/` and
`demo/sample_outputs/` are the measurements used in the paper, and
`python3 make_figures.py` redraws every figure in `latex/figs/` from them
(see the header of `make_figures.py` for the figure-to-table mapping).

| paper item | script | needs |
|---|---|---|
| Q1 planted-graph correctness, test suite | `bash run_test_case.sh --unit`; `experiments/exp1_correctness/run.sh` | nothing |
| Q1 pairing error vs. size (fig_pairing_error) | `python3 experiments/exp1_correctness/pairing_error_sweep.py --pairing both --anchor-lens 8,2,1 --private-fracs 0,0.1 --linker-every 4` | nothing (~15 min) |
| Q1 real HLA-DRB1 graphs | `bash run_test_case.sh --real`; score with `experiments/plot/hla_truth_check.py demo/work/hla/popA.vcf demo/work/hla/popB.vcf demo/work/hla/hla_diff.vcf` | vg, minimap2, paftools (downloaded by `demo/01_setup_tools.sh`) |
| Finding 1 (fig_matching) | `build/release/bench_matching > results/matching_cost_vs_k_cpp.csv` | nothing |
| Finding 2 (fig_workload) | `cd python && python3 profile_workload.py` (Gini); allele histogram as in `results/allele_count_hist.csv` | networkx |
| Q2 parallel scaling (fig_scaling) | `demo/07_parallel_scaling.sh` (the paper's job: `experiments/ibex/scaling_64.slurm`) | a many-core node |
| Q2 throughput vs. linearize-and-map (fig_throughput) | `demo/10_throughput_baseline.sh` | vg, bcftools, minimap2 + chr22 inputs (below) |
| reference-anchored scaling (fig_refdiffscaling) | `demo/05_scaling_demo.sh`, `demo/06_scaling_limit.sh` (the paper's job: `experiments/ibex/scaling_8gbp.slurm`) | a large-memory node for the 8-Gbp point |
| Q3 ablations (fig_ablation) | `demo/08_ablation.sh` | nothing |
| Q4 biology (fig_biology, F_ST sweep) | `demo/09_fetch_biology_data.sh` then `demo/09_hprc_biology.sh`; sweep: `demo/hprc_enrich.py --fst-sweep 0.05,0.1,0.2,0.3,0.5` | vg, bcftools, tabix + ~1.5 GB of downloads |
| all figures | `python3 make_figures.py` | numpy, matplotlib |

Every harness writes the TSV/CSV that `make_figures.py` reads; pass the
matching `--*-csv` flag to plot a fresh run.  `demo/README.md` describes
each demo and its expected output; `experiments/README.md` describes the
Slurm experiment suite.  `experiments/accl/` and `experiments/ibex/` hold
the job scripts used for the paper's measurements on the authors' clusters
(the camera-ready re-measurements on 48-CPU / 192 GB jobs, and the
64/128-thread scaling, ablation and 8-Gbp results on whole nodes).

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

## 7. Repository layout

```
run_test_case.sh       build + the test case of section 4 (--unit, --real, --quick)
CMakeLists.txt         build (options: G2G_MARCH, G2G_BUILD_TESTS, the five ablation flags)
Dockerfile, container/ container recipes (Docker; Apptainer/Singularity) and the in-image self-test
.github/workflows/     CI: builds both container images and runs the test case inside them
include/g2g/, src/     C++20 engine (static library + the `g2g` command)
tests/                 GoogleTest suite (fetched by CMake; -DG2G_BUILD_TESTS=OFF to skip)
tools/                 bench_matching: matching cost vs. allele count (fig_matching)
demo/                  end-to-end demos and harnesses: synthetic truth, real HLA-DRB1
                       graphs, scaling sweeps, parallel scaling, ablation, Q4 biology,
                       throughput; demo/README.md documents each with expected output
demo/sample_outputs/   result tables and figures of the real-data runs
experiments/           Slurm experiment suite (exp1-exp6), the test-case job (slurm/),
                       the job scripts used for the paper's measurements (accl/, ibex/),
                       generators and plot scripts; experiments/README.md
results/               measured result tables behind the paper's figures
make_figures.py        regenerates every figure in latex/figs/ from results/ + demo/sample_outputs/
python/g2g/            Python reference prototype (correctness oracle; slow)
python/                workload characterisation scripts (profile_workload.py, synth_generator.py)
latex/                 the paper
logs/                  Slurm output and test-case logs (git-ignored)
```

## 8. Troubleshooting

| symptom | cause / fix |
|---|---|
| `cmake: command not found`, or the compiler rejects C++20 | load your site's toolchain modules (`module avail gcc cmake python` lists them) or install GCC >= 9 and CMake >= 3.18 |
| `/usr/bin/env: 'bash\r': No such file or directory` | the files were copied with Windows line endings; run `find . -name '*.sh' -o -name '*.py' -o -name '*.slurm' \| xargs sed -i 's/\r$//'` (a fresh `git clone` is already correct: `.gitattributes` forces LF) |
| `Illegal instruction` when running a binary built elsewhere | it was built with `-march=native` on a different CPU; rebuild with `G2G_MARCH=x86-64-v3` (or in the job that uses it) |
| CMake fails fetching googletest (no network) | build without tests: `-DG2G_BUILD_TESTS=OFF` (what `run_test_case.sh` does without `--unit`), or point `G2G_GTEST_DIR` / `-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST` at an unpacked googletest-1.14.0 |
| `fatal error: span: No such file or directory` | should not happen: GCC 9's libstdc++ lacks `<span>` and `include/g2g/span.hpp` substitutes a shim; report it if it does |
| `ERROR: reference path 'X' not found in both graphs` | the graphs' path names differ; `grep '^[PW]' A.gfa \| cut -f2 \| head` shows them; `--ref` accepts a P-line name, or the sample or contig of a PanSN W-line |
| `Graph contains cycles - DAG required` | the graph has a cycle (e.g. from `vg msga` or a PGGB build); PanG2G needs an acyclic, reference-anchored graph such as `vg construct` or Minigraph-Cactus output |
| exit status 2, `snarl_pairs=0` | nothing could be paired: wrong `--ref` name, or in graph-to-graph mode the inputs share no boundary sequences; try `--ref`, `--pairing interval`, or `--dump-snarls` to inspect |
| positions look off by one against a VCF | `--vcf` output uses 0-based POS (section 5); add 1 |
| `sbatch: error: Unable to open output file` | submit from the repository root so `logs/` exists |

## 9. Citation and license

MIT license (see `LICENSE`).  Please cite the BIBM 2026 paper
(`CITATION.cff`); the camera-ready release is archived on Zenodo (DOI in the
paper's Availability paragraph once minted).
