# PanG2G: graph-to-graph alignment of pangenome variation graphs

[![Docker Hub](https://img.shields.io/docker/pulls/ahmadshakleya/g2g?logo=docker&label=Docker%20Hub)](https://hub.docker.com/r/ahmadshakleya/g2g)
[![image size](https://img.shields.io/docker/image-size/ahmadshakleya/g2g/latest?logo=docker&label=image%20size)](https://hub.docker.com/r/ahmadshakleya/g2g)
[![container build](https://github.com/accl-kaust/pang2g/actions/workflows/container.yml/badge.svg)](https://github.com/accl-kaust/pang2g/actions/workflows/container.yml)

PanG2G compares two pangenome variation graphs directly, graph against graph,
and reports what we call the delta graph: every allele (SNP, indel or
structural variant) that is present in one graph but absent from the other,
with its sequence and, when both graphs embed a common reference, its
reference coordinate. For two graphs it does what `bcftools isec` does for
two VCFs, without flattening either graph to a linear reference first.

```bash
g2g popA.gfa popB.gfa --ref GRCh38 --vcf delta.vcf
```

## 1. What it does

A pangenome variation graph, whether it comes from `vg construct`,
Minigraph-Cactus or PGGB, stores a set of haplotypes as a graph. Sequence
shared by everyone is a chain of nodes, and every variant site is a bubble
(a snarl, in vg's vocabulary) with one branch per allele. Given two such
graphs, say one built from African and one from East Asian haplotypes,
PanG2G decomposes each into its bubbles, pairs up the bubbles that describe
the same site in both graphs, matches the alleles inside each pair by
sequence (an exact bipartite assignment, solved with LAPJV), and then emits
every allele that found no partner. Those unmatched alleles, with their
sequences, are the delta graph: the variation private to one population.

There are two ways to pair the bubbles.

The reference-anchored diff (`--ref`) pairs them by the absolute coordinate
of the site on a reference path that both graphs embed. Each allele is first
broken down into atomic SNPs and indels, so the method does not care whether
the two graphs cut the reference into different nodes, or whether one graph
spans two adjacent SNPs with a single bubble where the other has two. This is
the mode to use whenever both graphs were built on the same reference, which
is the normal case.

The graph-to-graph aligner (no `--ref`) is for graphs with no shared
reference path. It pairs bubbles by their boundary sequences, or, with
`--pairing interval`, by an order-preserving dynamic program between shared
anchors. This is the mode that the parallel work-stealing runtime (OpenMP)
accelerates.

Both modes write a VCF-style file with one record per private allele; the
columns are described in section 5, and section 4 shows a complete example.

## 2. Quick start

You need GCC 9 or newer (or Clang 10), CMake 3.18 or newer, `python3` (3.6
or newer) and bash. Nothing else is required for the build or the test
case. The Python scripts only use the standard library; numpy speeds one of
them up if it happens to be installed.

```bash
git clone https://github.com/accl-kaust/pang2g.git     # or unpack the archive you received
cd pang2g
bash run_test_case.sh            # build + three checks with known answers, about 2 minutes
```

The run ends with a summary like this one:

```
==== summary ====
  PASS  build g2g (build/release)                                     62s
  PASS  1. planted ground truth (demo/03_synthetic_demo.sh)            0s
  PASS  2. enrichment self-test (demo/hprc_enrich.py --self-test)      1s
  PASS  3. biological validation on synthetic data (demo/09_hprc_biology.sh --synthetic)     4s

ALL CHECKS PASSED
```

and leaves the binary at `build/release/g2g`. `--unit` also runs the 161
GoogleTest unit tests (CMake downloads googletest, so that needs network
access), and `--real` downloads `vg` and `minimap2` and diffs two pangenomes
built from real HLA-DRB1 haplotypes. `bash run_test_case.sh --help` lists
the options.

### Container image

The same binary, demo scripts and test case are also packaged as a Docker
image on Docker Hub:

```bash
docker pull ahmadshakleya/g2g
docker run --rm ahmadshakleya/g2g g2g --help
docker run --rm -v "$PWD":/data ahmadshakleya/g2g g2g A.gfa B.gfa --ref GRCh38 --vcf delta.vcf
docker run --rm ahmadshakleya/g2g /opt/g2g/container/selftest.sh     # the test case of section 4
```

`Dockerfile` builds it (`docker build -t g2g .`), and `container/g2g.def` is
the same recipe for Apptainer and Singularity. The engine inside is compiled
for x86-64-v3, so one image runs on any x86-64 node.

Clusters usually provide Apptainer or Singularity rather than Docker.
Building from the definition file needs root or fakeroot, which clusters
rarely allow; converting the Docker Hub image needs no privileges:

```bash
apptainer pull g2g.sif docker://ahmadshakleya/g2g      # or: apptainer build g2g.sif container/g2g.def
apptainer exec g2g.sif g2g A.gfa B.gfa --ref GRCh38 --vcf delta.vcf
apptainer exec g2g.sif /opt/g2g/container/selftest.sh
```

`singularity` takes the same commands. Both runtimes mount the home
directory inside the image but usually nothing else, so data on a scratch
file system has to be bound explicitly, and the self-test's writable copy
goes to `/tmp` unless `G2G_SELFTEST_DIR` says otherwise:

```bash
singularity run --bind /scratch/$USER g2g.sif /scratch/$USER/A.gfa /scratch/$USER/B.gfa --ref GRCh38 --vcf /scratch/$USER/delta.vcf
G2G_SELFTEST_DIR=/scratch/$USER/g2g_selftest singularity exec --bind /scratch/$USER g2g.sif /opt/g2g/container/selftest.sh
```

The images are rebuilt by `.github/workflows/container.yml` on every push:
it builds both, runs the test case inside each, pushes the Docker image to
Docker Hub and keeps `g2g.sif` as an artifact of
the run.

## 3. Running on a cluster

PanG2G needs nothing from a cluster beyond a compiler, CMake and Python, so
the usual Slurm workflow applies without any site-specific setup.

On the login node, load whatever provides the toolchain at your site and run
the test case; compiling takes about a minute, which is light enough for a
login node:

```bash
git clone https://github.com/accl-kaust/pang2g.git && cd pang2g    # or unpack the archive
module load gcc cmake python     # GCC >= 9, CMake >= 3.18 and python3, by whatever names your site uses
bash run_test_case.sh --unit     # about 3 minutes
```

To run it as a batch job instead, submit `experiments/slurm/test_case.slurm`
from the repository root. The job loads the modules named in `G2G_MODULES`
(leave it unset if your cluster needs none), builds into `build/slurm` with
`-DG2G_MARCH=x86-64-v3` and runs `run_test_case.sh`; its output lands in
`logs/`. Set `--partition` and, if your site needs it, `--account` in the
file first.

```bash
G2G_MODULES="gcc cmake python" sbatch experiments/slurm/test_case.slurm
squeue --me
tail -f logs/test_case_<jobid>.out
```

The default build uses `-march=native`, so a binary compiled on the login node may not run on a
compute node of a different CPU generation; build with
`G2G_MARCH=x86-64-v3 bash run_test_case.sh` (or pass
`-DG2G_MARCH=x86-64-v3` to CMake) for a binary that runs on any AVX2 node, or
build inside the job that will use it, as the whole-node job scripts do. Compute
nodes often have no outbound network, so anything that downloads
(`demo/01_setup_tools.sh`, `demo/04_real_hla_demo.sh`,
`demo/09_fetch_biology_data.sh`, the googletest fetch) should run on the
login node; to include the unit tests in the batch job, unpack googletest
once on the login node and point `G2G_GTEST_DIR` at it:

```bash
curl -fsSL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz | tar xz -C /path/to/scratch
G2G_GTEST_DIR=/path/to/scratch/googletest-1.14.0 sbatch experiments/slurm/test_case.slurm
```

Keep large inputs and generated graphs on scratch rather than under a home
quota; the demo harnesses write under `demo/work/`, so clone the repository
on scratch or symlink `demo/work` there. The parallel aligner uses all the
CPUs Slurm gives the job (`--cpus-per-task`); pass `--threads N` to use
fewer. Figures need numpy and matplotlib, for example in a venv:
`python3 -m venv ~/g2g-venv && . ~/g2g-venv/bin/activate && pip install numpy matplotlib`.

For the bigger runs, `experiments/ibex/` holds whole-node job scripts
(parallel scaling to 128 threads, the ablation, the 8-Gbp sweep), and
`experiments/exp*/slurm.sh` is the broader experiment suite, which expects
the conda environment described in `experiments/README.md`. All of them are
submitted from the repository root.

## 4. The test case

`run_test_case.sh` builds the engine and runs three checks whose correct
answer is known in advance. None of them needs network access or third-party
tools.

The first check (`demo/03_synthetic_demo.sh`, well under a second) is a
planted ground truth. `demo/make_synthetic_pair.py` writes two GFA graphs
that embed the same 400 bp reference path, called `GRCh38`, but cut it into
different node sizes (25 bp in A, 17 bp in B), and plants 5 SNPs in both
graphs, 3 only in A and 2 only in B. `g2g --ref GRCh38` has to report
exactly the 3 A-only and 2 B-only sites at the planted coordinates, and a
graph diffed against itself has to report nothing:

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

The second check (`demo/hprc_enrich.py --self-test`, about a second) covers
the statistic behind the biological validation. That analysis asks
whether the delta alleles lie closer than chance to structural variants that
are strongly differentiated between the two populations (high F_ST in
gnomAD-SV) and to ClinVar pathogenic sites, using a permutation null. The
self-test plants an enrichment and has to recover it (fold well above 1,
p < 0.01), and has to find none for a random control (fold close to 1).

The third check (`demo/09_hprc_biology.sh --synthetic`, about 5 seconds)
runs that whole pipeline end to end on synthetic data. A 10 Mb pair with 300
planted population-specific sites is diffed with the real binary, a gnomAD-SV
and ClinVar fixture is synthesised around the recovered deltas, and the
pipeline has to find the signal and must not find one for shuffled deltas:

```
  [PASS] delta vs high-F_ST (signal)        fold=20.84x  p=0.000999
  [PASS] delta vs ClinVar pathogenic        fold=41.14x  p=0.000999
  [PASS] shuffled control (null)            fold=1.30x   p=0.23
```

With `--unit`, the 161 GoogleTests run as well. They cover the GFA reader,
snarl decomposition (checked against `vg snarls` when `vg` is on the PATH),
allele matching (LAPJV against brute force), both pairing modes, the parallel
runtime (its results must be identical to the sequential aligner at every
thread count), the reference-anchored diff (a merged bubble must equal the
split bubbles), and planted-graph recall and precision. The 9 tests that need
`vg` or `minimap2` are skipped when those tools are absent.

With `--real`, `demo/04_real_hla_demo.sh` downloads 12 real HLA-DRB1
haplotypes (the PGGB test data), builds two population graphs with `minimap2`
and `vg construct` the way real pangenomes are built, diffs them and
self-diffs one of them. Expect `a_only=349 b_only=12` and an empty self-diff;
`demo/logs/04_real_hla_demo.log` is a captured run.

Every step's full output is kept in `logs/test_case/`. `demo/README.md`
describes all the demos (scaling to 8 Gbp, parallel speed-up, the ablation,
throughput against linearize-and-map) with their expected numbers.

## 5. Using PanG2G on your own graphs

### Input

Two GFA 1.0 or 1.1 files; `S`, `L`, `P` and `W` lines are read. For the
reference-anchored diff both graphs must embed the reference as a path with
the same name: a `P` line named for example `GRCh38` or `chr22`, or a `W`
line whose PanSN name has that sample or contig (`--ref GRCh38` matches
`GRCh38#0#chr22`). A reference that spans several contigs is diffed contig
by contig. The graphs must be acyclic, which is true of reference-anchored
graphs from `vg construct` and Minigraph-Cactus; `vg msga` and some PGGB
graphs contain cycles.

To build two population graphs from a phased panel VCF, for one chromosome:

```bash
bcftools view -S afr_samples.txt -r 22 --force-samples -c1 -Oz panel.vcf.gz > afr.vcf.gz && tabix afr.vcf.gz
vg construct -r ref.fa -v afr.vcf.gz -R 22 | vg convert -f - > afr.gfa      # likewise eas.gfa
g2g afr.gfa eas.gfa --ref 22 --vcf delta.vcf --gaf delta.gaf
```

`demo/09_hprc_biology.sh` and `demo/10_throughput_baseline.sh` do exactly
this for the chromosome-22 runs, and `demo/09_fetch_biology_data.sh`
downloads the public inputs.

### Command line

```
g2g <graph_a.gfa> <graph_b.gfa> [options]
```

| option | meaning |
|---|---|
| `--ref <name>` | reference-anchored diff along the path `<name>`; use it whenever both graphs embed a reference |
| `--vcf <file>` | write the delta (private alleles) as a VCF-style file |
| `--gaf <file>` | write the alignment (snarl pairs and allele matches) in a GAF-like text format; default stdout |
| `--threads <n>` | worker threads for the graph-to-graph aligner (default: all available; `G2G_NUM_THREADS` works too) |
| `--pairing boundary\|interval` | bubble pairing for the graph-to-graph aligner (no effect with `--ref`): boundary-sequence key (default) or anchor-interval dynamic program |
| `--snarls builtin\|vg` | bubble finder: the built-in superbubble scan (default) or `vg snarls` (`vg` on the PATH or `$G2G_VG`) |
| `--sequential` | single-threaded greedy aligner (no effect with `--ref`) |
| `--verbose-timing` | per-phase timings, peak memory and scheduler statistics on stderr |
| `--dump-snarls <file>` | write every bubble of both graphs, for debugging |
| `--alpha --beta --gamma --gap` | scoring weights of the graph-to-graph aligner |

The exit status is 0 on success, 1 for bad arguments, unreadable input, a
cyclic graph or a reference path missing from a graph, and 2 when no bubble
could be paired at all, which usually means a wrong `--ref` name or
unrelated inputs. Progress and the counts (`snarl_pairs`, `a_only`,
`b_only`) go to stderr.

### Output format

`--vcf` writes one record per private allele. The header line
`##fileformat=G2GdeltaVCF-0.2` marks it as VCF-style rather than VCF: it has
the standard eight columns, but in `--ref` mode POS is 0-based (the
reference offset of the site, as in BED), so add 1 before comparing with a
1-based VCF.

| column | `--ref` mode | graph-to-graph mode |
|---|---|---|
| CHROM | reference path (contig) name | snarl id (`source>sink` node ids) |
| POS | 0-based reference offset of the site | 0 |
| ID | `A_<i>` / `B_<i>` | same |
| REF | reference sequence of the atomic variant; `.` for an insertion | `.` |
| ALT | alternate sequence; `.` for a deletion | the sequence private to one graph |
| INFO | `GRAPH=A` or `GRAPH=B` (which input the allele is private to), `SNARL=<contig>:<pos>`, `LEN=` ALT length | `GRAPH=`, `SNARL=`, `LEN=` |

Alleles are decomposed into atomic SNPs and indels before comparison, so a
bubble that spans two adjacent SNPs in one graph matches two separate bubbles
in the other. Indels are not left-normalised, so in repeats a call may differ
by a base from what `bcftools norm` would give.

### Limits

Only acyclic graphs are accepted; a cyclic input stops with
`Graph contains cycles - DAG required`. Memory grows linearly with graph
size: roughly 30 GB per gigabase of reference in `--ref` mode at one variant
per 200 bp (we measured 242 GB for an 8 Gbp pair with 971 million nodes).
Human chromosome 22 takes seconds. The `--ref` diff itself is
single-threaded; the parallel runtime accelerates the graph-to-graph
aligner.

## 6. Larger experiments

The scripts below scale the checks of the test case up, and each one writes
a TSV or CSV with its measurements into `demo/work/` or `results/phase5/`.

| experiment | script | needs |
|---|---|---|
| planted-graph correctness sweep | `bash run_test_case.sh --unit`; `experiments/exp1_correctness/run.sh` | nothing |
| pairing error against graph size | `python3 experiments/exp1_correctness/pairing_error_sweep.py --pairing both --anchor-lens 8,2,1 --private-fracs 0,0.1 --linker-every 4` | nothing (about 15 min) |
| real HLA-DRB1 graphs | `bash run_test_case.sh --real`; score with `experiments/analysis/hla_truth_check.py demo/work/hla/popA.vcf demo/work/hla/popB.vcf demo/work/hla/hla_diff.vcf` | vg, minimap2, paftools (downloaded by `demo/01_setup_tools.sh`) |
| matching cost against allele count | `build/release/bench_matching > matching_cost_vs_k.csv` | nothing |
| workload characterisation (prototype) | `cd python && python3 -u profile_workload.py` (about 25 min, pure Python) | networkx |
| parallel scaling | `demo/07_parallel_scaling.sh` (whole node: `experiments/ibex/scaling_64.slurm`) | a many-core node |
| throughput against linearize-and-map | `demo/10_throughput_baseline.sh` | vg, bcftools, minimap2 and the chr22 inputs below |
| reference-anchored scaling | `demo/05_scaling_demo.sh`, `demo/06_scaling_limit.sh` (whole node: `experiments/ibex/scaling_8gbp.slurm`) | a large-memory node for the 8-Gbp point |
| ablations | `demo/08_ablation.sh` | nothing |
| biological validation and F_ST sweep | `demo/09_fetch_biology_data.sh` then `demo/09_hprc_biology.sh`; sweep: `demo/hprc_enrich.py --fst-sweep 0.05,0.1,0.2,0.3,0.5` | vg, bcftools, tabix and about 1.5 GB of downloads |

The demo scripts use the binary that `run_test_case.sh` built
(`build/release/g2g`; `build/g2g` from `demo/02_build_g2g.sh` is found too,
and `G2G=/path/to/g2g` overrides both). `demo/README.md` describes each demo
and its expected output, `experiments/README.md` the Slurm experiment suite,
and `experiments/ibex/` keeps the whole-node job scripts for the biggest
runs.

### Public inputs for the real-data runs

| data | version / URL |
|---|---|
| 1000 Genomes phase 3 | release 20130502 (`ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`, `integrated_call_samples_v3.20130502.ALL.panel`), https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ |
| GRCh37 reference | `human_g1k_v37.fasta` (1000 Genomes technical/reference), or chromosome 22 from Ensembl GRCh37 (`Homo_sapiens.GRCh37.dna.chromosome.22.fa.gz`, identical sequence) |
| gnomAD-SV | v2.1 sites VCF (`gnomad_v2.1_sv.sites.vcf.gz`, GRCh37), https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/ |
| ClinVar | GRCh37 VCF, weekly release of 28 September 2026 (`##fileDate=2026-09-28`), https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz |
| HLA-DRB1 haplotypes | PGGB test data `DRB1-3123.fa.gz` (12 GRCh38 sequences) |
| tools | vg 1.75 to 1.76, minimap2 2.28 to 2.31, bcftools/samtools 1.24 (conda-forge/bioconda) |

The HPRC assembly pipeline under `experiments/exp6_hprc/` (Minigraph-Cactus
population graphs) is provided as scripts only; we have not run it to
completion ourselves.

## 7. Repository layout

```
run_test_case.sh       build + the test case of section 4 (--unit, --real, --quick)
CMakeLists.txt         build (options: G2G_MARCH, G2G_BUILD_TESTS, the five ablation flags)
Dockerfile, container/ container recipes (Docker; Apptainer/Singularity) and the in-image self-test
.github/workflows/     CI: builds both container images, runs the test case inside them, publishes
include/g2g/, src/     C++20 engine (static library + the `g2g` command)
tests/                 GoogleTest suite (fetched by CMake; -DG2G_BUILD_TESTS=OFF to skip)
tools/                 bench_matching: matching cost vs. allele count
demo/                  end-to-end demos and harnesses: synthetic truth, real HLA-DRB1
                       graphs, scaling sweeps, parallel scaling, ablation, biological
                       validation, throughput; demo/README.md documents each with
                       expected output
demo/sample_outputs/   example outputs of the synthetic and HLA demos
experiments/           Slurm experiment suite (exp1-exp6), the test-case job (slurm/),
                       whole-node job scripts for the biggest runs (ibex/), generators
                       and analysis scripts (analysis/); experiments/README.md
python/g2g/            Python reference prototype (correctness oracle; slow)
python/                workload characterisation scripts (profile_workload.py, synth_generator.py)
logs/                  Slurm output and test-case logs (git-ignored)
```

## 8. Troubleshooting

`cmake: command not found`, or the compiler rejects C++20: load your site's
toolchain modules (`module avail gcc cmake python` lists them) or install
GCC 9 or newer and CMake 3.18 or newer.

`/usr/bin/env: 'bash\r': No such file or directory`: the files were copied
with Windows line endings. A fresh `git clone` is already correct (the
repository forces LF); for a copied tree run
`find . -name '*.sh' -o -name '*.py' -o -name '*.slurm' | xargs sed -i 's/\r$//'`.

`Illegal instruction` when running a binary built elsewhere: it was built
with `-march=native` on a different CPU. Rebuild with `G2G_MARCH=x86-64-v3`,
or build in the job that uses it.

CMake fails to fetch googletest (no network): build without the tests
(`-DG2G_BUILD_TESTS=OFF`, which is what `run_test_case.sh` does without
`--unit`), or point `G2G_GTEST_DIR` or `-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST`
at an unpacked googletest-1.14.0.

`fatal error: span: No such file or directory`: should not happen. GCC 9's
libstdc++ has no `<span>`, and `include/g2g/span.hpp` substitutes a small
shim.

`ERROR: reference path 'X' not found in both graphs`: the two graphs name
their paths differently. `grep '^[PW]' A.gfa | cut -f2 | head` shows the
names; `--ref` accepts a P-line name, or the sample or contig of a PanSN
W-line.

`Graph contains cycles - DAG required`: the graph has a cycle, for example
from `vg msga` or a PGGB build. PanG2G needs an acyclic, reference-anchored
graph such as `vg construct` or Minigraph-Cactus output.

Exit status 2 with `snarl_pairs=0`: nothing could be paired. Check the
`--ref` name; in graph-to-graph mode the inputs may share no boundary
sequences, in which case try `--ref`, `--pairing interval`, or
`--dump-snarls` to look at the bubbles.

Positions look off by one against a VCF: the `--vcf` output uses 0-based POS
(section 5); add 1.

`sbatch: error: Unable to open output file`: submit from the repository
root, so that `logs/` exists.

## 9. Citation

If you use PanG2G, please cite the paper:

```bibtex
@inproceedings{shakleya2026pang2g,
  author    = {Shakleya, Ahmad and Fahmy, Suhaib A.},
  title     = {{PanG2G}: Scalable Base-Pair-Resolution Graph-to-Graph Alignment of Pangenome Graphs},
  booktitle = {Proceedings of the IEEE International Conference on Bioinformatics and Biomedicine (BIBM)},
  year      = {2026},
  publisher = {IEEE},
  address   = {Dallas, TX, USA}
}
```
