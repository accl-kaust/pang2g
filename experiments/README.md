# The experiment suite

Six experiments for the Graph-to-Graph (G2G) pangenome aligner: correctness,
baseline consistency, throughput, parallel scaling, an ablation study, and
biological validation on HPRC data.

Everything here runs on a Slurm cluster with conda or mamba and no sudo. It
was developed on a 64-core AMD EPYC 7763 node; set `--partition` and
`--account` in the `slurm.sh` files for your site. All jobs record wall time
and peak memory.

If you are on a cluster for the first time, start with the top-level test
case rather than this suite: `sbatch experiments/slurm/test_case.slurm`, as
described in the root `README.md`. The scripts under `experiments/accl/` and
`experiments/ibex/` are the whole-node jobs behind the paper's measurements
on the two clusters we used; they need no conda environment, only a
compiler, CMake and Python.

## 0. One-time setup

Run these once from the repository root, on a login node (compute nodes
usually have no outbound network):

```bash
# 1. Create the conda env (cmake, compilers, ninja, minimap2, bedtools,
#    samtools, bcftools, seqkit, k8, plus the Python analysis stack)
bash experiments/00_setup_env.sh
conda activate g2g_eval

# 2. Build the production binary + 5 ablation variants
bash experiments/01_build_all.sh
```

`01_build_all.sh` activates `g2g_eval` itself and checks the toolchain
first, so it also works if you forget to activate. It produces:

| Build dir | Binary | Purpose |
|---|---|---|
| `build/release` | `g2g` | Production (all optimizations) |
| `build/abl_bfs_layout` | `g2g` | BFS layout (no van Emde Boas) |
| `build/abl_steal_random` | `g2g` | Random steal (no largest-first) |
| `build/abl_naive_anchor` | `g2g` | Naive anchor scan (no minimizer index) |
| `build/abl_brute_matching` | `g2g` | Brute-force matching (no Hungarian) |
| `build/abl_static_partition` | `g2g` | Static partition (no work-stealing) |

### Optional third-party tools

Experiment 3 compares G2G against a PanGraph baseline, and experiment 6
needs minigraph-cactus and gfatools. Install them with:

```bash
bash experiments/download/install_tools.sh
```

PanGraph is fetched as the prebuilt Rust binary (v1.3.0 or later) from the
`neherlab/pangraph` GitHub releases; it is not on bioconda, and the old
Julia package no longer builds. It lands at `$G2G_DATA_ROOT/tools/pangraph`,
and that directory is appended to your PATH in `~/.bashrc`, so refresh your
current shell afterwards (`export PATH="$G2G_DATA_ROOT/tools:$PATH"`, or
open a new one). PanGraph is optional: without it, experiment 3 still runs
and records its PanGraph rows with `exit_code=-1`. That baseline can be
filled in later without rerunning the whole sweep (section 3.1).

## 1. Submitting experiments

The three scripts under `experiments/accl/` are the exact jobs used for the
BIBM 2026 camera-ready re-measurements on our cluster (48-CPU / 192 GB job
limit); they chain the demo harnesses (`demo/0[5-9]*.sh`, `demo/10_*.sh`)
and the exp1b sweep.

Always submit from the repository root: the Slurm scripts set
`G2G_ROOT=${SLURM_SUBMIT_DIR:-$PWD}`, which is the directory you call
`sbatch` from.

```bash
cd /path/to/G2G
conda activate g2g_eval
```

| Exp | Command | Partition | CPUs / RAM / Time | Needs data? |
|---|---|---|---|---|
| 1, correctness | `sbatch experiments/exp1_correctness/slurm.sh` | batch | 16 / 32G / 2h | no (synthetic) |
| 2, minimap2 consistency | `sbatch experiments/exp2_minimap2/slurm.sh` | batch | 16 / 64G / 4h | no (synthetic) |
| 3, throughput | `sbatch experiments/exp3_throughput/slurm.sh` | batch | 64 / 256G / 14h | no (synthetic) |
| 4, parallel scaling | `sbatch experiments/exp4_parallel_scaling/slurm.sh` | batch | 128 / 400G / 12h | no (synthetic) |
| 5, ablation | `sbatch experiments/exp5_ablation/slurm.sh` | batch | 64 / 256G / 12h | no (synthetic) |
| 6, HPRC biology | see section 4 | batch | 64 / 400G / 14h | yes |

Experiments 1 to 5 are fully synthetic and need no downloads, so they can be
submitted right away.

### Wider runs on a large node

For more throughput headroom and a thread sweep past 128 cores:

```bash
sbatch experiments/exp3_throughput/slurm_large_node.sh        # 128 CPU / 512G
sbatch experiments/exp4_parallel_scaling/slurm_large_node.sh  # 256 CPU / 768G
```

## 2. Monitoring jobs

```bash
squeue --me                 # queue / running state
sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,State   # resource summary
tail -f logs/exp1_<jobid>.out
```

Logs land in `logs/`; results in `results/phase5/<exp>/` as timestamped
CSVs.

## 3. Plotting results

### 3.1 Filling in the PanGraph baseline (experiment 3)

If you ran experiment 3 before installing PanGraph (its rows show
`exit_code=-1`), there is no need to rerun the whole sweep. Once `pangraph`
is on the PATH:

```bash
export PATH="$G2G_DATA_ROOT/tools:$PATH"   # if not already
ONLY_PANGRAPH=1 bash experiments/exp3_throughput/run.sh
```

This reuses the GFAs already under `$G2G_DATA_ROOT/synthetic_gfa`, skips
minimap2 and both G2G runs, and writes only the PanGraph rows to
`results/phase5/exp3/exp3_throughput_pangraph_<ts>.csv`. The plot command
below globs `exp3_throughput_*.csv`, so it picks up both the original sweep
and this backfill without a manual merge.

After an experiment finishes, render its figure (PDF and PNG go to
`results/phase5/figures/`):

```bash
python experiments/plot/plot_exp3_throughput.py --csv 'results/phase5/exp3/exp3_throughput_*.csv'
python experiments/plot/plot_exp4_parallel.py   --csv 'results/phase5/exp4/exp4_parallel_scaling_*.csv'
python experiments/plot/plot_exp5_ablation.py \
    --csv-main 'results/phase5/exp5/exp5_ablation_main_*.csv' \
    --csv-a4   'results/phase5/exp5/exp5_ablation_anchor_*.csv' \
    --csv-a5   'results/phase5/exp5/exp5_ablation_matching_*.csv'
python experiments/plot/plot_exp6_biology.py \
    --gnomad-all results/phase5/exp6/g2g_delta_vs_gnomad_all.bed \
    --gnomad-fst results/phase5/exp6/g2g_delta_vs_gnomad_highfst.bed \
    --clinvar    results/phase5/exp6/g2g_shared_vs_clinvar_patho.bed
```

## 4. Experiment 6 (HPRC biological validation)

This one needs real data. By default everything lands under
`$G2G_DATA_ROOT`, which defaults to `$HOME/g2g_data`; if you have a larger
scratch or shared volume, export `G2G_DATA_ROOT=/path/to/share` before
submitting and every stage will use it. Run the stages in order, since each
later stage depends on the previous one.

Compute nodes have no internet, so stages 0 and 1, which download, must run
on the login node and not as Slurm jobs (a batch job would hang or fail);
only the compute-heavy stages 2 and 3 are submitted with `sbatch`.

```bash
# -- On the LOGIN NODE, inside tmux ----------------------

# Stage 0: install third-party tools (PanGraph, minigraph, cactus, gfatools)
bash experiments/download/install_tools.sh

# Stage 1: download datasets (about 400 GB, hours), so use tmux
tmux new -s g2g_dl
conda activate g2g_eval
bash experiments/download/download_all_login.sh
#   HPRC assemblies (balanced AFR/EAS, classified from the 1KGP 3202 panel),
#   GRCh38 reference (no-alt set), gnomAD-SV v4.1, ClinVar (GRCh38)
#   detach with Ctrl-b d ; reattach with `tmux attach -t g2g_dl`
#
#   HPRC Year 1 has only about 6 East Asian assemblies, so the AFR/EAS panels
#   are balanced to the EAS count (about 5 per population) for a size-matched
#   comparison. ClinVar SV extraction is an optional layer; if it yields 0 SVs
#   the experiment 6 job skips it and still produces the core gnomAD-SV result.

# -- On the COMPUTE NODES (as Slurm jobs) --------------------------

# Stage 2: build the AFR and EAS pangenomes with minigraph-cactus
#   (small balanced panel, about 5 samples per population: a few hours each; 14 h walltime)
sbatch experiments/exp6_hprc/slurm_build_pangenomes.sh

# Stage 3: G2G align + gnomAD-SV enrichment + ClinVar overlap
sbatch experiments/exp6_hprc/slurm.sh
```

`download_all_login.sh` runs all three download scripts in sequence and can
be rerun: completed files are skipped, so if a transfer drops you just start
it again. The old `slurm_download.sh` is deprecated and now exits with an
error pointing here.

`slurm.sh` chains the analysis steps: stratify samples, run G2G, intersect
gnomAD-SV, intersect ClinVar. The individual scripts (`01_` through `05_`)
can also be run by hand for debugging.

Compute-node jobs are capped at 14 hours, and every Slurm script here asks
for `--time=14:00:00` or less. The experiment 6 build fits easily because
the balanced panel is small; if a build ever hits the limit, split AFR and
EAS into two separate jobs.

## 5. Data and storage layout

| Path | What |
|---|---|
| `$G2G_DATA_ROOT` (default `$HOME/g2g_data`) | HPRC assemblies, gnomAD/ClinVar, built pangenomes, synthetic GFA scratch |
| `results/phase5/<exp>/` | per-experiment CSVs and timing files (in the repo) |
| `results/phase5/figures/` | rendered figures |
| `logs/` | Slurm stdout/stderr |

Datasets and pangenome builds go under your home by default; if you have a
larger scratch or shared volume, export `G2G_DATA_ROOT` to point there
before submitting any job. The synthetic experiments (3 to 5) only use it as
scratch for generated GFAs, while experiment 6 stores the real HPRC, gnomAD
and ClinVar data and the built pangenomes there.

## 6. Ablation flags (experiment 5)

The five variants are compile-time CMake options (see the top-level
`CMakeLists.txt`), each disabling one optimization:

| Flag | Disables |
|---|---|
| `G2G_BFS_LAYOUT` | cache-oblivious van Emde Boas snarl-tree layout |
| `G2G_STEAL_RANDOM` | steal-largest-first work-stealing policy |
| `G2G_NAIVE_ANCHOR` | minimizer anchor index |
| `G2G_BRUTE_MATCHING` | Hungarian (LAPJV) bipartite matching |
| `G2G_STATIC_PARTITION` | work-stealing scheduler |

`01_build_all.sh` builds all of them into separate `build/abl_*` directories,
and `exp5_ablation/run.sh` drives the comparison.
