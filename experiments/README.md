# G2G Phase 5 — Evaluation

Experiment suite for the Graph-to-Graph (G2G) pangenome aligner. Six
experiments cover correctness, baseline consistency, throughput, parallel
scaling, an ablation study, and biological validation on HPRC data.

Everything here runs on a Slurm cluster with conda/mamba and no sudo (it was
developed on a 64-core AMD EPYC 7763 node; set `--partition`/`--account` in the
`slurm.sh` files for your site, e.g. `batch` on KAUST Ibex, where CPU jobs are
submitted from `ilogin`). All jobs record wall time and peak memory.

For a first run on IBEX use the top-level test case instead of this suite:
`sbatch experiments/ibex/test_case.slurm` (see the root `README.md`); the
other scripts in `experiments/ibex/` are the whole-node jobs behind the
paper's IBEX measurements and need no conda environment, only the modules
named in `experiments/ibex/common.sh`.

---

## 0. One-time setup

Run these once from the **repo root**, on a login node (compute nodes usually
have no outbound network).

```bash
# 1. Create the conda env (cmake, compilers, ninja, minimap2, bedtools,
#    samtools, bcftools, seqkit, k8, plus Python analysis stack)
bash experiments/00_setup_env.sh
conda activate g2g_eval

# 2. Build the production binary + 5 ablation variants
bash experiments/01_build_all.sh
```

`01_build_all.sh` auto-activates `g2g_eval` and preflights the toolchain, so
it also works if you forget to activate first. Outputs:

| Build dir | Binary | Purpose |
|---|---|---|
| `build/release` | `g2g` | Production (all optimizations) |
| `build/abl_bfs_layout` | `g2g` | BFS layout (no van Emde Boas) |
| `build/abl_steal_random` | `g2g` | Random steal (no largest-first) |
| `build/abl_naive_anchor` | `g2g` | Naïve anchor scan (no minimizer index) |
| `build/abl_brute_matching` | `g2g` | Brute-force matching (no Hungarian) |
| `build/abl_static_partition` | `g2g` | Static partition (no work-stealing) |

### Optional third-party tools

Exp3 compares G2G against a **PanGraph** baseline, and Exp6 needs
minigraph-cactus / gfatools. Install them with:

```bash
bash experiments/download/install_tools.sh
```

- **PanGraph** is fetched as the prebuilt Rust binary (v1.3.0+) from the
  `neherlab/pangraph` GitHub releases — it is *not* on bioconda, and the old
  Julia package no longer builds. It lands at `$G2G_DATA_ROOT/tools/pangraph`
  and that dir is appended to your PATH in `~/.bashrc`.
- After running the installer, refresh your current shell so `pangraph`
  resolves: `export PATH="$G2G_DATA_ROOT/tools:$PATH"` (or open a new shell).
- PanGraph is **optional**: if it isn't installed, Exp3 still runs and records
  its PanGraph rows with `exit_code=-1` (skipped sentinel). You can backfill
  that baseline later without rerunning the whole sweep — see
  [§3.1](#31-backfilling-the-pangraph-baseline-exp3).

---

## 1. Submitting experiments

> The three scripts under `experiments/accl/` are the exact jobs used for the BIBM 2026
> camera-ready re-measurements on the authors' cluster (48-CPU / 192 GB job limit); they
> chain the demo harnesses (`demo/0[5-9]*.sh`, `demo/10_*.sh`) and the exp1b sweep.


**Always submit from the repo root** — the SLURM scripts set
`G2G_ROOT=${SLURM_SUBMIT_DIR:-$PWD}`, i.e. the directory you call `sbatch`
from.

```bash
cd /path/to/G2G
conda activate g2g_eval
```

| Exp | Command | Partition | CPUs / RAM / Time | Needs data? |
|---|---|---|---|---|
| 1 — Correctness | `sbatch experiments/exp1_correctness/slurm.sh` | batch | 16 / 32G / 2h | no (synthetic) |
| 2 — minimap2 consistency | `sbatch experiments/exp2_minimap2/slurm.sh` | batch | 16 / 64G / 4h | no (synthetic) |
| 3 — Throughput | `sbatch experiments/exp3_throughput/slurm.sh` | batch | 64 / 256G / 14h | no (synthetic) |
| 4 — Parallel scaling | `sbatch experiments/exp4_parallel_scaling/slurm.sh` | batch | 128 / 400G / 12h | no (synthetic) |
| 5 — Ablation | `sbatch experiments/exp5_ablation/slurm.sh` | batch | 64 / 256G / 12h | no (synthetic) |
| 6 — HPRC biology | see [§4](#4-experiment-6-hprc-biological-validation) | batch | 64 / 400G / 14h | **yes** |

Experiments 1–5 are fully synthetic and need no downloads — submit them right
away.

### Wider runs on `synode` (384 CPU / 1.5 TB)

For more throughput headroom and a thread sweep past 128 cores:

```bash
sbatch experiments/exp3_throughput/slurm_synode.sh        # 128 CPU / 512G
sbatch experiments/exp4_parallel_scaling/slurm_synode.sh  # 256 CPU / 768G
```

---

## 2. Monitoring jobs

```bash
squeue --me                 # queue / running state
sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,State   # resource summary
tail -f logs/exp1_<jobid>.out
```

Logs land in `logs/`; results in `results/phase5/<exp>/` as timestamped CSVs.

---

## 3. Plotting results

### 3.1 Backfilling the PanGraph baseline (Exp3)

If you ran Exp3 before installing PanGraph (its rows show `exit_code=-1`),
you don't need to rerun the whole sweep. Once `pangraph` is on PATH:

```bash
export PATH="$G2G_DATA_ROOT/tools:$PATH"   # if not already
ONLY_PANGRAPH=1 bash experiments/exp3_throughput/run.sh
```

This reuses the GFAs already under `$G2G_DATA_ROOT/synthetic_gfa`, skips
minimap2 and both G2G runs, and writes only the PanGraph rows to
`results/phase5/exp3/exp3_throughput_pangraph_<ts>.csv`. The plot command
below globs `exp3_throughput_*.csv`, so it picks up both the original sweep
and this backfill automatically — no manual merge needed.

After an experiment finishes, render its figure (writes PDF+PNG to
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

---

## 4. Experiment 6 (HPRC biological validation)

This one needs real data. By default everything lands under
`$G2G_DATA_ROOT` (defaults to `$HOME/g2g_data`, since the shared
`/data/fat` volume is not writable by individual user accounts). If you
do get write access to a larger shared volume, just export
`G2G_DATA_ROOT=/path/to/share` before submitting and every stage will
use it. Run the stages in order; each later stage depends on the previous.

> **Compute nodes have no internet.** Stages 0 and 1 (anything that
> downloads) must run on the **login node**, *not* as a Slurm job —
> a batch job has no outbound network and will hang/fail. Only the
> compute-heavy stages (2 and 3) are submitted with `sbatch`.

```bash
# ── On the LOGIN NODE, inside tmux ──────────────────────

# Stage 0 — install third-party tools (PanGraph, minigraph, cactus, gfatools)
bash experiments/download/install_tools.sh

# Stage 1 — download datasets (~400 GB, hours). Long-running → use tmux.
tmux new -s g2g_dl
conda activate g2g_eval
bash experiments/download/download_all_login.sh
#   → HPRC assemblies (balanced AFR/EAS, classified from the 1KGP 3202 panel),
#     GRCh38 reference (no-alt set), gnomAD-SV v4.1, ClinVar (GRCh38)
#   detach with Ctrl-b d ; reattach with `tmux attach -t g2g_dl`
#
#   NB: HPRC Year-1 has only ~6 East-Asian assemblies, so the AFR/EAS panels
#   are auto-balanced to the EAS count (≈5 per pop) for a size-matched compare.
#   ClinVar SV extraction is an optional layer; if it yields 0 SVs the Exp6
#   Slurm job skips it and still produces the core gnomAD-SV result.

# ── On the COMPUTE NODES (as Slurm jobs) ──────────────────────────

# Stage 2 — build the AFR and EAS pangenomes with minigraph-cactus
#   (small balanced panel ~5 samples/pop → a few hours each; 14h walltime)
sbatch experiments/exp6_hprc/slurm_build_pangenomes.sh

# Stage 3 — G2G align + gnomAD-SV enrichment + ClinVar overlap
sbatch experiments/exp6_hprc/slurm.sh
```

`download_all_login.sh` runs all three download scripts in sequence and is
re-runnable — completed files are skipped, so if a transfer drops you can
just start it again. The old `slurm_download.sh` is deprecated and now exits
with an error pointing here.

`slurm.sh` chains the analysis steps: stratify samples → run G2G → intersect
gnomAD-SV → intersect ClinVar. The individual scripts
(`01_…` through `05_…`) can also be run by hand for debugging.

> **Walltime cap:** compute-node jobs are limited to **14 hours**. Every Slurm
> script here is set to `--time=14:00:00` or less. The Exp6 build fits easily
> because the balanced panel is small (~5 samples/pop). If a build ever bumps
> the limit, split AFR and EAS into two separate jobs.

---

## 5. Data & storage layout

| Path | What |
|---|---|
| `$G2G_DATA_ROOT` (default `$HOME/g2g_data`) | HPRC assemblies, gnomAD/ClinVar, built pangenomes, synthetic GFA scratch |
| `results/phase5/<exp>/` | Per-experiment CSVs + timing files (in repo, on `/home`) |
| `results/phase5/figures/` | Rendered figures |
| `logs/` | Slurm stdout/stderr |

`G2G_DATA_ROOT` defaults to `$HOME/g2g_data`. The cluster's `/data/fat`
share is not writable by individual user accounts, so datasets and
pangenome builds go under your home by default. If you are granted a quota
on a larger shared volume, export `G2G_DATA_ROOT` to point there before
submitting any job — the synthetic experiments (3–5) only use it as scratch
for generated GFAs, while Exp6 stores the real HPRC/gnomAD/ClinVar data and
built pangenomes there.

---

## 6. Ablation flags (Experiment 5)

The five variants are compile-time CMake options (see top-level
`CMakeLists.txt`), each disabling one optimization:

| Flag | Disables |
|---|---|
| `G2G_BFS_LAYOUT` | cache-oblivious van Emde Boas snarl-tree layout |
| `G2G_STEAL_RANDOM` | steal-largest-first work-stealing policy |
| `G2G_NAIVE_ANCHOR` | minimizer anchor index |
| `G2G_BRUTE_MATCHING` | Hungarian (LAPJV) bipartite matching |
| `G2G_STATIC_PARTITION` | work-stealing scheduler |

`01_build_all.sh` builds all of them into separate `build/abl_*` dirs;
`exp5_ablation/run.sh` drives the comparison.
