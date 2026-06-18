# G2G Phase 5 — Evaluation

Experiment suite for the Graph-to-Graph (G2G) pangenome aligner. Six
experiments cover correctness, baseline consistency, throughput, parallel
scaling, an ablation study, and biological validation on HPRC data.

Everything here is tuned for the **KAUST ROCS / SANDS testbed** (Slurm head
node `mcmgt01`, partitions `mc` and `synode`, conda/mamba, no sudo). All jobs
record wall time and peak memory.

---

## 0. One-time setup

Run these once from the **repo root**, on the head node `mcmgt01`.

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

---

## 1. Submitting experiments

**Always submit from the repo root** — the SLURM scripts set
`G2G_ROOT=${SLURM_SUBMIT_DIR:-$PWD}`, i.e. the directory you call `sbatch`
from.

```bash
cd ~/Research/G2G/G2G_Project_Summer_2026
conda activate g2g_eval
```

| Exp | Command | Partition | CPUs / RAM / Time | Needs data? |
|---|---|---|---|---|
| 1 — Correctness | `sbatch experiments/exp1_correctness/slurm.sh` | mc | 16 / 32G / 2h | no (synthetic) |
| 2 — minimap2 consistency | `sbatch experiments/exp2_minimap2/slurm.sh` | mc | 16 / 64G / 4h | no (synthetic) |
| 3 — Throughput | `sbatch experiments/exp3_throughput/slurm.sh` | mc | 64 / 256G / 24h | no (synthetic) |
| 4 — Parallel scaling | `sbatch experiments/exp4_parallel_scaling/slurm.sh` | mc | 128 / 400G / 12h | no (synthetic) |
| 5 — Ablation | `sbatch experiments/exp5_ablation/slurm.sh` | mc | 64 / 256G / 12h | no (synthetic) |
| 6 — HPRC biology | see [§4](#4-experiment-6-hprc-biological-validation) | mc | 64 / 400G / 48h | **yes** |

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
jobstats <jobid>            # resource-utilization summary (CPU/mem/peak)
jobstats <jobid> -g         # Grafana dashboard URL (login rocs/rocs)
tail -f logs/exp1_<jobid>.out
```

Logs land in `logs/`; results in `results/phase5/<exp>/` as timestamped CSVs.

---

## 3. Plotting results

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

This one needs real data on the shared project volume `/data/fat/g2g`. Run the
stages in order; each later stage depends on the previous.

```bash
# Stage 0 — install third-party tools (Julia/PanGraph, minigraph, cactus, gfatools)
bash experiments/download/install_tools.sh

# Stage 1 — download datasets (~hours; large). Runs as a Slurm job.
sbatch experiments/exp6_hprc/slurm_download.sh
#   → HPRC assemblies (20 AFR + 20 EAS), gnomAD-SV v4.1, ClinVar (GRCh38)

# Stage 2 — build the AFR and EAS pangenomes with minigraph-cactus
#   (~48 CPU-h each; 120h walltime ceiling)
sbatch experiments/exp6_hprc/slurm_build_pangenomes.sh

# Stage 3 — G2G align + gnomAD-SV enrichment + ClinVar overlap
sbatch experiments/exp6_hprc/slurm.sh
```

`slurm.sh` chains the analysis steps: stratify samples → run G2G → intersect
gnomAD-SV → intersect ClinVar. The individual scripts
(`01_…` through `05_…`) can also be run by hand for debugging.

> **Tip:** `slurm_build_pangenomes.sh` is 5 days of walltime. A future cluster
> policy will require jobs >3 days to be preemptible (`--qos=spot`); it's not
> enforced yet, but if it starts failing to schedule, add `#SBATCH --qos=spot`.

---

## 5. Data & storage layout

| Path | What |
|---|---|
| `/data/fat/g2g` | `$G2G_DATA_ROOT` — HPRC assemblies, gnomAD/ClinVar, built pangenomes |
| `results/phase5/<exp>/` | Per-experiment CSVs + timing files (in repo, on `/home`) |
| `results/phase5/figures/` | Rendered figures |
| `logs/` | Slurm stdout/stderr |

The `/data/fat/` share is the multi-user project volume (65 TB, RAIDZ2); keep
large datasets and pangenome builds there, not in `/home`.

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
