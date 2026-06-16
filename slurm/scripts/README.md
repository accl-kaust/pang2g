# G2G Phase 5 — Evaluation Campaign (SLURM)

Publication-grade experiment scripts for the G2G graph-to-graph pangenome
alignment paper (target venue: **IEEE BIBM**). This directory drives the entire
Phase 5 evaluation: throughput, parallel scaling, ablations, and biological
validation, with full logging, wall-clock + peak-memory measurement, per-phase
timing, and reproducibility manifests.

> **Scale.** These scripts default to **HPRC-realistic** graph sizes — up to
> ~10M snarls (~36M nodes, ~3.3 GB per GFA), matching real human pangenomes.
> The whole synthetic sweep (1K…10M, ×2 graphs, ×3 seeds) is ~25 GB on disk and
> fits the 512 GB / 128-core / >30 TB testbed. Drop the top size if RAM-limited
> (see `G2G_N_TREATMENT` / `G2G_SIZES` below).

---

## What's here

| File | Purpose |
|------|---------|
| `lib_measure.sh` | Shared bash measurement library. Wall-clock + peak-RSS + CPU% via GNU `time -v`, optional `perf` L3-miss rate, environment capture, SHA manifests. Sourced by every experiment. |
| `setup_cluster.sh` | One-time environment setup: conda env `phase5`, **OpenMP** toolchain, Release build (`-O3 -march=native`), full `ctest` gate, perf/large-graph readiness checks. |
| `exp3_throughput.sh` | **Q2 + Q3.** Throughput sweep 1K…10M snarls; G2G sequential + parallel + baseline panel (vg giraffe/map, GraphAligner, minigraph, minimap2 via linearise-and-map). Includes the haplotype-panel sweep (the headline figure) and a baseline failure-mode log. |
| `exp4_parallel_scaling.sh` | **Q2.** Thread sweep 1→128 at multiple graph sizes. Reports **two** speedup curves — end-to-end and DP-phase-only — plus an Amdahl sequential-fraction analysis. |
| `exp5_ablation.sh` | **Q4.** Five ablations (vEB vs BFS layout, steal-largest vs random ×2, minimizer-index vs naïve, Hungarian vs brute-force) run at both a control size (fits L3) and HPRC scale (exceeds L3). |
| `exp6_biological.sh` | **Q5.** HPRC download → minigraph-cactus pangenomes (AFR vs EAS) → G2G alignment → delta-allele intersection with gnomAD-SV + ClinVar → permutation enrichment + vg-giraffe comparison. |
| `submit_all.sh` | Master submission. Submits all four jobs with correct dependencies (exp4/exp5 wait on exp3's graph cache; exp6 runs independently). Supports `--dry-run`, `--local`, `--no-bio`, and named-subset submission. |
| `panels/afr_samples.txt`, `panels/eas_samples.txt` | HPRC sample IDs (15 each → 30 haplotypes) for the two population panels in exp6. |
| `../synth/gen_hprc_scale_gfa.py` | HPRC-calibrated synthetic GFA generator (size mode + haplotype-panel mode; Zipf allele distribution fit to HPRC snarl statistics). Emits ground-truth TSVs. |

---

## Quick start

```bash
# 0. One-time setup (build + env + sanity checks)
bash slurm/scripts/setup_cluster.sh

# 1. Launch the whole campaign
bash slurm/scripts/submit_all.sh

#    …or preview without submitting
bash slurm/scripts/submit_all.sh --dry-run

#    …or run a single experiment
sbatch slurm/scripts/exp3_throughput.sh
```

Results land under `results/{exp3,exp4,exp5,exp6}/`, each with its own
`env.txt` (full hardware/software capture) and `manifest.sha256`.

---

## Run order & dependencies

```
exp6 ───────────────────────────────────────────►   independent, long pole
exp3 ──┬─► exp4                                      afterok (shares graph cache)
       └─► exp5                                      afterok
```

* **exp6 first.** HPRC download (~1 TB) and two minigraph-cactus builds
  (~96 CPU-h total) are the long pole. It shares nothing with the synthetic
  jobs, so launch it first and let it run in the background. Every stage is
  **resumable** — re-running after pre-emption re-uses cached downloads and
  built graphs.
* **exp3 before exp4/exp5.** exp3 generates the synthetic graph cache in
  `$G2G_GRAPHDIR`; exp4 and exp5 reuse it. `submit_all.sh` wires the
  `afterok` dependency automatically.

---

## Configuration (environment variables)

All paths and sweep parameters are env-configurable; defaults assume the repo
layout. Override per-run or export once before `submit_all.sh`.

**Paths**
| Var | Default | Meaning |
|-----|---------|---------|
| `G2G_REPO` | auto (two levels up) | Repository root |
| `G2G_BIN` | `$REPO/build/g2g` | G2G binary |
| `G2G_ENV` | `phase5` | conda env name |
| `CONDA_BASE` | auto-detect | conda install prefix |
| `G2G_GRAPHDIR` | `$REPO/results/exp3/graphs` | synthetic graph cache (shared) |
| `G2G_DATA_ROOT` | `$REPO/data/hprc` | HPRC assemblies + pangenomes (point at scratch!) |
| `G2G_SCRATCH` | `$DATA_ROOT/work` | cactus scratch (huge, ephemeral) |

**Sweep sizing**
| Var | Default | Used by |
|-----|---------|---------|
| `G2G_SIZES` | `1000 10000 100000 1000000 10000000` | exp3 throughput points |
| `G2G_RUNS` | `5` | repetitions per measurement |
| `G2G_THREADS` | `1 2 4 8 16 32 64 128` | exp4 thread sweep |
| `G2G_N_TREATMENT` | `10000000` | exp5 HPRC-scale ablation size |
| `G2G_N_CONTROL` | `10000` | exp5 control size (fits L3) |

**exp6 truth sets** (supply these or stage them under `$REF_DIR`)
| Var | Meaning |
|-----|---------|
| `G2G_GRCH38` | GRCh38 FASTA (for coordinate lifting + genome sizes) |
| `G2G_GNOMAD_SV` | gnomAD-SV v4 BED (`chrom start end id fst`) |
| `G2G_CLINVAR_SV` | ClinVar SV BED (`chrom start end id`) |
| `G2G_PERM_N` | permutations for the enrichment null (default 1000) |

---

## Output instrumentation (important)

The scaling/ablation scripts scrape two stderr lines from the G2G binary:

```
[timing]   io=<s> index=<s> decompose=<s> dp=<s>
[runstats] tasks=<n> steals=<n> imbalance=<r>
```

These drive the per-phase columns and the Amdahl analysis. **If `main.cpp` /
`parallel.cpp` do not yet emit them**, the scripts degrade gracefully (phase
columns become `NA`) — but you lose the decomposition that explains *why* the
parallel speedup behaves as it does. Adding these two `fprintf(stderr, …)`
lines behind a `--verbose-timing` flag is the single highest-value
instrumentation change for the paper. (Listed in the roadmap gap analysis.)

---

## Notes on measurement fidelity

* Wall-clock and peak RSS come from GNU `time -v`, **not** the bash builtin —
  the builtin can't report memory and rounds to 1 ms.
* "Memory used" in every table is **Maximum resident set size** (the physical
  high-water mark reviewers expect), reported in GB in the summaries.
* CPU-seconds (user+sys) are recorded alongside wall-clock so parallel
  efficiency is reasoned about honestly, never conflated.
* `perf` L3-miss measurement degrades to `NA` when `perf_event_paranoid`
  restricts access — check `results/*/env.txt` to see whether it was available
  on the node.

---

## Toolchain

The runtime is **OpenMP** (migrated from oneTBB). `setup_cluster.sh` installs no
TBB package; it only needs a C++17 compiler with OpenMP, cmake, and the
bioinformatics tools for exp6 (`minigraph`, `cactus-pangenome`, `bedtools`,
`vg`, `samtools`). Missing exp6 tools are detected at runtime and the affected
stages skip rather than abort.
