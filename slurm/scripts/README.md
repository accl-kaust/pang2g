# G2G — exp3 & exp6 Fix + Setup Guide

This guide makes **exp3** (throughput scaling) and **exp6** (biological validation)
run correctly on the SANDS/ROCS SLURM testbed. It covers what was broken, the
scripts that fix it, the exact order to stage data and submit jobs, and the
copy-paste command sequence.

> **One environment variable governs everything below.** Your real data lives at
> `~/Research/G2G/hprc_data`, but the scripts default to `$REPO/data/hprc`. So in
> **every** shell you use for staging or submitting, export:
>
> ```bash
> export G2G_DATA_ROOT=~/Research/G2G/hprc_data
> ```
>
> (Add it to `~/.bashrc` to avoid forgetting.) Skip this and jobs will look in the
> wrong directory and report missing data.

---

## 0. TL;DR — the whole sequence

Run on the **login node** (`mcmgt01`), from the repo root
`~/Research/G2G/G2G_Project_Summer_2026`:

```bash
export G2G_DATA_ROOT=~/Research/G2G/hprc_data

# --- 1. Stop the runaway exp3 and its blocked dependents ---
scancel 33651 33652 33653

# --- 2. Stage data (login node has internet; compute nodes do NOT) ---
bash slurm/scripts/prestage_hprc.sh          # 60 HPRC assemblies (already done for you)
bash slurm/scripts/prestage_refs.sh          # GRCh38 + gnomAD-SV + ClinVar truth sets

# --- 3. Install the graph-build tools into the phase5 env (one time) ---
conda install -n phase5 -c bioconda -c conda-forge cactus minigraph vg

# --- 4. Resubmit, patched. exp3 with a sane first-pass config: ---
sbatch --export=ALL,G2G_REPO=$PWD,G2G_DATA_ROOT=$G2G_DATA_ROOT,G2G_SIZES="1000 10000 100000 1000000",G2G_RUNS=3 \
       slurm/scripts/exp3_throughput.sh

# --- 5. Resubmit exp6 (reference data now present): ---
sbatch --export=ALL,G2G_REPO=$PWD,G2G_DATA_ROOT=$G2G_DATA_ROOT \
       slurm/scripts/exp6_biological.sh

squeue --me
```

Everything past this point explains *why*, and how to scale up.

---

## 1. What was broken, and what fixes it

### exp3 — ran 21+ h and would never finish
The single-threaded baseline aligner is **O(n_A·n_B)** — quadratic in snarl count.
Measured: n=10⁴ takes ~80 min/run. Extrapolated: n=10⁵ ≈ 5.6 days/run, n=10⁶ ≈
1.5 years/run. The default size sweep goes to **10⁷**, so the sequential baseline
never completes, monopolising the only `mc` node and blocking exp4/exp5 (which wait
on a dependency).

**Fix (`exp3_throughput.sh`):** a sequential-baseline size gate, `G2G_SEQ_MAX_N`
(default **10000**). The sequential baseline runs **only** at n ≤ this threshold —
where it's needed as the speedup *denominator* and is tractable. The **parallel**
runs still execute at **every** size, so the headline throughput-scaling curve is
preserved up to 10⁷. Above the gate you get one clear log line and a
`skipped_seq_max_n` row instead of a multi-day hang.

### exp6 — four separate blockers
1. **Wrong HPRC download URLs** → fixed (correct S3 layout, driven by the official
   HPRC index). *(done in earlier session)*
2. **Compute nodes are air-gapped** → downloads moved to the login node via
   `prestage_hprc.sh`; the compute job only *verifies* pre-staged files.
3. **SLURM path resolution** — submitting directly left `$G2G_REPO` unset, so the
   job tried to write to `/var/lib/slurm` (permission denied). Fixed: REPO resolves
   from `$G2G_REPO` → `$SLURM_SUBMIT_DIR` → script location, with a loud guard. The
   reliable invocation is `--export=ALL,G2G_REPO=$PWD`.
4. **cactus reference misconfiguration** — the build used `--reference _MINIGRAPH_`,
   but that's cactus's *internal* label for minigraph nodes, **not** a valid
   reference sample. cactus needs a real genome (GRCh38) as the uncollapsed
   coordinate system. Fixed: Stage 1 writes a `GRCh38` row in the seqfile, Stage 2
   passes `--reference GRCh38`, plus a preflight guard that fails with a clear
   message if GRCh38 is missing.

**Remaining requirement:** the three reference/truth files must be staged (next
section). Without them Stage 2 (graph build) and Stage 5 (enrichment) cannot run,
and exp6 exits 3.

---

## 2. Staging the data (login node only)

The testbed compute nodes have **no outbound internet**; the login node does. Both
staging scripts therefore run on `mcmgt01`, are **idempotent** (skip files already
present and verified), and fail fast with a clear message if run somewhere
air-gapped.

### 2a. HPRC assemblies — `prestage_hprc.sh`
Downloads the 60 HPRC assemblies (30 AFR + 30 EAS haplotypes), driven by the
official HPRC assembly index, SHA-256 verifying each file. Writes to
`$G2G_DATA_ROOT/assemblies/`. *(You've already run this successfully.)*

> **Don't delete `$G2G_DATA_ROOT/index/`.** `prestage_hprc.sh` also stores a
> small (33K) HPRC index file at `$G2G_DATA_ROOT/index/Year1_assemblies_v2_genbank.index`.
> exp6 needs it on the **compute node** to map each sample to its URL+SHA — and
> the compute node can't re-fetch it (air-gapped). A cleanup that empties `index/`
> while leaving `assemblies/` intact will make exp6 fail Stage 1 with
> `HPRC index NOT available`, even though all 60 assemblies are present. If that
> happens, re-fetch just the index (see the troubleshooting table) or re-run
> `prestage_hprc.sh`.

```bash
export G2G_DATA_ROOT=~/Research/G2G/hprc_data
bash slurm/scripts/prestage_hprc.sh
```

### 2b. Reference + truth sets — `prestage_refs.sh`
Downloads and converts the three files exp6 needs into the **exact** BED layouts
Stage 5 expects. Writes to `$G2G_DATA_ROOT/reference/`:

| File | Source | Format produced |
|---|---|---|
| `GRCh38.fa` | NCBI no-alt analysis set (HPRC/1000G lineage) | FASTA (+ `.fai` if samtools present) |
| `gnomad_sv_v4.bed` | gnomAD v4.1 SV sites (GCS public bucket) | `chrom  start  end  id  fst` |
| `clinvar_sv.bed` | dbVar **nstd102** (ClinVar SVs), GRCh38 | `chrom  start  end  id` |

```bash
export G2G_DATA_ROOT=~/Research/G2G/hprc_data
bash slurm/scripts/prestage_refs.sh
```

**About the `fst` column (read this).** gnomAD-SV v4 does not publish a per-pair
F_ST, so the script derives an **AFR-vs-EAS differentiation proxy** from the
per-population allele frequencies it does publish:

```
fst = |AF_afr − AF_eas|      (0..1; higher = more differentiated)
```

This is a transparent stand-in that matches the AFR/EAS panels exp6 compares — it
is **not** true Weir–Cockerham F_ST. State it that way in the paper. If you obtain a
real F_ST track, point `G2G_GNOMAD_URL` at it (or replace `gnomad_sv_v4.bed`
directly, keeping column 5 = F_ST). SVs shorter than 50 bp are filtered out
(`G2G_MIN_SV_LEN`).

**If you must stage from your laptop instead** (e.g. a future node has no
internet), run `prestage_refs.sh` anywhere with internet, then copy the directory
over:

```bash
rsync -avP ~/Research/G2G/hprc_data/reference/  shaklea@mcmgt01:~/Research/G2G/hprc_data/reference/
```

---

## 3. Installing the graph-build tools (one time)

exp6 Stage 2 builds pangenomes with **cactus** (and uses **vg**/**minigraph**).
These confirmed present in your `phase5` env already (cactus 3.1.4, minigraph 0.21,
vg 1.73.0). If you rebuild the env:

```bash
conda install -n phase5 -c bioconda -c conda-forge cactus minigraph vg
```

---

## 4. Submitting the jobs

### Always pass `G2G_REPO=$PWD`
Submit from the repo root and **always** include `--export=ALL,G2G_REPO=$PWD`. This
is what fixed the `/var/lib/slurm` permission error. (`SLURM_SUBMIT_DIR` now also
resolves it as a fallback, but passing `G2G_REPO` explicitly is the reliable path.)

### exp3 — choose a config to match your time budget
The patched script keeps the **full parallel curve** and gates only the sequential
baseline. Runtime is dominated by the largest parallel sizes.

**Recommended first pass** (completes in a few hours, full scaling story up to 10⁶):
```bash
sbatch --export=ALL,G2G_REPO=$PWD,G2G_DATA_ROOT=$G2G_DATA_ROOT,G2G_SIZES="1000 10000 100000 1000000",G2G_RUNS=3 \
       slurm/scripts/exp3_throughput.sh
```

**Full sweep to 10⁷** (≈24 h — bumps the `--time=24:00:00` wall; raise `--time` or
drop `G2G_RUNS` to 2 if it's tight):
```bash
sbatch --export=ALL,G2G_REPO=$PWD,G2G_DATA_ROOT=$G2G_DATA_ROOT \
       slurm/scripts/exp3_throughput.sh
```

Knobs: `G2G_SIZES` (sweep), `G2G_RUNS` (repeats/point, default 3),
`G2G_SEQ_MAX_N` (sequential cap, default 10000 — raising it adds days, don't unless
you have a reason).

### exp6 — once reference data is staged
```bash
sbatch --export=ALL,G2G_REPO=$PWD,G2G_DATA_ROOT=$G2G_DATA_ROOT \
       slurm/scripts/exp6_biological.sh
```

**Heads-up on cost:** the cactus build (Stage 2) is the long pole — roughly **48
CPU-hours per panel** (AFR and EAS), so expect a long wall time on a 128-core node.
It's checkpointed by panel: a finished panel is skipped on resubmit. Watch progress:
```bash
tail -f slurm/logs/exp6-<JOBID>.err
```
Stage status `S1..S5` is printed near the end; all five = a complete run.

---

## 5. Cancelling the stuck jobs

The runaway exp3 and its two dependents:
```bash
scancel 33651          # runaway exp3
scancel 33652 33653    # exp4, exp5 (waiting on exp3's dependency)
```
Confirm with `squeue --me`. Then resubmit per §4 (and exp4/exp5 via your usual
`submit_all.sh`, which will re-create the dependency chain against the new exp3).

---

## 6. Housekeeping

Stale logs from the failed exp6 attempts can be removed so you tail the right one:
```bash
rm -f slurm/logs/exp6-33690* slurm/logs/exp6-33706* slurm/logs/exp6-33707*
```
Raw reference downloads can be deleted after conversion to reclaim space:
```bash
rm -rf ~/Research/G2G/hprc_data/reference/_downloads
```

---

## 7. Quick troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `mkdir: cannot create '/var/lib/slurm/...'` | `G2G_REPO` unset | resubmit with `--export=ALL,G2G_REPO=$PWD` |
| Stage 1 `missing_prestage` | assemblies not staged / wrong `G2G_DATA_ROOT` | export `G2G_DATA_ROOT`, run `prestage_hprc.sh` |
| Stage 1 `HPRC index NOT available` | `index/` dir emptied (e.g. by a cleanup) while assemblies remain | re-fetch the 33K index: `wget -O $G2G_DATA_ROOT/index/Year1_assemblies_v2_genbank.index https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index` (or re-run `prestage_hprc.sh`) |
| `PREFLIGHT FAILED` with a numbered list | one or more data inputs missing | do the one command shown next to each numbered item, then resubmit |
| Stage 2 `GRCh38 reference FASTA missing` | truth data not staged | run `prestage_refs.sh` |
| Stage 2 `neither cactus nor minigraph` | tools not in env | `conda install -n phase5 -c bioconda -c conda-forge cactus minigraph vg` |
| exp3 still runs for days | `G2G_SEQ_MAX_N` raised, or huge `G2G_SIZES` | keep `SEQ_MAX_N=10000`; trim `G2G_SIZES` |
| `no outbound internet from this machine` | staging on a compute node | run staging on `mcmgt01`, or `rsync` from a machine with internet |

---

## Files in this delivery
- `exp3_throughput.sh` — patched (sequential size gate + REPO resolution).
- `exp6_biological.sh` — patched (index-verify Stage 1, GRCh38 reference, Stage 2 guards, REPO resolution).
- `prestage_hprc.sh` — login-node HPRC assembly downloader (SHA-verified, resumable).
- `prestage_refs.sh` — login-node GRCh38 + gnomAD-SV + ClinVar downloader/converter.
- `README.md` — this guide.
