# BIBM 2026 camera-ready TODO — status after the 2026-10-04 revision pass

Branch: `claude/camera-ready-paper-e5f755` (worktree; not merged to main).
R1/R2/R3 = reviewer #1/#2/#3. Camera-ready deadline: **25 Oct 2026**; limit **8 pages
including references** (BIBM 2026 CFP). The compiled draft is `latex/PanG2G_camera_ready_draft.pdf`
(8 pages). Items marked **USER** need an action only you can take.

## 1. Code and data release (R1.3, R2.1)

- [ ] **USER** Make `github.com/ahmadshakleya/G2G` public (still 404 on 2026-10-04).
- [x] Secrets / cluster paths scrubbed: personal e-mail removed from all `#SBATCH --mail-user` lines,
      `--partition=mc/synode` replaced by a generic `batch` with a comment, `mcmgt01`/ROCS/SANDS references
      removed from `experiments/README.md` and `download_all_login.sh`. No credentials were found.
- [x] `LICENSE` (MIT), root `README.md` (build, run, reproduction table, data versions), `CITATION.cff`,
      `.zenodo.json` added.
- [ ] **USER** Mint the Zenodo DOI for the tagged release and replace `ZENODO-DOI-PLACEHOLDER` in
      `Sections/5_Conclusion.tex` (and `README.md`).
- [x] Availability claim verified and corrected:
  - [x] The AFR/EAS graphs are built from the **1000 Genomes phase-3 panel + vg construct**
        (`demo/09_fetch_biology_data.sh`, `demo/09_hprc_biology.sh`); the paper now says so and names the
        releases. `experiments/exp6_hprc/` (HPRC + Minigraph-Cactus) is described as scripts only.
  - [x] Reference-anchored sweep scripts (`demo/05_scaling_demo.sh`, `demo/06_scaling_limit.sh`) and the
        F_ST sweep (`demo/hprc_enrich.py --fst-sweep`) are included and documented.
  - [x] `make_figures.py` (root) regenerates **every** figure in `latex/figs/` from committed tables
        (`results/*.tsv|csv`, `demo/sample_outputs/*.tsv`); `demo/make_figures.py` is a shim to it.
- [x] Data versions in README and paper: 1000G release 20130502, GRCh37, gnomAD-SV v2.1,
      ClinVar GRCh37 weekly release 2026-09-28, PGGB HLA-DRB1 haplotypes, tool versions.

## 2. Code/paper mismatches

- [x] **Snarl pairing.** Paper now describes the implemented boundary-key FIFO rule (default, linear) and
      the new `--pairing interval` anchor-interval order-preserving DP (`src/align.cpp`,
      `correspond_snarls_interval`). Algorithm 1 and the "linear" claim match the code. Boundary-key drift
      (R2.2) is explained and measured (Q1, Fig. 4).
- [x] **Which code produced the ≥80%?** It was the gtest assertion threshold, not a measurement; the paper
      now reports the measured values (P=1.000, R=0.999–1.000 up to 10^5 snarls) and says so. `exp1/run.sh`'s
      inline sweep (which parsed the wrong VCF column) was replaced by
      `experiments/exp1_correctness/pairing_error_sweep.py` (5–10^5 snarls, both deltas).
- [x] **Anchors.** Used for the α term and, in interval mode, to pin the coordinate systems; text and
      Fig. 1 caption corrected; degree-rule statistics on real vg construct graphs added
      (`experiments/plot/graph_degree_stats.py`).
- [x] **Unpaired snarls** now contribute all alleles to the delta (`build_deltas`, shared by sequential and
      parallel aligners; `AlignmentResult::unpaired_snarls_a/b`; tests in `tests/test_pairing.cpp`).
- [x] **Shared graph** → "shared alignment (GAF-style file)" consistently in abstract/Fig. 1/Algorithm 1.
- [x] **Anchor degree rule on real graphs** (HLA vg construct): 50% of nodes pass, >80% of them 1-bp SNP
      alleles, ~15% ≥15 bp; 10–13% of bubbles share a boundary key. vg construct default (32-bp chopping)
      confirmed.
- [x] **`vg snarls` backend** works (`--snarls vg`, vg 1.76) and is byte-identical to the built-in finder
      on the HLA graphs after the decomposer fix (below). Text no longer claims the built-in finder "fuses"
      decomposition with allele enumeration; both backends share it.
- [x] **HPRC/MC runs:** none were ever run; the paper keeps "no results on nested Minigraph-Cactus graphs".
- [x] **NEW — decomposer bug fixed.** The old first-common-successor sink rule dropped the second of two
      adjacent variants on vg construct graphs (HLA A-only recall 0.74). Replaced by a superbubble scan
      (Onodera et al. 2013) in `src/decompose.cpp`: recall 0.99, identical to `vg snarls`, ~20× faster.
      All 161 tests pass; synthetic exactness and linear scaling unchanged. Q4 was re-run (below).

## 3. Experiments

- [x] Pairing error vs. graph size (5–10^5 snarls, 1/2/8-bp boundaries, 0/10% private snarls, 3 seeds,
      both pairing rules): `results/exp1b_pairing_error.csv`, Fig. 4 (`fig_pairing_error`).
- [x] Interval-restricted, order-preserving pairing implemented and evaluated (F1 ≥ 0.999 everywhere).
- [ ] Semi-linear baseline (Giraffe): not run; paper still says "not benchmarked".
- [ ] Hybrid LPT + chunked work-stealing scheduler at 128 threads: not run (no 128-core node reachable;
      IBEX needs key-based SSH, see below); still future work in the paper.
- [x] Real-data check of both modes on HLA-DRB1 vg construct graphs (new Q1 paragraph;
      `experiments/plot/hla_truth_check.py`).
- [x] Q4 chr22 re-run with the fixed decomposer (local WSL, vg 1.76/bcftools 1.24): delta 608,718 alleles
      (427,685 AFR-only, 181,033 EAS-only; was 587,010). Enrichment numbers updated in Q4, Fig. 8 and
      `demo/sample_outputs/*.tsv`.

## 4. Numbers and claims

- [x] 78% DP-phase efficiency = 50.1/64 (DP phase) — kept.
- [x] "~67% k=2" was the prototype workload; the scaling workload has 50.4% k=2 tasks
      (`results/scaling_workload_k_hist.csv`) — Q2 corrected.
- [x] Semi-linear cost "should also grow linearly in h" kept as an expectation.
- [x] Related-work claims checked (POA aligns two POGs; Kececioglu & Starrett 2004 NP-complete; Jain et al.
      2020 NP-complete with graph edits; ODGI path similarity within one graph; HPRC used vg deconstruct).

## 5. Bibliography

- [x] Six new entries checked (DOIs/pages). `xu2025parsga` author list completed
      (Banerjee, Gibney, Xu, Aluru). `onodera2013superbubbles` added. URL/abstract fields stripped for space.
- [x] Dropped for space (uncited now): santa_paper, segram_damla_onur, Secomandi2025, roberts2004minimizers,
      dagum1998openmp. Uncited legacy entries left in the .bib (BibTeX ignores them).

## 6. Build and layout

- [x] Compiles (pdflatex + bibtex) without undefined references; **8 pages** including references.
- [x] `definition` environment renders with IEEEtran; cross-references resolve.
- [x] Algorithm 1 fits its column.
- [x] Author block uncommented. ⚠ **USER**: add the IEEE copyright notice / PDF eXpress step once BIBM sends
      the camera-ready instructions (the CFP says the format will be specified at acceptance).
- [x] Figures merged for space: Gini + matching cost (Fig. 3), parallel scaling + throughput (Fig. 5);
      vEB schematic, F_ST-sweep figure and the scaling table removed (content kept in text).

## 7. Text polish

- [x] British spellings, stray parentheses, "this are", "Limiations", "first system to address it" fixed;
      consistency pass with Definition 2 done.

## 8. Response to reviewers

- [x] `latex/response_to_reviewers.md` (point-by-point; reviewer comments are paraphrased — paste the
      originals before sending).

## Cluster re-measurements (ACCL Slurm cluster, EPYC 7763 nodes; `experiments/accl/`, `results/accl/`)

IBEX accepts only password SSH for the local key and the ROCS head node `mcmgt01` needs the lab VPN
profile, so the re-measurements ran on the ACCL cluster (`acclhead1`, key auth works), which is the node
class of the submitted paper. Its per-job limit is 48 CPUs / 192 GB, so the 64/128-thread and 8-Gbp points
could not be repeated; everything else was:

- [x] Build + 161 tests with GCC 9.4 (needed a span-shim fix for temporaries).
- [x] Pairing-error sweep: identical to the laptop run; both rules 1.8–2.0 s at 10^5 snarls (Q1 updated).
- [x] Parallel scaling 1–48 threads with the new decomposer: serial set-up 7.8 → 3.5 s, single-thread
      20.6 → 15.1 s, end-to-end speedup at 32 threads 2.59 → 3.59× (Q2 sentence added; the figure keeps the
      original dedicated-node 1–64-core curve because the 48-CPU allocation is hardware threads on a
      shared node: DP-phase 21.9× at 32, 26.3× at 48).
- [x] Ablation: 71×, 1773×, 1.82×, 1.15×, 1.00× (Q3 and Fig. 7 updated; submitted values kept in the text).
- [x] Throughput vs linearize-and-map, h = 2…200: PanG2G 7.6 → 12.3 s, baseline 7.0 → 624 s, **51× at
      h = 200** (was 12.7×; the old decomposer dominated PanG2G's time). Q2, Fig. 5b, Conclusion updated;
      the submitted table is kept as `demo/sample_outputs/throughput_chr22_submitted.tsv`.
- [x] Q4 chr22 on the cluster reproduces the laptop numbers (608,718 deltas; 1.46× / 1.61×; flat F_ST sweep).
- [ ] Single-core scaling sweep to 4 Gbp with the new decomposer: the first run shared its node with the
      throughput job (4 Gbp in 1481 s vs 1880 s; decomposition 54 s vs 745 s; memory unchanged); a clean
      single-tenant re-run is queued (`cr_scaling_limit.slurm -w acclnode04`) and will replace that sentence
      in §IV-C if the numbers move.
- [ ] Not possible under the 48-CPU limit: the 128-thread hybrid-scheduler experiment (still future work).
