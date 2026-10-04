# BIBM 2026 camera-ready TODO

Branch: `claude/bibm-reviewer-comments-7fb1c9` (commits 1d90708, f2d4622, dba29ca).
R1/R2/R3 = reviewer #1/#2/#3.

## 1. Blocking: code and data release (R1.3, R2.1, R2's score depends on it)

- [ ] Make `github.com/ahmadshakleya/G2G` public (it currently returns 404).
- [ ] Check that the repo contains no secrets, cluster paths or credentials before making it public (e.g. KAUST ROCS/SANDS paths in `experiments/`).
- [ ] Add a LICENSE and a README with build and run instructions (the root `README.md` only says `# G2G`).
- [ ] Archive a tagged release on Zenodo and put the DOI where `% TODO` sits in `Sections/5_Conclusion.tex`.
- [ ] Verify the availability claim that the repo "reruns every experiment and figure reported here". In particular:
  - [ ] The paper says the AFR/EAS graphs come from **1000 Genomes phase-3 + vg construct**, but `experiments/exp6_hprc/` appears to build them from **HPRC assemblies**. Add the matching scripts or fix the paper.
  - [ ] Scripts for the reference-anchored diff sweep (8 Mbp to 8.2 Gbp) and the F_ST sweep are included.
  - [ ] `make_figures.py` regenerates every figure in `latex/figs/`.
- [ ] Mention the 1000 Genomes panel URL/version, gnomAD-SV version and ClinVar release date in the repo (or the paper).

## 2. Code/paper mismatches (would be visible to anyone reading the public repo)

- [ ] **Snarl pairing.** The paper (Methodology, "Snarl decomposition and pairing") says greedy "highest-scoring unused B-partner, in topological order, anchor-guided". The C++ engine (`src/align.cpp`, `correspond_snarls`) pairs by **exact source/sink sequence key, FIFO among collisions** (since commit e05a276, 2026-06-27). The Python prototype is greedy but compares every pair of snarls, so it is quadratic, not linear as the paper says.
  - [ ] Decide which version the paper describes, and rewrite the pairing paragraph, Algorithm 1 (`PairSnarls` "greedy, anchor-guided") and the "linear in the number of snarls" claim to match.
  - [ ] Note that boundary-key FIFO pairing can drift: one extra snarl with a colliding key shifts every later pair with that key (R2.2). Short 1-bp boundary nodes make collisions frequent.
- [ ] **Which code produced the ≥80% precision/recall (Q1)?** Before or after e05a276, C++ or Python? Also, `experiments/exp1_correctness/run.sh` sweeps 5–200 snarls, but the paper says 5–50.
- [ ] **Anchors are unused for pairing.** `align_snarls` ignores the anchors argument; they only add α·sim to W. Either wire them in or stop claiming they "pin the coordinate systems" (Introduction, Methodology stage 1, Fig. 1 caption, Algorithm 1).
- [ ] **Unpaired snarls are dropped from the delta.** `build_deltas` only emits unmatched alleles of paired snarls. This contradicts Definition 2 (Problem Definition) and the "every allele present in A but absent from B" claim (Introduction, Output subsection). Fix the code or the text.
- [ ] **Shared graph.** It is not written as a graph, only as anchor/snarl lines in the GAF. The paper now says that in the Algorithm 1 caption; check the abstract, pipeline figure and Fig. 1 caption still read consistently.
- [ ] **Anchor degree rule on real graphs.** Check how many nodes pass `in/out-degree ≤ 1` on the vg construct chr22 graphs. The new text assumes vg construct's 32-bp node chopping; confirm the default was used.
- [ ] **`vg snarls` backend.** Methodology now says the built-in finder "fuses decomposition with allele enumeration". Confirm, and confirm that the `G2G_SNARLS=vg` path works.
- [ ] **HPRC/MC runs.** If PanG2G was ever run on real Minigraph-Cactus (HPRC) graphs (`exp6_hprc/03_run_g2g_on_hprc.sh`), report it (R3 asked directly). Otherwise keep "this paper reports no results on nested Minigraph-Cactus graphs".

## 3. Experiments worth running (strengthen answers to reviewers)

- [ ] Pairing error vs. graph size: run exp1 at 10^3–10^5 snarls (R2.2, R3 "20% is high"). Add a sentence or small plot with the trend.
- [ ] Optional: interval-restricted, order-preserving pairing between anchors. The paper proposes it as the fix; even a prototype result would answer R2.2 well.
- [ ] Optional: semi-linear baseline (Giraffe), since R3 asked. The paper currently says "not benchmarked".
- [ ] Optional: hybrid LPT + chunked work-stealing scheduler at 128 threads (R2.3). The paper currently says it's future work.

## 4. Verify numbers and claims added in the revision

- [ ] "Up to 64 cores (78% DP-phase efficiency)": computed as 50.1/64. Confirm it is the DP-phase number.
- [ ] "~67% of tasks are k=2 bubbles" used in Q2: taken from Finding 2 (synthetic, HPRC-calibrated). Confirm it applies to the scaling workload too.
- [ ] Semi-linear paragraph: "cost should also grow linearly in h" is an expectation, not a measurement. Keep it phrased that way.
- [ ] New Related Work claims:
  - [ ] POA "can align two partial-order graphs" (Lee 2002).
  - [ ] Aligning alignments is NP-complete (Kececioglu & Starrett 2004).
  - [ ] Sequence-to-graph alignment is NP-complete with edits on the graph (Jain 2020).
  - [ ] ODGI computes path-level similarity within one graph.
  - [ ] HPRC release used `vg deconstruct` for VCFs.

## 5. Bibliography

- [ ] Check the six new entries (DOIs, pages, author lists): `lee2002poa`, `kececioglu2004aligning`, `jain2020complexity`, `riesen2009approximate`, `singh2008isorank`, `guarracino2022odgi`.
- [ ] `xu2025parsga` key vs. author (Banerjee et al.). Fine for BibTeX but check the author list.
- [ ] Uncited entries (BibTeX ignores them; delete for tidiness if you like): `SMITH1981195`, `li2009bwa`, `langmead2012fast`, `rautiainen2020graphaligner`, `ebler2022pangenie`, `kim2019hisat2`, `minkin2013sibelia`, `Hyun2022`, `Liao2023`.

## 6. Build and layout

- [ ] Compile (pdflatex + bibtex ×2). The edits were never compiled.
- [ ] Check the page count against the BIBM camera-ready limit (Problem Definition, Related Work and Results all grew). Trim if needed.
- [ ] Check that the new `definition` environment (amsthm, `\theoremstyle{definition}`) renders correctly with IEEEtran.
- [ ] Check that the new cross-references resolve: `sec:problem`, `sec:refdiff`, `eq:score`/`eq:dp` (now defined in a different subsection), `fig:delta` from the anchor paragraph.
- [ ] Algorithm 1 still fits its column after the extra lines and longer caption.
- [ ] Uncomment the author block if needed for camera-ready; check the IEEE copyright notice and PDF eXpress requirements.

## 7. Text polish

- [ ] British spellings left in Results Q2: "linearise", "centralised", "serialises" (one changed to "serializes"), "characterises" (Methodology).
- [ ] Results Q2: stray `)` in "...with minimap2 \cite{li2018minimap2})."
- [ ] Methodology decomposition paragraph: "this are the large SV snarls" → "these are".
- [ ] Discussion: subsection title typo "Limiations".
- [ ] Results Q1: "whether the parallel engine reproduces the sequential result)" has a stray `)`.
- [ ] "the first system to solve it" (Introduction contributions): pairing is not optimal, so consider "to address it".
- [ ] Read the whole paper once for consistency between the new Problem Definition (Definition 2) and the rest.

## 8. Response to reviewers

- [ ] Write the point-by-point response letter (R1: 3 points, R2: 3, R3: 6 major + 6 minor).
