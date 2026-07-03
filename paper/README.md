# G2G — IEEE BIBM 2026 manuscript

`g2g_bibm.tex` is the manuscript (IEEEtran `conference` class).

## Building

`refs.bib` is committed alongside this file, but it was assembled from knowledge
and **must be verified** before submission — confirm volume/page/DOI against the
primary sources, and resolve the entries marked `% VERIFY` (notably
`xu2025parsga`'s author list and `li2024minigraph`). If you maintain your own
bibliography, reconcile against it.

The figures are **not** committed — regenerate `figs/*.pdf` from the demo result
TSVs:

  ```bash
  python3 ../demo/make_figures.py \
      --ablation-csv ../demo/work/ablation/ablation_results.tsv \
      --scaling-csv  ../demo/work/parallel/parallel_scaling.tsv \
      --biology-csv  ../demo/work/biology/biology_results.tsv
  # writes figs/ next to your cwd; move/symlink it next to the .tex
  ```

Then:

```bash
latexmk -pdf g2g_bibm.tex      # or: pdflatex; bibtex; pdflatex x2
```

## Before submission

Every value still being confirmed is wrapped in `\PH{...}` (renders red). Search
for `\PH` and `\TBD` and resolve each:

- parallel speedups and ablation slowdowns — final numbers from the cluster runs;
- `fig_throughput` is still a sketch (the linearise-and-map baseline sweep is the
  one experiment not yet run);
- the real-HPRC Q4 run, the bib keys flagged at the top of the `.tex`, and the
  author/acknowledgment placeholders.

Remove the `\PH`/`\TBD` macros (or redefine them to no-ops) for the camera-ready.
