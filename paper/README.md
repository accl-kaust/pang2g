# G2G — IEEE BIBM 2026 manuscript

`g2g_bibm.tex` is the manuscript (IEEEtran `conference` class).

## Building

Needs the standard IEEE LaTeX toolchain plus two local inputs that are **not**
committed here:

- `refs.bib` — the bibliography (keep your existing one alongside this file).
- `figs/*.pdf` — regenerate from the demo result TSVs:

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
