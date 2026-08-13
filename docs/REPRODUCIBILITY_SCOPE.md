# Reproducibility scope

The public package supports two distinct workflows:

1. **Quantitative redraw:** regenerate Figures 2–5 from committed derived CSV
   tables with `python scripts/plot_figures.py`.
2. **Analysis recomputation:** recompute registered chains from separately
   acquired source datasets with `python scripts/recompute_results.py
   --truth-root <path>`.

The redraw path does not claim to reproduce final author-composed Adobe boards.
The analysis path fails closed when required data, schema, or provenance are
missing. Figure 1 is conceptual and has no numeric producer.

The release is generated one way from a controlled upstream workbench. Its Git
history is independent, and it is not an input authority for frozen values.
