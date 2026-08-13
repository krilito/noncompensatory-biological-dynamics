# Source Data classification

This directory contains redistributable derived numeric tables for quantitative
figure redraw and supplementary table snapshots. It contains no raw RNA-seq or
single-cell matrix and no private identifier map.

- `figure2*.csv` through `figure5*.csv`: active main-figure quantitative tables.
- `supplementary_tables/table1.csv` through `table6.csv`: active supplementary
  tables.
- `supplementary_panels/*.csv`: active derived supplementary panel tables.

Use `python scripts/plot_figures.py` for the main quantitative redraw. Final
author-composed artwork is not included.
