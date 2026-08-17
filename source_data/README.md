# Source Data classification

This directory contains redistributable derived numeric tables for quantitative
figure redraw and supplementary table snapshots. It contains no raw RNA-seq or
single-cell matrix and no private identifier map.

- `figure2*.csv` through `figure5*.csv`: active main-figure quantitative tables.
  IMvigor210 is present only as the committed summary row in `figure5.csv`.
  There is no sample-level IMvigor table.
- `supplementary_tables/table1.csv` through `table6.csv`: active supplementary
  tables.
- `supplementary_panels/*.csv`: active derived supplementary panel tables.

Use `python scripts/plot_figures.py` for the main quantitative redraw and
`python scripts/reproduce_core_results.py` to compare committed numbers with
the preprint contract. Final author-composed artwork is not included.

Redrawing figures from these tables is LEVEL 2 reproduction. It is not a
LEVEL 3 rebuild from raw transcriptomic accessions.
