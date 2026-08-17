# Data access and redistribution

The repository commits only derived numeric Source Data classified for
redistribution. It does not distribute raw bulk RNA-seq matrices, single-cell
count matrices, restricted clinical objects, or private identifier maps.

Use `manifests/datasets.tsv` for the accession, publication, required local
objects, schema, and role of each cohort. `scripts/acquire_data.py` writes an
acquisition checklist; it deliberately does not automate downloads or infer
permission from a public URL. Place obtained inputs under the ignored
`data/external/` tree or provide the declared truth root to an analysis runner.

Every dataset remains governed by its original source and terms. A public
accession does not by itself authorize redistribution. Missing or unresolved
inputs must remain HOLD.

Committed `source_data/*.csv` files are derived numeric tables for figure
and table reproduction. They are not raw expression matrices.

IMvigor210 is accession-only. The repository keeps the paper-reported
summary AUC, CI, P and n in `source_data/figure5.csv`. It does not
redistribute sample-level IMvigor scores, responses or phenotypes.
Figure 5 panel c is therefore not a LEVEL 2 sample-level redraw.

Do not add raw matrices, BAM/FASTQ, private identifier maps, or
institutional paths to this repository.
