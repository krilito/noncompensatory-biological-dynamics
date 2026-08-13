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
