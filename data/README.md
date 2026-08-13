# Local data directory

Raw and restricted datasets are not versioned. Obtain them from the accessions
and sources in `../manifests/datasets.tsv`, retain their original terms, and
place local inputs under `data/external/` or a separately declared truth root.

Committed redistributable derived tables live in `../source_data/`. Generated
local derived objects belong under `data/permitted_derived/`; both local data
directories are ignored by Git.
