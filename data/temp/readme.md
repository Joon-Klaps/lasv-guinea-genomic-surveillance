# Cached intermediate files

These files are written by `../../analysis/data-preprocessing/metadata-curation.ipynb` and are kept so that the notebook can be re-run without MMseqs2 and without a GeoNames account. The notebook only recomputes a file when it is missing.

- `sequences_for_mmseqs.fasta` holds the sequences longer than 500 nucleotides that were compared with each other.
- `mmseqs_results.tsv` holds the MMseqs2 `easy-search` hits between those sequences (minimum identity 80%, coverage 60%).
- `geocode_results.csv` caches the GeoNames lookups of locations.
