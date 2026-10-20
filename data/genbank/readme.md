# GenBank data

- `LASV_GENBANK_20250725.csv` is the NCBI Virus download of all _Mammarenavirus lassaense_ nucleotide records (taxid 3052310) on 25 July 2025.
- `LASV_GENBANK_subset4annotations.acc.csv` lists the accessions that `../../analysis/data-preprocessing/metadata-curation.ipynb` selected for a second download with the full record.
- `LASV_GENBANK_subset4annotations.gbc.xml` holds the complete GenBank records (INSDSeq XML) of those accessions. It contains the sequences, which the notebook extracts with the Rust parser in `../../analysis/data-preprocessing/xml-parser/` to compare them with MMseqs2. It is stored with Git LFS.

To rebuild the XML file, upload the accession list to the GenBank batch Entrez page, choose "Complete Record", and download it as INSDSeq XML.
