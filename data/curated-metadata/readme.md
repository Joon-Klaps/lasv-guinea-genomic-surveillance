# Curated metadata

One row per sequence record (one segment of one isolate). The tables hold all sequences from GenBank (25 July 2025) and the earlier curated datasets, with the reason that a sequence was left out, so nothing was deleted. Filter on the columns below to get the subset you need.

| File | Rows | Content |
| --- | --- | --- |
| `LASV_GENBANK_20250725_curated.raw.xlsx` | 3,278 | The output of `../../analysis/data-preprocessing/metadata-curation.ipynb`, before manual review and before any alignment. Contains the nucleotide sequence of each record (unaligned) |
| `LASV_GENBANK_20250725_curated.manual.xlsx` | 3,278 | The same table after manual review. It is the input of `../../analysis/added-background-validation/excelToFasta.py` |
| `LASV_GENBANK_20250725_curated.newseqsremoved.xlsx` | 3,222 | The reviewed table without the 56 records of the 28 new Guinean genomes |

## Columns

- `Accession`, `Strain`, `Segment`, `genes`, `Length`: identification of the record. `genes` lists the genes found in the sequence (`GPC`, `NP`, `Z`, `L`).
- `kept_in_alignment`: `TRUE` if the sequence is used in the analyses and `FALSE` if it is excluded. It is empty when there was no earlier decision about the sequence, and `excelToFasta.py` treats an empty value as kept.
- `discarded_reason`: why a sequence is excluded. The main reasons are a location that is only known at country level, a length of 500 nucleotides or fewer, duplicates, laboratory strains, a missing collection date, a sample from outside West Africa, patent, modified or vaccine sequences, a composition failure in IQ-TREE, and low quality.
- `duplicate_of`: the record that a duplicate was matched with.
- `checked_in_initial_curation`, `last_alignment`: the history in the earlier alignments (see `../archive/alignments/infoAln.md`). `last_alignment` is the number of the latest earlier alignment that contained the sequence.
- `Country`, `admin1`, `location`, `Geo_Location`: the location as reported in GenBank and standardised.
- `latitude`, `longitude`, `precision`: the coordinates that are used, and how precise they are. `precision` records where they come from, for example `coordinates_GenBank`, `manually_verified`, a city, `admin-1`, `admin-2` or `country`. `latitude-old` and `longitude-old` are the coordinates of the earlier curation.
- `Collection_Date`, `Release_Date`, `Year`, `Host`, `Tissue_Specimen_Source`: sampling information from GenBank.
- `Submitters`, `Organization`, `Publications`, `reference`, `source` (usually a DOI), `BioProject`, `BioSample`, `SRA_Accession`: provenance.
- `sequence`: the nucleotide sequence of the record, unaligned.
- `alignment_name`, `old_alignment_name`: the sequence name in the alignments (`accession|strain|host|country|date`) and the name in earlier alignments.
- `to_remove_from_starting_alignment`: `TRUE` for sequences that are removed from the starting alignment. `excelToFasta.py` writes them to `../alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background_filter.yml`.
- `comments`: free-text remarks of the curation.

## manual/

Files that were used during the manual review.

- `duplicates_for_manual_verification.csv`: pairs of records (`Accession`, `duplicate_of`) that were checked by hand.
- `genes_for_manual_verification.csv`: records for which the genes that the sequence covers were set by hand (`Accession`, `Segment`, `genes`). The notebook reads this file.
- `iqtree_composition_fail.csv`: records that failed the composition test in IQ-TREE and were excluded.
- `locations_for_manual_verification.xlsx`: locations as written in GenBank (`Original_Location`) with coordinates that were looked up by hand and the source of the coordinates in `Notes`. The notebook does not read this file.
