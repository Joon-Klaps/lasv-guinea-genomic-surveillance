# Data

Sequence-level data for the analyses in `../analysis/`. The files in `data/` hold no patient-level clinical data. Large text files (`*.fasta`, `*.fa`, `*.log`, `*.tsv`) and a few large inputs are stored with [Git LFS](https://git-lfs.com), so install it before cloning.

| Directory | Content |
| --- | --- |
| [genbank/](genbank/readme.md) | The NCBI Virus download of 25 July 2025, and the GenBank records of the sequences that were deduplicated |
| [curated-metadata/](curated-metadata/readme.md) | One row per sequence record, with curation flags, the reason a sequence was excluded and coordinate-annotated locations. **Start here if you want the dataset** |
| [new-guinea-data/](new-guinea-data/readme.md) | The 28 new Guinean genomes and the pairwise comparison of the Conakry transmission chain |
| [alignments/](alignments/readme.md) | Starting, intermediate and final alignments of the L and S segments, and the alignments of the strains that have both segments |
| [trees/](trees/readme.md) | Maximum-likelihood trees that were used to check the alignments and to test for reassortment |
| [temp/](temp/readme.md) | Intermediate files of the curation notebook that are cached so that it can be re-run |
| [archive/](archive/readme.md) | Earlier curated metadata, GenBank downloads and alignment lists that the curation notebook reads as input |

## How the files connect

1. `genbank/LASV_GENBANK_20250725.csv` and the files in `archive/` are curated by `../analysis/data-preprocessing/metadata-curation.ipynb` into `curated-metadata/LASV_GENBANK_20250725_curated.raw.xlsx`.
2. That table is reviewed by hand into `curated-metadata/LASV_GENBANK_20250725_curated.manual.xlsx`.
3. `../analysis/added-background-validation/` takes the sequences from the reviewed table, aligns them and writes the alignments in `alignments/`, with the trees in `trees/` used as starting trees.
4. `alignments/final/` is the input of the BEAST analyses in `../analysis/phylogeographic/`.
5. `alignments/matchingL-S/` and `trees/matchingL-S/` are the input of the reassortment tests in `../analysis/reassortment/`.
