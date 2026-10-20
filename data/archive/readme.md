# Archive of earlier data

Earlier curated metadata, GenBank downloads and alignment lists. They are not new data of this study, but `../../analysis/data-preprocessing/metadata-curation.ipynb` reads some of them to carry the earlier curation decisions over.

- `alignments/` lists the sequence identifiers of the earlier alignments 1 to 3 per segment (`LASV_{L,S}_alignment_{1,2,3}.csv`). `infoAln.md` explains how the three versions differ. The notebook reads these lists.
- `genbank/` holds NCBI Virus downloads from 2021 (`All-*-sequences-downloaded-NCBI-Virus-*`) and the curated workbooks per country (`LASV_ALL_NCBI_VIRUS_{COTDIVOIRE,GUINEA,LIBERIA,SIERRALEONE}_V*.xlsx`). The notebook reads the four curated workbooks, and the raw downloads are kept for reference.
- `metadata/` holds the curated metadata of the earlier work. The notebook reads `LASV_all_the_metadata-allinfo_LK_v01.xlsx`, which is stored with Git LFS. `LASV_all_the_metadata.xlsx` is a version with another layout, and `LASV_to_be_added_Liana_metadata.xlsx` is a working table with a `to_be_added` and a `reason_for_exclusion` column per sequence. The notebook does not read these two.

The archived trees and alignments are in `../trees/archive/` and `../alignments/starting-alignment/archive/`.
