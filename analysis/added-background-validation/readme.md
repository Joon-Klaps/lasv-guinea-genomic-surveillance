# Background Sequencing Validation

This document outlines the background sequencing validation process for the LASV GenBank dataset, focusing on the curation and validation of sequences from West Africa. The goal is to ensure high-quality, geographically accurate, and scientifically valid data for further analysis.

## Overview

1. Extract sequences from excel file, [`LASV_GENBANK_20250725_curated.manual.xlsx`](../../data/curated-metadata/LASV_GENBANK_20250725_curated.manual.xlsx), using a mini python script.

```bash
python excelToFasta.py
```

2. run pipeline [lasvphylo](https://github.com/Joon-Klaps/lasvphylo):
   - extract genes from sequences, aligning them with mafft to a set of reference sequences
   - concantenated the genes to each other in the same orientation (usefull for specific codon models within BEAST)
   - add them to the alignment with muscle
   - make a ML tree with IQTree ( for spotting outliers)

3. Validate the alignment & tree.
4. Update the excel file ([`LASV_GENBANK_20250725_curated.manual.xlsx`](../../data/curated-metadata/LASV_GENBANK_20250725_curated.manual.xlsx)) and patch those that have been deemed invalid.
