# Working notes on the curation of the GenBank dataset

These notes predate the current README and the thesis. They are kept for the curation decisions they record and may be out of date. The current description of the workflow is in [readme.md](readme.md).

Goal is to continue the project that describes the LASV sporadic outbreaks in Guinea from 2020 to 2024.

--Might be outdated--

## Data Preparation

**Updated July 25, 2025**: NCBI Virus GenBank was queried for new LASV sequences using the following query:

```url
https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Mammarenavirus%20lassaense,%20taxid:3052310
```

The downloaded data is stored in [`data/genbank/LASV_GENBANK_20250725.csv`](data/genbank/LASV_GENBANK_20250725.csv).

**It seems like some sequences from [Wiley et al., 2019](manuscripts/A%20Sporadic%20and%20Lethal%20Lassa%20Fever%20Case%20in%20Forest%20Guinea%202019.pdf) were removed in this step.**: MG812630-MG812685 and MH215278-MH215289, most of these sequences were being removed from alignment 2, meaning duplicates. However these samples are from wide range of dates, hence uncertain why they weren't kept. We reincluded some of them after reinterpreting the metadata & joining sequences.

Samples from [Bangura et al. 2024](manuscripts/Bangura%20et%20al%202024%20Emer%20Mic%20Inf%20-%20Elisabeth%20-%20Spatio-temporal%20spread%20of%20Lassa%20virus.pdf) were included.

## Data Processing workflow (notebook)

A comprehensive data preprocessing workflow has been implemented in [`analysis/data-preprocessing/metadata-curation.ipynb`](analysis/data-preprocessing/metadata-curation.ipynb) that:

### 1. Data Integration and Merging

- Merges newly downloaded GenBank data (2025) with previously curated sequences from 2021
- Integrates historical alignment data from versions 1-11
- Combines country-specific curated datasets (Guinea, Sierra Leone, Liberia, Côte d'Ivoire)
- Handles duplicate columns and missing data through intelligent coalescing

### 2. Advanced Sequence Deduplication

- Uses MMseqs2 for sensitive %ANI-based duplicate detection (≥65% identity, ≥40% coverage)
- Combines sequence similarity with metadata matching (isolate names, references, submitters)
- Prioritizes sequences from established alignments when resolving duplicates
- Identifies and removes near-identical sequences with detailed reasoning

### 3. Quality Filtering and Curation

- **Length filtering**: Removes sequences ≤500bp
- **Content filtering**: Excludes sequences with "patent", "modified", or "vaccine" in titles
- **Geographic filtering**: Limits to West African countries (GIN, LBR, SLE, NGA, MLI, CIV, GHA, BFA, TGO, BEN)
- **Precision filtering**: Removes sequences with only country-level geocoding accuracy
- **Historical preservation**: Maintains exclusion decisions from previous alignment iterations

### 4. Comprehensive Geocoding System

- Multi-tiered location annotation:
  1. Manual corrections for known problematic entries
  2. Coordinate borrowing from related sequences (same submitter/location)
  3. GeoNames API geocoding for missing coordinates
- Assigns precision levels: coordinates_GenBank, manually_verified, city_geonames, region_geonames, country_geonames
- Generates detailed correction comments and validation

### 5. Quality Control and Validation

- Geographic validation through interactive maps
- Visualization of sequence distributions across West Africa
- Coordinate accuracy validation and outlier detection

## Key Outputs

**Primary Dataset**: [`data/curated-metadata/LASV_GENBANK_20250725_curated.raw.xlsx`](data/curated-metadata/LASV_GENBANK_20250725_curated.raw.xlsx)

- Complete curated dataset with all sequences, curation flags, and detailed metadata
- Includes `kept_in_alignment` flags indicating sequences suitable for phylogenetic analysis
- Contains `discarded_reason` column with detailed exclusion justifications
- Preserves all historical curation decisions and alignment version tracking

**Intermediate Files**:

- [`data/temp/sequences_for_mmseqs.fasta`](data/temp/sequences_for_mmseqs.fasta) - Sequences used for deduplication
- [`data/temp/mmseqs_results.tsv`](data/temp/mmseqs_results.tsv) - Similarity search results

## Technical Implementation

### Dependencies

- **R packages**: tidyverse, readxl, openxlsx, countrycode, xml2, geonames, RColorBrewer, data.table
- **External tools**: MMseqs2 (conda environment), custom Rust XML parser
- **APIs**: GeoNames for geocoding (requires username configuration)

### Key Innovations

- **Rust-based XML parser**: Fast extraction of sequences from GenBank XML format
- **Intelligent duplicate resolution**: Preserves sequences from established alignments
- **Coordinate borrowing**: Leverages related sequences for missing geographic data
- **Comprehensive traceability**: Maintains full history of curation decisions

## Historical Context

### Previous Work (2020-2021)

The original curation work by Liana and collaborators established important filtering criteria:

- **Sim & Raph 2020**: Focused on human sequences only with initial deduplication
- **Liana 2020**: Comprehensive curation but noted gaps:
  - Nigerian sequences needed rechecking
  - Benin sequences required validation
  - Guinea sequences were removed due to insufficient similarity
  - Location precision was a major limitation

### Current Improvements (2025)

- **Expanded geographic scope**: Systematic inclusion of all West African countries
- **Enhanced deduplication**: MMseqs2-based approach with metadata validation
- **Improved geocoding**: Multi-tiered system with precision tracking
- **Better integration**: Seamless merging of historical and new data
- **Quality assurance**: Geographic validation and comprehensive documentation

## Alignment checking:

- Within the L alignment, there tend to be some missaligned regions through large gaps, need to add the original reference genes to them in order to check if the codon positions are still correct.

- Within the S alignment, some regions at the end need to be checked and removed, they tend to be misaligned as well, when looking it up is the NP region in the GPC2 region so that genome has NP - GPC1.

Also excel To fasta isnt removing form L segment original sequences from Guinea.

## Phylogeographic analysis with BEAST

check the [xml template](analysis/phylogeographic/xml-templates/template_codon_partitioned.xml) to see which variables we've been using.

()
