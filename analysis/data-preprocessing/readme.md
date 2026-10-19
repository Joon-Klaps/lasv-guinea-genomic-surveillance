# Data Preprocessing

This directory contains the comprehensive data preprocessing pipeline for LASV sequences. The workflow merges newly downloaded GenBank data (2025) with previously curated sequences from 2021, applies rigorous filtering and deduplication, and performs geocoding for location annotation.

## Overview

The preprocessing pipeline is implemented in [`filtering.ipynb`](filtering.ipynb) and follows these major steps:

1. **Data Loading & Integration** - Merge new GenBank data with historical curated datasets
2. **Sequence Deduplication** - Use MMseqs2 for %ANI-based duplicate detection
3. **Quality Filtering** - Apply length, content, and geographic filters
4. **Geocoding & Location Annotation** - Standardize and annotate sequence locations
5. **Quality Control** - Generate maps and validation outputs
6. **Final Export** - Create curated datasets for downstream analysis

## Input Files

### Primary Data Sources

- `../../data/genbank/LASV_GENBANK_20250725.csv` - Recently downloaded GenBank data (July 2025)
- `../../data/archive/metadata/LASV_all_the_metadata-allinfo_LK_v01.xlsx` - Previously curated metadata from 2021

### Historical Alignment Data

- `../../data/archive/alignments/LASV_{S,L}_alignment_{1,2,3}.csv` - Previous alignment iterations

### Country-Specific Curated Data

- `../../data/archive/genbank/LASV_ALL_NCBI_VIRUS_GUINEA_V08.xlsx` - Guinea sequences
- `../../data/archive/genbank/LASV_ALL_NCBI_VIRUS_SIERRALEONE_V02.xlsx` - Sierra Leone sequences
- `../../data/archive/genbank/LASV_ALL_NCBI_VIRUS_LIBERIA_V01.xlsx` - Liberia sequences
- `../../data/archive/genbank/LASV_ALL_NCBI_VIRUS_COTDIVOIRE_V01.xlsx` - Côte d'Ivoire sequences

### Sequence Data for Deduplication

- `../../data/genbank/LASV_GENBANK_subset4annotations.gbc.xml` - GenBank XML with full sequence data
- `../../data/genbank/LASV_GENBANK_subset4annotations.acc.csv` - Accession list for GenBank batch download

## Processing Steps

### 1. Data Loading and Merging

- Load GenBank 2025 data and visualize release patterns by country and quarter
- Import historical curated metadata and country-specific datasets
- Extract sequence IDs from existing FASTA alignments (versions 1-11)
- Merge all data sources with proper handling of duplicate columns and data coalescing

### 2. Sequence Deduplication Using MMseqs2

- Extract full sequences from GenBank XML annotations using custom Rust parser
- Create FASTA file with sequences >500bp for MMseqs2 analysis
- Run MMseqs2 easy-search with parameters:
  - Minimum identity: 65%
  - Coverage threshold: 40%
  - Global alignment mode (both strands)
- Identify duplicates based on:
  - High sequence similarity (>99.9% coverage, <2 mismatches)
  - Matching metadata (isolate name, reference, or submitters)
- Prioritize sequences already in alignments when selecting representatives

### 3. Quality Filtering and Curation

Apply filtering rules with detailed exclusion reasons:

- **Length filter**: Remove sequences ≤500bp
- **Content filter**: Exclude sequences with "patent", "modified", or "vaccine" in title
- **Geographic filter**: Limit to West African countries (GIN, LBR, SLE, NGA, MLI, CIV, GHA, BFA, TGO, BEN)
- **Historical exclusions**: Preserve exclusion decisions from previous alignment iterations
- **Precision filter**: Remove sequences with only country-level geocoding precision

### 4. Geocoding and Location Annotation

- Apply manual location corrections for known problematic entries
- Use geonames API for automated geocoding of missing coordinates
- Implement coordinate borrowing from related sequences (same submitter/location)
- Assign precision levels: coordinates_GenBank, manually_verified, city_geonames, region_geonames, country_geonames
- Handle special cases and generate correction comments

### 5. Quality Control and Validation

- Generate maps showing sequence distributions across West Africa
- Separate visualizations for sequences considered for inclusion vs. all available sequences
- Validate coordinate accuracy and identify potential geographic outliers

## Output Files

### Primary Outputs

- `../../data/curated-metadata/LASV_GENBANK_20250725_curated.raw.xlsx` - **Complete curated dataset** with all sequences, flags, and metadata

### Intermediate Files

- `../../data/temp/sequences_for_mmseqs.fasta` - FASTA file for deduplication analysis
- `../../data/temp/mmseqs_results.tsv` - MMseqs2 similarity search results

### Planned Outputs (not yet generated)

- `../../data/curated-metadata/LASV_GENBANK_20250725_curated.final.csv` - Final filtered dataset for analysis
- `../../data/curated-metadata/LASV_GENBANK_20250725_curated.manual.xlsx` - Sequences requiring manual review

## Key Features and Innovations

### Advanced Deduplication Pipeline

- Uses MMseqs2 for sensitive sequence similarity detection
- Combines sequence similarity with metadata matching
- Preserves sequences from established alignments when resolving duplicates

### Comprehensive Geocoding System

- Multi-tiered approach: manual corrections → coordinate borrowing → API geocoding
- Detailed precision tracking and validation
- Handles country code standardization and coordinate system consistency

### Robust Data Integration

- Merges multiple historical datasets with proper column coalescing
- Preserves historical curation decisions and exclusion reasons
- Maintains traceability through alignment version tracking

### Quality Assurance

- Geographic validation through mapping
- Detailed exclusion reason tracking
- Comprehensive metadata validation and standardization

## Dependencies

- **R packages**: tidyverse, readxl, openxlsx, countrycode, xml2, geonames, RColorBrewer, data.table, processx
- **External tools**: MMseqs2 (via conda environment), custom Rust XML parser
- **APIs**: GeoNames for location geocoding

## Notes

- The Rust XML parser (`xml-parser/target/release/genbank_parser`) must be compiled before running the sequence extraction step
- GeoNames API requires username configuration: `options(geonamesUsername="your_username")`
- MMseqs2 should be available in a conda environment named "mmseqs2-env"
