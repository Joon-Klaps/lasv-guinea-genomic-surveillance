# Lassa virus in Guinea, 2020-2024: phylogeography and reassortment

Code and sequence-level data for the analyses in _Genomic surveillance of Lassa virus in Guinea through in-country sequencing_ (Camara, Annibaldis, Klaps, Ifono _et al._, Emerging Infectious Diseases 2026, [10.3201/eid3205.260386](https://doi.org/10.3201/eid3205.260386)), which is also Chapter 4 of the PhD thesis of Joon Klaps (KU Leuven, 2026).

The study adds 28 Lassa virus (LASV) genomes from human cases in Guinea, sequenced in-country by metagenomic nanopore sequencing, to all public West African LASV sequences. It estimates time-calibrated phylogenies and continuous phylogeography for the L and S segments, and tests the concatenated segments for reassortment. The new genomes are in GenBank under accessions PV847661 to PV847668 and PX115263 to PX115312.

## What this repository contains, and what it does not

The repository holds sequences, alignments, trees, curated GenBank metadata, BEAST X inputs and logs, and the scripts and notebooks that connect them. It holds no patient-level clinical data. The metadata of the new genomes are limited to sequence identifier, sampling date, location with a precision class, and coverage.

Posterior tree files (`*.trees`, several GB per run) are not in the repository. The maximum clade credibility trees summarised from them are (`analysis/phylogeographic/runs/{L,S}/*.combined.ann.tree`). Large text files (`*.log`, `*.tsv`, `*.fa`, `*.fasta`) are stored with [Git LFS](https://git-lfs.com), so install it before cloning.

## Layout

| Path                                                                                                                      | Content                                                                                                                                               |
| ------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `data/genbank/`                                                                                                           | NCBI Virus download of _Mammarenavirus lassaense_ (taxid 3052310) of 25 July 2025                                                                     |
| `data/curated-metadata/`                                                                                                  | Curated metadata per sequence, with the reason each excluded sequence was dropped, and the files used for manual checks                               |
| `data/new-guinea-data/`                                                                                                   | The 28 new Guinean genomes (identifier, sampling date, location, coverage, GenBank IDs) and the pairwise statistics of the Conakry transmission chain |
| `data/alignments/`                                                                                                        | Alignments of the L and S segments, including the matched and concatenated L and S alignment                                                          |
| `data/trees/`                                                                                                             | Maximum-likelihood trees used to check the alignments                                                                                                 |
| `analysis/data-preprocessing/`                                                                                            | Curation of the GenBank metadata and deduplication (`metadata-curation.ipynb`, with a Rust parser for GenBank XML)                                    |
| `analysis/added-background-validation/`                                                                                   | Alignment and validation of the public background sequences with the [lasvphylo](https://github.com/Joon-Klaps/lasvphylo) pipeline                    |
| `analysis/phylogeographic/`                                                                                               | BEAST X XML templates, the scripts that build and run the analyses, and the final runs in `runs/`                                                     |
| `analysis/reassortment/`                                                                                                  | Recombination and reassortment tests (RDP, IQ-TREE approximately unbiased test)                                                                       |
| `analysis/tree-visualisation/`                                                                                            | Figures of the time-calibrated trees and dispersal history, made with [ggphylogeo](https://github.com/Joon-Klaps/ggphylogeo)                          |
| `data/archive/`, `data/trees/archive/`, `data/alignments/starting-alignment/archive/`, `analysis/archive/run-Simon/*.xml` | Earlier curated metadata, alignments and trees that the curation notebook and `excelToFasta.py` read as inputs                                        |

## Workflow

1. **Collect and curate the public sequences.** `data/genbank/LASV_GENBANK_20250725.csv` is the NCBI Virus download. `analysis/data-preprocessing/metadata-curation.ipynb` removes duplicates with MMseqs2 (the identity, coverage and mismatch thresholds are set in the notebook), drops sequences shorter than 500 bp, sequences from patents or vaccines and sequences without sufficient geographic precision, keeps Guinea, Sierra Leone and Liberia for the main analysis, and geocodes the rest (GeoNames, with manual corrections). The result is `data/curated-metadata/LASV_GENBANK_20250725_curated.*.xlsx`. See [analysis/data-preprocessing/readme.md](analysis/data-preprocessing/readme.md).
2. **Align.** `analysis/added-background-validation/` extracts the sequences, aligns the genes against reference sequences with MAFFT, concatenates them in a consistent orientation for codon partitioning, and builds an IQ-TREE tree to spot outliers. Misaligned regions were corrected by hand. The final alignments are in `data/alignments/final/` and `data/alignments/matchingL-S/`.
3. **Build the BEAST X inputs.** `analysis/phylogeographic/make-xmls.sh` calls `preparexml.py` to fill `xml-templates/template_codon_partitioned.xml` (GTR+Gamma with codon partitions, uncorrelated lognormal relaxed clock, Bayesian SkyGrid, relaxed random walk for continuous space) for each segment.
4. **Run BEAST X.** `run-beast.sh` runs one XML on a GPU node. Each segment ran twice from different seeds (112345 and 212345) for 500 000 000 states, logging every 100 000 states. The final inputs, logs and operator summaries are in `analysis/phylogeographic/runs/{L,S}/run_{1,2}/`, and the XML in `runs/{L,S}/` is the file the two runs share. The runs used BEAST `v10.5.0-beta6-dev` (prerelease #9837a1190f) with BEAGLE on GPU, as recorded in the header of each log.
5. **Summarise.** `runs/logCombiner.sh` combines the two runs after a burn-in of 10 000 000 states and `runs/treeAnnotater.sh` summarises the combined posterior as a maximum clade credibility tree (`runs/L/L.12.combined.ann.tree`, `runs/S/S.12.combined.ann.tree`).
6. **Visualise.** `analysis/tree-visualisation/` draws the trees and the dispersal history with ggphylogeo and ggtree. The dispersal-history analysis with the `seraphim` R package was run by S. Dellicour on the exported run files (`analysis/phylogeographic/simon-export.sh`).
7. **Test for reassortment.** `analysis/reassortment/` extracts the strains with complete L and S segments, scans the concatenated alignment with RDP, and compares the placement of the two segments with the approximately unbiased test in IQ-TREE. See [analysis/reassortment/readme.md](analysis/reassortment/readme.md).

## Software

MMseqs2 v14.7e284, MAFFT v7.508, IQ-TREE v2.1.4, BEAST X (see step 4 for the exact build), TreeAnnotator v10.5.0, RDP, ggphylogeo v0.1.2 ([doi:10.5281/zenodo.18507142](https://doi.org/10.5281/zenodo.18507142)) and ggtree v4.0.4. The consensus genomes were generated with [ViMOP](https://github.com/opr-group-bnitm/vimop). The Nextflow runs used to build the alignments are described in the [lasvphylo](https://github.com/Joon-Klaps/lasvphylo) repository.

## Third-party data

The river layer in `analysis/tree-visualisation/Environmental_files/LASV_rivers` comes from the FAO GeoNetwork and the base map from Natural Earth. The public sequences remain the property of their submitters. Cite the original GenBank records, listed in `data/curated-metadata/`, when you reuse them.

## Acknowledgements

We thank our collaborators at the Bernhard Nocht Institute for Tropical Medicine (BNITM), Germany, and their partners at the Irrua Specialist Teaching Hospital (ISTH), Nigeria, for their valuable support.

## Citation

Camara J, Annibaldis G, Klaps J, Ifono K, Koundouno F, Sidibé Y, et al. Genomic Surveillance of Lassa Virus through In-Country Sequencing, Guinea. Emerg Infect Dis. 2026;32(5):819-824. https://doi.org/10.3201/eid3205.260386
