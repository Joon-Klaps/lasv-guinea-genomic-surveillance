# Tree visualisation

[`ggphylogeo.ipynb`](ggphylogeo.ipynb) draws the time-calibrated trees and the dispersal history of the L and S segments from the maximum clade credibility trees, using [ggphylogeo](https://github.com/Joon-Klaps/ggphylogeo) (v0.1.2, [doi:10.5281/zenodo.18507142](https://doi.org/10.5281/zenodo.18507142)) and ggtree. It is an R notebook and has to be run from this directory, because all paths are relative.

## Inputs

- `L.12.combined.ann.tree` and `S.12.combined.ann.tree` are symbolic links to the maximum clade credibility trees in `../phylogeographic/runs/{L,S}/`.
- `../phylogeographic/runs/{L,S}/*.metadata_{L,S}.csv` hold the sampling dates and coordinates that were given to BEAST.
- `../../data/new-guinea-data/GIN-LASV-final-ViMOP-seq-for-Leuven_V02.xlsx` links the sequencing identifiers of the new genomes to their GenBank accessions.
- `new-seqs.tsv` flags the sequences that are new in this study.
- `Environmental_files/` holds the base map (Natural Earth) and the river layer (FAO GeoNetwork). `FAO_GeoNetwork.shp` and `.dbf` are stored with Git LFS.

## Outputs

The notebook writes `figures/phylo_LS.pdf` and `figures/phylogeo_LS_with_S_zooms.pdf`.

## The notebook is not the final version

The figures in `figures/` were made after the last edit of this notebook (the PDFs on 20 February 2026 and the PNG and JPG files on 12 March 2026, against 9 February 2026 for the notebook). Running the notebook reproduces the same trees, topologies and dispersal paths, but not the exact appearance of the published figures. Tip labels are GenBank accessions instead of the sequencing identifiers (for example `PX115271` instead of `G0780`), the dispersal map uses another colour scale, and the tree figure has no clade shading or host symbols. The notebook is saved here without its outputs.

## R packages

ggphylogeo, ggtree, phangorn, dplyr, stringr, purrr, tidyr, readxl, rnaturalearth, sf, terra, ggplot2, ggrepel, patchwork and ggforce.
