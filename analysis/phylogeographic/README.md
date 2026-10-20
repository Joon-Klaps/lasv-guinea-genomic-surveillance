# Phylogeographic BEAST analysis of LASV samples

## Prepare XML for BEAST

jinja templates have been made in order to easily modify the xml files for beast without having to call beauti every time.

There are two templates:

- **`xml-templates/template_codon_partitioned.xml`: for codon partitioned analysis**
- **`xml-templates/template_jinja.xml`: for standard analysis (no codon partitioning)**

It's the `xml-templates/template_codon_partitioned.xml` that we will be using for our analyses as it has better weights, uses hamiltonianMonteCarloOperator for the paritioned codon models for faster convergence.

Afterwards we run

```bash
make-xmls.sh # uses preparexml.py
```

to generate the xml files for each segment (L and S). The final analyses leave out the sequences from Nigeria, Côte d'Ivoire, Benin, Togo and Mali (`--skip-countries NIG NGA CIV BEN TGO MLI` in `make-xmls.sh`), because the analysis was hard to get to converge when Nigeria was included. The calls for the other variants are commented out in the script.

## Run BEAST

Each segment ran twice with different seeds, on the `gpu_p100_long` partition with 3 GPUs per run.

```bash
sbatch run-beast.sh <path/to/xmlfile.xml> <output directory>
```

## Final runs

The final XML, logs and operator summaries of the two runs per segment are in `runs/{L,S}/run_{1,2}/`, and `runs/{L,S}/` holds the XML that the two runs share and the metadata (sampling dates and coordinates) that it was built from. The posterior tree files (`*.trees`) are several GB per run and are not in the repository. The summaries are:

- `runs/logCombiner.sh` combines the two runs after a burn-in of 10,000,000 states.
- `runs/treeAnnotater.sh` summarises the combined posterior as a maximum clade credibility tree (`runs/L/L.12.combined.ann.tree` and `runs/S/S.12.combined.ann.tree`), which `../tree-visualisation/` draws.

## Analyze BEAST output

An export of the output files was sent to S. Dellicour for analysis in R with the `seraphim` package.

```bash
simon-export.sh
```

`seraphim-scripts/` holds the R scripts of an earlier seraphim analysis, which read the output of an earlier run (`RRW_phylogeography/LASV_*_alignment_4`). The dispersal history of the final runs is not included in this repository. The dispersal paths in the figures of `../tree-visualisation/` come from the maximum clade credibility trees.
