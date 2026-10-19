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

to generate the xml files for each segment (L and S). Note that we also generate xml for all seqs and those without Nigeria, as when Nigeria was included, it got tricky to converge.

## Run BEAST

Ran it on gpu_p100_long with 4 gpus.

```bash
run-beast.sh
```

## Analyze BEAST output

An export of the output files was sent to S. Dellicour for analysis in R with the `seraphim` package.

```bash
simon-export.sh
```

He gave us this image back:

![seraphim-dispersal-history](Dispersal_history.pdf)
