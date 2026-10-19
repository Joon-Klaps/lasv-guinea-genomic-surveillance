# Reasortment testing

The samples from the human-to-human transmission cluster had the L & S segment phylo placements in completely different clades, which could be a sign of reassortment.

To test this, we used:

- the RDP5 software to look for evidence of recombination in the L and S segments.
- IQTRee with an approximately unbiased (AU) test (Shimodaira, 2002) via -au option:
  `iqtree -s example.phy -z example.treels -n 0 -zb 1000 -zw -au`

However to do this, we need a tree topology of samples that have both complete L & S segments.

## Data preparation

To do this, we used the `extractCompleteStrains.py` script to extract the strains that have both complete L and S segments, and then generated the corresponding fasta files for these strains.

```bash
python extractCompleteStrains.py
```

## IQTREE

We then ran IQTREE with the AU test to compare the phylogenetic placements of the L and S segments for these strains.

```bash
sbatch iqtree.slurm
# iqtree -B 1000 -m GTR+G4 -T AUTO
```

### RDP5

For RDP5 we need to the sequences concatenate first

```bash
seqkit concat \
    ../../data/alignments/matchingL-S/LASV_GENBANK_20250725_curated.added_background.renamed.L.msa \
    ../../data/alignments/matchingL-S/LASV_GENBANK_20250725_curated.added_background.renamed.S.msa \
    > ../../data/alignments/matchingL-S/LASV_GENBANK_20250725_curated.added_background.renamed.CONCAT.LS.msa
```

quick check

```bash
rg -c ">" ../../data/alignments/matchingL-S/LASV_GENBANK_20250725_curated.added_background.renamed.CONCAT.LS.msa
# 153
```

Then we can load the concatenated fasta file into RDP5 and run the recombination analysis.
