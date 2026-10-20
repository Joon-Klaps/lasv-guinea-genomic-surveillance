# Maximum-likelihood trees

IQ-TREE trees that were used to check the alignments and to test for reassortment. The time-calibrated trees of the BEAST analyses are in [../../analysis/phylogeographic/runs/](../../analysis/phylogeographic/runs/).

- `*.added_background_{L,S}_concat.fa.cleaned.treefile` are the trees of the starting alignments after cleaning. The lasvphylo pipeline takes them as starting trees (`--tree_L` and `--tree_S` in `../../analysis/added-background-validation/run-pipeline.sh`).
- `archive/` holds the trees before cleaning and the trees of the earlier version 10 alignments.
- `matchingL-S/` holds the IQ-TREE runs (GTR+G4 with 1,000 ultrafast bootstraps, `../../analysis/reassortment/iqtree.slurm`) of the L and S alignments of the 153 strains that have both segments, which are in `../alignments/matchingL-S/`. Each run has the maximum-likelihood tree (`.treefile`), the bootstrap consensus tree (`.contree`), the report (`.iqtree`), the log (`.log`), the distance matrix (`.mldist`), the starting trees (`.bionj` and `.parstree`), the splits (`.splits.nex`) and the checkpoint (`.ckp.gz`).
