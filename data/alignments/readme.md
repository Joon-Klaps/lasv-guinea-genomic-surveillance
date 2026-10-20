# Alignments

Nucleotide alignments of the L and S segments. File names carry the dataset version `LASV_GENBANK_20250725_curated`, followed by the segment. Sequence names have the form `accession|strain|host|country|date`, and the new Guinean genomes are named `|ID-segment||GIN|date`. The workflow that makes them is described in [../../analysis/added-background-validation/readme.md](../../analysis/added-background-validation/readme.md).

## starting-alignment/

The alignments of the sequences that were already curated before this study (493 L and 1,015 S sequences), after cleaning them (misaligned regions were corrected by hand). They are the background to which the new sequences are added.

- `*.background_{L,S}.cleaned.msa` are the alignments that the [lasvphylo](https://github.com/Joon-Klaps/lasvphylo) pipeline starts from.
- `archive/` holds earlier versions: the uncleaned alignments, the version 10 alignments (`LASV_{L,S}_alignment_10_muscle_msa.afa`) and the concatenated-gene versions (`*_concat.maned*.fa`).

## added-background-seqs/

- `*.added_background_{L,S}.fasta` hold the 28 new Guinean genomes per segment, as written by `excelToFasta.py` from the reviewed metadata table.
- `*.added_background_filter.yml` lists, per gene (`gpc`, `np`, `l` and `z`), the sequences that are removed from the starting alignment.
- `lasvphylo-from-scratch/` is the output of a separate lasvphylo run that aligned all sequences from scratch (560 L and 1,062 S sequences) instead of adding the new sequences to the starting alignment. It holds the gene sequences (`subseq/`), the gene alignments (`mafft_orient/`), the concatenated genes (`seqkit_concat/`) and the IQ-TREE run (`iqtree/`).

## final/

The alignments used for the BEAST analyses.

- `*.added_background_{L,S}.renamed.msa` are the final alignments (521 L and 1,043 S sequences).
- `lasvphylo/` is the output of the lasvphylo run that produced them: the gene sequences of the new genomes (`subseq/`), their alignment against reference sequences (`mafft_orient/`), the concatenated genes (`seqkit_concat/`), the starting alignment with the new sequences added (`combined/`, with `.renamed.msa` using the harmonised sequence names) and the IQ-TREE run that was used to spot outliers (`iqtree/`).

## matchingL-S/

The strains that have both a complete L and a complete S segment (153 strains), made by `../../analysis/reassortment/extractCompleteStrains.py`. Sequences are named by strain. The L and S files list the strains in the same order, and the concatenated file lists them in another order.

- `*.renamed.L.msa` and `*.renamed.S.msa` are the two segments.
- `*.renamed.CONCAT.LS.msa` is the concatenation of the two, which is the input of RDP.
