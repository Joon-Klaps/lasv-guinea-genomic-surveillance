#! /bin/bash -l

conda activate nf-core

# Running the excel conversion script
python excelToFasta.py --output ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background

# Check if everything went all right
if ! [ -f ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background_L.fasta ] && [ -f ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background_S.fasta ]; then
    echo "Something went wrong with file creation."
    exit 1
fi

mkdir -p ../../data/alignments/added-background-seqs/lasvphylo

# Running the nextflow pipeline Joon-Klaps/lasvphylo
nextflow run /Users/joonklaps/Desktop/School/PhD/side_projects/archive/guinea/LVE-PRJ02-supplementary_V02/Suppl/Files_used_in_pipeline_JK/nf-core-lasvphylo \
    -profile docker,arm \
    --max_cpus   10 \
    --max_memory '30.GB' \
    --max_time   '6.h' \
    -bg \
    -resume \
    --modify_list $(realpath ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background_filter.yml) \
    --input_id_S LASV_GENBANK_20250725_curated.added_background_S \
    --input_id_L LASV_GENBANK_20250725_curated.added_background_L \
    --input_S $(realpath ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background_S.fasta) \
    --input_L $(realpath ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background_L.fasta) \
    --alignment_S $(realpath ../../data/alignments/starting-alignment/LASV_GENBANK_20250725_curated.background_S.cleaned.msa) \
    --alignment_L $(realpath ../../data/alignments/starting-alignment/LASV_GENBANK_20250725_curated.background_L.cleaned.msa) \
    --tree_S $(realpath ../../data/trees/LASV_GENBANK_20250725_curated.added_background_S_concat.fa.cleaned.treefile) \
    --tree_L $(realpath ../../data/trees/LASV_GENBANK_20250725_curated.added_background_L_concat.fa.cleaned.treefile) \
    --outdir $(realpath ../../data/alignments/added-background-seqs/lasvphylo)
