#!/bin/bash -l
#SBATCH --clusters=genius
#SBATCH --account="lp_phylogeo_inf_gpu"
#SBATCH --mail-type="END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="joon.klaps@kuleuven.be"
#SBATCH -c 9
#SBATCH --mem 45000M
#SBATCH --nodes 1
#SBATCH -p gpu_p100
#SBATCH --gres=gpu:1
#SBATCH --time="24:00:00"

conda activate nf-core

# Running the excel conversion script
python excelToFasta.py --output ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background

# Check if everything went all right
if ! [ -f ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background_L.fasta ] && [ -f ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background_S.fasta ]; then
    echo "Something went wrong with file creation."
    exit 1
fi

# Running the nextflow pipeline Joon-Klaps/lasvphylo

nextflow run /user/leuven/344/vsc34477/staging/jklaps/lasvphylo \
    -c /user/leuven/344/vsc34477/data/nf-core/configs/conf/vsc_kul_uhasselt.config \
    -ansi-log false \
    -profile singularity,wice \
    -resume \
    --modify_list $(realpath ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background_filter.yml) \
    --input_id_S LASV_GENBANK_20250725_curated.added_background_S \
    --input_id_L LASV_GENBANK_20250725_curated.added_background_L \
    --input_S $(realpath ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background_S.fasta) \
    --input_L $(realpath ../../data/alignments/added-background-seqs/LASV_GENBANK_20250725_curated.added_background_L.fasta) \
    --alignment_S $(realpath ../../data/alignments/starting-alignment/LASV_GENBANK_20250725_curated.added_background_S_concat.maned.cleaned.fa) \
    --alignment_L $(realpath ../../data/alignments/starting-alignment/LASV_GENBANK_20250725_curated.added_background_L_concat.maned.cleaned.fa) \
    --tree_S $(realpath ../../data/trees/LASV_GENBANK_20250725_curated.added_background_S_concat.fa.cleaned.treefile) \
    --tree_L $(realpath ../../data/trees/LASV_GENBANK_20250725_curated.added_background_L_concat.fa.cleaned.treefile) \
    --outdir $(realpath ../../data/alignments/added-background-seqs/lasvphylo)

