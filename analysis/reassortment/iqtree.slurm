#!/bin/bash
#SBATCH --clusters=wice
#SBATCH --account="lp_phylogeo_inf_gpu"
#SBATCH --mail-type="END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="joon.klaps@kuleuven.be"
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time="12:00:00"
#SBATCH --array=1-2%2

####################################
##### Run a iqtree2 analysis ######
####################################

#### How to run: ####
# sbatch run-iqtree2.sh

# write a function that catches any errors
abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

set -e
trap 'abort' 0

locations=(
    "../../data/alignments/matchingL-S/LASV_GENBANK_20250725_curated.added_background.renamed.L.msa"
    "../../data/alignments/matchingL-S/LASV_GENBANK_20250725_curated.added_background.renamed.S.msa"
)

current_location="$(realpath ${locations[$SLURM_ARRAY_TASK_ID-1]})"
echo "Processing file: $current_location"


# write some cluster metadata to the error and output files
hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2



mkdir -p ../../data/trees/matchingL-S
DATA_DIR="$(realpath ../../data)"

singularity exec --nv -B $VSC_DATA:$VSC_DATA -B $VSC_STAGING:$VSC_STAGING -B $DATA_DIR:$DATA_DIR -B $VSC_SCRATCH:$VSC_SCRATCH -B $VSC_HOME:$VSC_HOME -B $PWD docker://quay.io/biocontainers/iqtree:2.3.6--hdbdd923_0 \
    iqtree \
    -B 1000 \
    -m GTR+G4 \
    -T AUTO \
    -s $current_location \
    -pre ../../data/trees/matchingL-S/$(basename $current_location .msa)

# write some cluster metadata to the error and output files
printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE ***
************
'


