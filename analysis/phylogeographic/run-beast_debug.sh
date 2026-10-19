#!/bin/bash
#SBATCH --clusters=genius
#SBATCH --account="lp_phylogeo_inf_gpu"
#SBATCH --mail-type="END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="joon.klaps@kuleuven.be"
#SBATCH --partition=gpu_p100_debug
#SBATCH --gres=gpu:3
#SBATCH --gpus-per-node=3
#SBATCH --nodes 1
#SBATCH --time="1:00:00"
#SBATCH --array=1-1

####################################
##### Run a beast analysis ######
####################################

# Accept XML file as command-line argument
XML_FILE="$1"
OUTDIR="$2"

if [ -z "$XML_FILE" ]; then
    echo "Error: No XML file provided"
    echo "Usage: sbatch run-beast.sh <path/to/xmlfile.xml> <outpath>"
    exit 1
fi

if [ ! -f "$XML_FILE" ]; then
    echo "Error: XML file not found: $XML_FILE"
    exit 1
fi

XML_FILE="$(realpath "$XML_FILE")"

if [ ! -d "$OUTDIR" ]; then
    echo "WARNING: Output directory not found: $OUTDIR"
    mkdir -p "$OUTDIR"
fi

STEM="$(basename ${XML_FILE%.xml})"

# Create a directory for this task
RUN_DIR="$OUTDIR/run_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

echo "Processing XML: $XML_FILE"
echo "Output stem: $STEM"
echo "Run directory: $PWD"

# Copy the XML file to the run directory to avoid mounting issues
cp "$XML_FILE" .
LOCAL_XML=$(basename "$XML_FILE")

abort() {
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

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

module load cluster/genius/gpu_p100_long
module load CUDA

singularity exec --nv \
    -B $VSC_DATA:$VSC_DATA \
    -B $VSC_STAGING:$VSC_STAGING \
    -B $VSC_SCRATCH:$VSC_SCRATCH \
    -B $VSC_HOME:$VSC_HOME \
    -B $PWD \
    docker://jklaps/beast-beagle-cuda:v10.5.0 \
    beast \
    -beagle_GPU \
    -beagle_order 1,2,3 \
    -overwrite \
    -working \
    -seed ${SLURM_ARRAY_TASK_ID}12345 \
    -save_stem "${STEM}.${SLURM_ARRAY_TASK_ID}" \
    "$LOCAL_XML" \
    > "${STEM}.${SLURM_ARRAY_TASK_ID}.out"

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE ***
************
'
