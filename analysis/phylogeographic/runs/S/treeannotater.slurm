#!/bin/bash
#SBATCH --clusters=genius
#SBATCH --account="lp_phylogeo_inf_gpu"
#SBATCH --mail-type="END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="joon.klaps@kuleuven.be"
#SBATCH --partition=batch
#SBATCH --exclusive
#SBATCH --time="24:00:00"

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
module load Beast/10.5.0-GCC-12.3.0-CUDA-12.1.1

java -Xmx180g -cp /apps/leuven/rocky8/skylake/2023a/software/Beast/10.5.0-GCC-12.3.0-CUDA-12.1.1/lib/beast.jar dr.app.tools.treeannotator.TreeAnnotator LASV.S.combined1-5.trees LASV.S.combined1-5.tree

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE ***
************
'