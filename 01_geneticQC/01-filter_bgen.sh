#!/bin/bash

#SBATCH --job-name=xmakebgen
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22

#Load relevant module(s)
module load apps/qctool/2.0.7

# on compute node, change directory to 'submission directory':
cd $SLURM_SUBMIT_DIR

BGEN_PATH=#


## array id
CHR=${SLURM_ARRAY_TASK_ID}

if ((${#CHR} == 1))
then
        CHR=0${CHR}
fi

qctool \
-threads 8 \
-incl-samples #.txt \
-g $BGEN_PATH/data.chr${CHR}.bgen \
-s $BGEN_PATH/data.chr${CHR}.sample \
-bgen-bits 8 \
-og #${CHR}.bgen \
-os #${CHR}.sample
