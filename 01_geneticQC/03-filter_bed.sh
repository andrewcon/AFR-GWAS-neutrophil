#!/bin/bash

#SBATCH --job-name=xmakebfile
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

# on compute node, change directory to 'submission directory':
cd $SLURM_SUBMIT_DIR

module load apps/plink/1.90

DIRECTLY_GENO_FILES=#
OUT_DIR=#
KEEP_FILE=#.txt

# Run plink
plink \
    --threads 8 \
    --memory 65536 \
    --bfile $DIRECTLY_GENO_FILES \
    --keep $KEEP_FILE \
    --make-bed \
    --out $OUT_DIR
