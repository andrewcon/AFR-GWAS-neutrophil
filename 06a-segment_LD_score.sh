#!/bin/bash

#SBATCH --job-name=gcta_ld_seg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22

i=${SLURM_ARRAY_TASK_ID}
if [ 10 -gt "$i" ]; then CHR=0$i; else CHR=$i; fi

GCTA_PATH=#
BED_PATH=#/data.chr${CHR}
OUT_PATH=#/step1-chr${CHR}

# Change working directory to submission directory
cd "${SLURM_SUBMIT_DIR}"

$GCTA_PATH --thread-num 8 --bfile $BED_PATH --ld-score-region 200 --out $OUT_PATH