#!/bin/bash
#SBATCH --partition short
#SBATCH --job-name=gcta_make_grm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=03:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8

# Change working directory to submission directory
cd "${SLURM_SUBMIT_DIR}"

GCTA_PATH=#
BED_PATH=#

OUT=#

echo $GCTA_PATH
echo $GRM_CHR_LIST
echo $OUT

# Run gcta64
$GCTA_PATH \
--bfile $BED_PATH \
--make-grm \
--thread-num 16 \
--out $OUT
