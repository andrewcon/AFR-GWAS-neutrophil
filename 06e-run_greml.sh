#!/bin/bash

#SBATCH --job-name=gcta_greml
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-3

# Change working directory to submission directory
cd "${SLURM_SUBMIT_DIR}"

i=${SLURM_ARRAY_TASK_ID}

# Load modules required for runtime e.g.
GCTA_PATH=#
GRM_PATH=#.txt
PHENO=#.pheno
COVAR=#.covar
QCOVAR=#.qcovar
OUT=#

echo $GCTA_PATH
echo $GRM_PATH
echo $PHENO
echo $COVAR
echo $QCOVAR
echo $GxE
echo $OUT

if [[ "$i" -ge "1" && "$i" -le "2" ]]
then
  # Run gcta64
    $GCTA_PATH \
    --reml \
    --reml-no-constrain \
    --mgrm $GRM_PATH \
    --pheno $PHENO \
    --covar $COVAR \
    --qcovar $QCOVAR \
    --thread-num 16 \
    --out $OUT/${i}-data.chr1-22
else
    GxE=#.gxe
  # Run gcta64
    $GCTA_PATH \
    --reml \
    --reml-no-constrain \
    --mgrm $GRM_PATH \
    --pheno $PHENO \
    --covar $COVAR \
    --qcovar $QCOVAR \
    --gxe $GxE \
    --thread-num 16 \
    --out $OUT/${i}-data.chr1-22
fi

