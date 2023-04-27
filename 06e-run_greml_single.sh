#!/bin/bash

#SBATCH --job-name=gcta_greml_single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8

# Change working directory to submission directory
cd "${SLURM_SUBMIT_DIR}"

# Load modules required for runtime e.g.
GCTA_PATH=#
GRM_PATH=#
PHENO=#.pheno
COVAR=#.covar
QCOVAR=#.qcovar
OUT=#

echo $GCTA_PATH
echo $GRM_PATH
echo $PHENO
echo $COVAR
echo $QCOVAR
echo $OUT


# Run gcta64
  $GCTA_PATH \
  --reml \
  --grm $GRM_PATH \
  --pheno $PHENO \
  --covar $COVAR \
  --qcovar $QCOVAR \
  --thread-num 16 \
  --out $OUT/wod_single-data.chr1-22

# Run gcta64
  $GCTA_PATH \
  --reml \
  --grm $GRM_PATH \
  --pheno $PHENO \
  --covar $COVAR \
  --qcovar #.qcovar \
  --thread-num 16 \
  --out $OUT/wd_single-data.chr1-22

# Run gcta64
  $GCTA_PATH \
  --reml \
  --grm $GRM_PATH \
  --pheno #.pheno \
  --covar #.covar \
  --qcovar #.qcovar \
  --gxe #.gxe \
  --thread-num 16 \
  --out $OUT/GxE_single-data.chr1-22
