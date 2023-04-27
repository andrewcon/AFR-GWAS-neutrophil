#!/bin/bash
#SBATCH --partition short
#SBATCH --job-name=xbolt_lmm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

# on compute node, change directory to 'submission directory':
cd $SLURM_SUBMIT_DIR

#/bolt \
--bfile=# \
--bgenFile=#.bgen \
--bgenFile=#.bgen \
--bgenMinMAF=0 \
--sampleFile=#.sample \
--geneticMapFile=#.txt.gz \
--phenoFile=#.txt \
--phenoCol=nc_log \
--covarFile=#.txt \
--covarMaxLevels=30 \
--covarCol=sex \
--covarCol=assessment_center \
--covarCol=neutrophil_device \
--covarCol=sample_year \
--covarCol=sample_month \
--qCovarCol=sample_day \
--qCovarCol=sample_day_minutes \
--qCovarCol=age \
--qCovarCol=PC{1:100} \
--lmm \
--LDscoresFile=#.ldscore.gz \
--LDscoresMatchBp \
--numThreads=16 \
--remove=#.txt \
--verboseStats \
--statsFileBgenSnps=#.bgen.gz \
--statsFile=#.gz
