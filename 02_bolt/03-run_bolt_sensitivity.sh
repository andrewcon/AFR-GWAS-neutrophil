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
--geneticMapFile=#/genetic_map_hg19_withX.txt.gz \
--phenoFile=#.txt \
--phenoCol=nc_log \
--covarFile=#.txt \
--covarMaxLevels=30 \
--covarCol=neutrophil_device \
--covarCol=sample_year \
--covarCol=sample_month \
--covarCol=cob_un_regions \
--covarCol=Kpop \
--covarCol=smoking_status \
--covarCol=drinker_status \
--covarCol=menstrual_status \
--covarCol=assessment_center \
--qCovarCol=sample_day \
--qCovarCol=sample_day_minutes \
--qCovarCol=age \
--qCovarCol=BMI \
--qCovarCol=PC{1:100} \
--lmm \
--LDscoresFile=#.ldscore.gz \
--LDscoresMatchBp \
--numThreads=16 \
--remove=#.txt \
--verboseStats \
--statsFileBgenSnps=#.bgen.gz \
--statsFile=#.gz
