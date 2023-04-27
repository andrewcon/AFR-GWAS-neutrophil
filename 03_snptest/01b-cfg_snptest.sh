#!/bin/bash
#SBATCH --partition short
#SBATCH --job-name=snptest_array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22

# DO NOT RUN THIS FILE WITH sbatch 
# It has to be run through ./run_snptest.sh
# This is the config file for running SNPTEST on multiple clusters
# The main script run_snptest.sh executes a for loop
# with N = 18 clusters
# Each cluster has 22 chromosomes to run on, so there will be 
# 22 array jobs for each 18 serial jobs

# on compute node, change directory to 'submission directory':
cd $SLURM_SUBMIT_DIR

# Load snptest module
module load lib/openblas/0.3.15

# For each chromosome in cluster, perform GWAS
for i in ${SLURM_ARRAY_TASK_ID}; do
	if [ 10 -gt "$i" ]; then CHR=0$i; else CHR=$i; fi
	GEN=#.bgen
	SAM=#.txt
	OUTDIR=#
    
	# Check if $OUTDIR and $OUTDIR/log exist
	# If not, create them
    [ ! -d "$OUTDIR" ] && mkdir $OUTDIR
	[ ! -d "$OUTDIR/log" ] && mkdir $OUTDIR/log

	#/snptest_v2.5.6 \
		-data ${GEN} ${SAM} \
		-o ${OUTDIR}/chr${CHR}.txt \
		-log ${OUTDIR}/log/chr${CHR}.log \
		-missing_code NA \
		-use_raw_phenotypes \
		-pheno ${TRAIT} \
		-frequentist 1 \
		-method expected \
		-chunk 200 \
		-hwe
	
	# Remove hashes from snptest output so they are 
	# ready to be read with META or any other program
	#.sh ${OUTDIR}/chr${CHR}.txt
done

