#!/bin/bash

# Load sample file to be used with SNPTEST
SAMPLE_PATH=#.txt

# Get column names corresponding to the 18 clusters the GWAS will be run on
COL_NAMES=$(head -n 1 $SAMPLE_PATH | cut -f 5-18)

# Assign the column names to an array TRAITS
declare -a TRAITS=($COL_NAMES)

# Run a sbatch serial job on each cluster,
# which will each have array jobs for each chromosome
# TRAIT var is current cluster and is 
# exported and passed onto sbatch with -V parameter
for x in ${TRAITS[@]}; do
    declare TRAIT=$x
    export "TRAIT"
    sbatch #.sh --export=TRAIT=$TRAIT
done
