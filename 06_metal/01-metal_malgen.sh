#!/bin/bash
#SBATCH --job-name=meta_array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

# on compute node, change directory to 'submission directory':
cd $SLURM_SUBMIT_DIR

module load apps/snptest/2.5.4

# If running on the filtered SNPTEST chr GWAS, add a .filtered at the end of each .txt 
# Also remove .unfiltered from the --output
# For each chromosome perform META
PATH/metal #.txt
PATH/metal #.txt
PATH/metal #.txt
PATH/metal #.txt