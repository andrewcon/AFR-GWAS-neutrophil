#!/bin/bash
#SBATCH --job-name=xcat_meta_malgen
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:05:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --array=0-3

declare -a SNPTEST_DIR=(
    "#" 
    "#" 
    "#" 
    "#"
)

i=${SLURM_ARRAY_TASK_ID}

cd ${SNPTEST_DIR[i]}

if [ "$i" -eq "0" ]; then
    head -n 1 meta_overall.txt > meta_overall_significant.txt
    awk '$12 < 5e-8' meta_overall.txt >> meta_overall_significant.txt
    head -n 1 meta_overall.txt > meta_overall_significant_5e-5.txt
    awk '$12 < 5e-5' meta_overall.txt >> meta_overall_significant_5e-5.txt
elif [ "$i" -eq "1"]; then
    head -n 1 meta_cm.txt > meta_cm_significant.txt
    awk '$12 < 5e-8' meta_cm.txt >> meta_cm_significant.txt
    head -n 1 meta_cm.txt > meta_cm_significant_5e-5.txt
    awk '$12 < 5e-5' meta_cm.txt >> meta_cm_significant_5e-5.txt
elif [ "$i" -eq "2"]; then
    head -n 1 meta_sma.txt > meta_sma_significant.txt
    awk '$12 < 5e-8' meta_sma.txt >> meta_sma_significant.txt
    head -n 1 meta_sma.txt > meta_sma_significant_5e-5.txt
    awk '$12 < 5e-5' meta_sma.txt >> meta_sma_significant_5e-5.txt
else
    head -n 1 meta_other.txt > meta_other_significant.txt
    awk '$12 < 5e-8' meta_other.txt >> meta_other_significant.txt
    head -n 1 meta_other.txt > meta_other_significant_5e-5.txt
    awk '$12 < 5e-5' meta_other.txt >> meta_other_significant_5e-5.txt
fi

