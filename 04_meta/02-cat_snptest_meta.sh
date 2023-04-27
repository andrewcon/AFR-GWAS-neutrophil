#!/bin/bash
#SBATCH --partition short
#SBATCH --job-name=xcat_snptest
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-17

declare -a SNPTEST_DIR=(
    "" 
    "" 
    ""
)

i=${SLURM_ARRAY_TASK_ID}

cd ${SNPTEST_DIR[i]}

if [ "$i" -le "15" ]; then
    cat chr01.txt > chr1-22_all.txt
    for x in {02..22}; do
        file=chr$x.txt
        echo "Removing column names from $file..."
        sed '1d' $file >> chr1-22_all.txt
    done
    head -n 1 chr1-22_all.txt > chr1-22_significant.txt
    awk '$22 < 5e-8' chr1-22_all.txt >> chr1-22_significant.txt
elif [ "$i" -ge "16" && "$i" -le "17"]; then
    cat meta_chr01.txt > meta_chr1-22_all.txt
    for x in {02..22}; do
        file=meta_chr$x.txt
        echo "Removing column names from $file..."
        sed '1d' $file >> meta_chr1-22_all.txt
    done
    head -n 1 meta_chr1-22_all.txt > meta_chr1-22_significant.txt
    awk '$6 < 5e-8' meta_chr1-22_all.txt >> meta_chr1-22_significant.txt
else
    cat meta_chr01.txt.unfiltered > meta_chr1-22_all.txt.unfiltered
    for x in {02..22}; do
        file=meta_chr$x.txt.unfiltered
        echo "Removing column names from $file..."
        sed '1d' $file >> meta_chr1-22_all.txt.unfiltered
    done
    head -n 1 meta_chr1-22_all.txt.unfiltered > meta_chr1-22_significant.txt.unfiltered
    awk '$6 < 5e-8' meta_chr1-22_all.txt.unfiltered >> meta_chr1-22_significant.txt.unfiltered
fi

