#!/bin/bash
#SBATCH --job-name=makebed
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22

# On compute node, change directory to 'submission directory':
cd $SLURM_SUBMIT_DIR

module load apps/plink/2.0.0
module load apps/qctool/2.0.7

BGEN_PATH=#
BED_PATH=#

# Run plink
for i in ${SLURM_ARRAY_TASK_ID}; do
	if [ 10 -gt "$i" ]; then CHR=0$i; else CHR=$i; fi

    if [ -f "#${CHR}.snpstats" ]; then
        echo "#${CHR}.snpstats exists, skipping qctool."
    else
        qctool \
        -g $BGEN_PATH/data.chr${CHR}.bgen \
        -s $BGEN_PATH/data.chr${CHR}.sample \
        -snp-stats \
        -snp-stats-columns "info" \
        -osnp #${CHR}.snpstats

        awk '$9 > 0.3' #${CHR}.snpstats | tail -n +2 | cut -f 1 > #${CHR}.snpid.filterInfo
        awk '$9 > 0.3' #${CHR}.snpstats | tail -n +2 | cut -f 2 > #${CHR}.rsid.filterInfo
    fi

    #/plink2 \
    --threads 4 \
    --memory 8192 \
    --keep #.txt \
    --rm-dup 'force-first' \
    --hwe 1e-10 'midp' 'keep-fewhet' \
    --mac 17 \
    --geno \
    --extract #${CHR}.rsid.filterInfo \
    --bgen $BGEN_PATH/data.chr${CHR}.bgen 'ref-first' \
    --sample $BGEN_PATH/data.chr${CHR}.sample \
    --make-bed \
    --out #${CHR}

done
