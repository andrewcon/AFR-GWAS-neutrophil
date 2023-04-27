#!/bin/bash
#SBATCH --partition short
#SBATCH --job-name=meta_array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22

# on compute node, change directory to 'submission directory':
cd $SLURM_SUBMIT_DIR

GWAS_DIR=#

# If running on the filtered SNPTEST chr GWAS, add a .filtered at the end of each .txt 
# Also remove .unfiltered from the --output
# For each chromosome perform META
i=${SLURM_ARRAY_TASK_ID}
if [ 10 -gt "$i" ]; then CHR=0$i; else CHR=$i; fi

WITH_DUFFY="$GWAS_DIR/res_wd_K1/chr${CHR}.txt $GWAS_DIR/res_wd_K2/chr${CHR}.txt  $GWAS_DIR/res_wd_K3/chr${CHR}.txt 
$GWAS_DIR/res_wd_K4/chr${CHR}.txt $GWAS_DIR/res_wd_K5/chr${CHR}.txt $GWAS_DIR/res_wd_K6/chr${CHR}.txt 
$GWAS_DIR/res_wd_K7/chr${CHR}.txt"
WITHOUT_DUFFY="$GWAS_DIR/res_wod_K1/chr${CHR}.txt $GWAS_DIR/res_wod_K2/chr${CHR}.txt  $GWAS_DIR/res_wod_K3/chr${CHR}.txt 
$GWAS_DIR/res_wod_K4/chr${CHR}.txt $GWAS_DIR/res_wod_K5/chr${CHR}.txt $GWAS_DIR/res_wod_K6/chr${CHR}.txt 
$GWAS_DIR/res_wod_K7/chr${CHR}.txt"

#/meta \
    --snptest \
    --method 1 \
    --use_info_col \
    --cohort $WITH_DUFFY \
    --threshold 0 \
    --output # \
    --snptest \
    --method 1 \
    --use_info_col \
    --cohort $WITHOUT_DUFFY \
    --threshold 0 \
    --output #.txt.unfiltered

WITH_DUFFY="$GWAS_DIR/res_wd_K1/chr${CHR}.txt.filtered $GWAS_DIR/res_wd_K2/chr${CHR}.txt.filtered  $GWAS_DIR/res_wd_K3/chr${CHR}.txt.filtered 
$GWAS_DIR/res_wd_K4/chr${CHR}.txt.filtered $GWAS_DIR/res_wd_K5/chr${CHR}.txt.filtered $GWAS_DIR/res_wd_K6/chr${CHR}.txt.filtered 
$GWAS_DIR/res_wd_K7/chr${CHR}.txt.filtered"
WITHOUT_DUFFY="$GWAS_DIR/res_wod_K1/chr${CHR}.txt.filtered $GWAS_DIR/res_wod_K2/chr${CHR}.txt.filtered  $GWAS_DIR/res_wod_K3/chr${CHR}.txt.filtered 
$GWAS_DIR/res_wod_K4/chr${CHR}.txt.filtered $GWAS_DIR/res_wod_K5/chr${CHR}.txt.filtered $GWAS_DIR/res_wod_K6/chr${CHR}.txt.filtered 
$GWAS_DIR/res_wod_K7/chr${CHR}.txt.filtered"

#/meta \
    --snptest \
    --method 1 \
    --use_info_col \
    --cohort $WITH_DUFFY \
    --threshold 0.3 \
    --output $GWAS_DIR/res_wd_meta/meta_chr${CHR}.txt
#/meta \
    --snptest \
    --method 1 \
    --use_info_col \
    --cohort $WITHOUT_DUFFY \
    --threshold 0.3 \
    --output $GWAS_DIR/res_wod_meta/meta_chr${CHR}.txt

