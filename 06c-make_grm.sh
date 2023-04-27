#!/bin/bash

#SBATCH --job-name=gcta_make_grm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:15:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-87

# Change working directory to submission directory
cd "${SLURM_SUBMIT_DIR}"

i=${SLURM_ARRAY_TASK_ID}
if [[ "$i" -ge "0" && "$i" -le "3" ]]; then CHR=01;
    elif [[ "$i" -gt "3" && "$i" -le "7" ]]; then CHR=02;
        elif [[ "$i" -gt "7" && "$i" -le "11" ]]; then CHR=03;
            elif [[ "$i" -gt "11" && "$i" -le "15" ]]; then CHR=04;
                elif [[ "$i" -gt "15" && "$i" -le "19" ]]; then CHR=05;
                    elif [[ "$i" -gt "19" && "$i" -le "23" ]]; then CHR=06;
                        elif [[ "$i" -gt "23" && "$i" -le "27" ]]; then CHR=07;
                            elif [[ "$i" -gt "27" && "$i" -le "31" ]]; then CHR=08;
                                elif [[ "$i" -gt "31" && "$i" -le "35" ]]; then CHR=09;
                                    elif [[ "$i" -gt "35" && "$i" -le "39" ]]; then CHR=10;
                                        elif [[ "$i" -gt "39" && "$i" -le "43" ]]; then CHR=11;
                                            elif [[ "$i" -gt "43" && "$i" -le "47" ]]; then CHR=12;
                                                elif [[ "$i" -gt "47" && "$i" -le "51" ]]; then CHR=13;
                                                    elif [[ "$i" -gt "51" && "$i" -le "55" ]]; then CHR=14;
                                                        elif [[ "$i" -gt "55" && "$i" -le "59" ]]; then CHR=15;
                                                            elif [[ "$i" -gt "59" && "$i" -le "63" ]]; then CHR=16;
                                                                elif [[ "$i" -gt "63" && "$i" -le "67" ]]; then CHR=17;
                                                                    elif [[ "$i" -gt "67" && "$i" -le "71" ]]; then CHR=18;
                                                                        elif [[ "$i" -gt "71" && "$i" -le "75" ]]; then CHR=19;
                                                                            elif [[ "$i" -gt "75" && "$i" -le "79" ]]; then CHR=20;
                                                                                elif [[ "$i" -gt "79" && "$i" -le "83" ]]; then CHR=21;
                                                                                    elif [[ "$i" -gt "83" && "$i" -le "87" ]]; then CHR=22; fi

GCTA_PATH=#
BED_PATH=#

shopt -s nullglob
GRP_CHR_LIST=(/#.txt)
GRP=$(( $i % 4 + 1))

OUT=#

echo $GCTA_PATH
echo $GRM_CHR_LIST
echo $OUT

# Run gcta64
$GCTA_PATH \
--bfile $BED_PATH \
--extract "${GRP_CHR_LIST[$i]}" \
--make-grm \
--thread-num 4 \
--out $OUT
