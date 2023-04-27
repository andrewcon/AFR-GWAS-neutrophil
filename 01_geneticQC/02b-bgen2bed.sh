#!/bin/bash
#SBATCH --job-name=xmerge
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:15:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8

module load apps/plink/1.90

plink \
--merge-list #.txt \
--out #.chr1-22
