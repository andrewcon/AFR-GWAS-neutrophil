# On compute node, change directory to 'submission directory':
cd #

# Load module
module load languages/python-3.7.7-anaconda-2020.20-gprMAX

# Activate ldsc dependencies
source activate ldsc

plink \
    --bfile # \
    --cm-map genetic_map_chr@_combined_b37.txt \
    --make-bed \
    --out #

# Compile the LD score file
python ldsc.py \
	--bfile fileset_with_cms \
	--l2 \
	--ld-wind-cm 1 \
	--out ldsc