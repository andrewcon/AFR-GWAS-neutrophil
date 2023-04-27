#!/bin/bash

declare -a SNPTEST_DIR=(
    "#" 
    "#" 
    "#"
)

for dir in "${SNPTEST_DIR[@]}"; do
    cd $dir
    for i in {1..22}; do
        if [ 10 -gt "$i" ]; then CHR=0$i; else CHR=$i; fi
        file=${dir}/chr${CHR}.txt
        echo "Removing hashes from $file..."
        sed -i '/^#/d' $file
    done
done

