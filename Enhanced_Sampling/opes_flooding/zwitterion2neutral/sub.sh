#!/bin/bash -l

seed=20240319

for (( i=101; i<151; i=i+1 )); do
    mkdir $i
    cp -r org/* $i 
    cd $i || exit 1
    sed -i "s/MYSEED/$(($seed+$i))/g" in.glycine
    sed -i "s/JOBNAME/n_$i/g" q.sh
    sbatch q.sh
    cd ..
done
