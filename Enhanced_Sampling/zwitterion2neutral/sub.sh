#!/bin/bash -l

seed=20240319

#for (( i=101; i<111; i=i+1 )); do
#for (( i=111; i<121; i=i+1 )); do
#for (( i=121; i<131; i=i+1 )); do
#for (( i=131; i<141; i=i+1 )); do
#for (( i=141; i<151; i=i+1 )); do
for (( i=151; i<161; i=i+1 )); do
    #mkdir $i
    #cp -r org/* $i 
	#cp ../anion/$i/in.glycine $i 		
    cd $i || exit 1
    #sed -i "s/MYSEED/$(($seed+$i))/g" in.glycine
    #sed -i "s/JOBNAME/n_$i/g" q.sh
    sbatch q.sh
    cd ..
done
