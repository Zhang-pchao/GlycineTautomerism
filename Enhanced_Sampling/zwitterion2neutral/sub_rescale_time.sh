#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 500:00:00
#SBATCH --job-name=rescale_time

cd $SLURM_SUBMIT_DIR

module load conda
conda activate dflow_ir

for (( i=151; i<161; i=i+1 )); do
    cd $i || exit 1
	python ../rescale_time.py
    cd ..
done
