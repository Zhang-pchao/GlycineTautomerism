#!/bin/bash -l
#PBS -q normal2
#PBS -N abacus
#PBS -l nodes=1:ppn=20
#PBS -l walltime=1000:00:00

cd $PBS_O_WORKDIR

module load anaconda/anaconda.2020.02
source activate abacus_3.0.5
ulimit -s unlimited
export OMP_NUM_THREADS=1

mpirun -np 20 abacus