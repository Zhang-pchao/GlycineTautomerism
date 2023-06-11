#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os, sys
import glob
import shutil
import random

import plumed
#reuse same plumed kernel, to avoid multiple warnings
PLUMED_KERNEL=plumed.Plumed()
from functools import partial
plumed.read_as_pandas = partial(plumed.read_as_pandas, kernel=PLUMED_KERNEL)

def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path) 

def get_input(file,info,foot_idx,subjob='pbs'):
    f = open(file, 'w')
    
    if subjob == 'pbs':
        f.write('#PBS -q normal2\n')
        f.write('#PBS -N deltaf\n')
        f.write('#PBS -l nodes=1:ppn=1\n')
        f.write('#PBS -l walltime=100:00:00\n')
        f.write('cd $PBS_O_WORKDIR\n')
        f.write('module load anaconda/anaconda.2020.02\n')
        f.write('source activate plumed_2.8.1\n')
    else: #slurm
        f.write('#!/bin/bash\n')
        f.write('#SBATCH -N 1\n')
        f.write('#SBATCH -n 1\n')
        f.write('#SBATCH -t 100:00:00\n')
        f.write('#SBATCH --mem-per-cpu=16GB\n')        
        f.write('#SBATCH --job-name=deltaf\n')
        f.write('cd $SLURM_SUBMIT_DIR\n')
        f.write('module load conda\n')
        f.write('conda activate plumed_v2.8.1\n')

    f.write('python ../FES_from_Reweighting4Gly_skipfooter.py ')
    f.write('-o ./fes-rew.dat ')
    f.write('-f ../../COLVAR_tmp ')
    f.write('--cv %s '%info['cv'])
    f.write('--min %f '%info['mins'])
    f.write('--max %f '%info['maxs'])
    f.write('--deltaFat %f '%info['deltaFat'])
    f.write('--sigma %f '%float(info['sigma']))
    f.write('--totalrow %d '%info['tot_time'])
    f.write('--skiprows %d '%info['skiprows'])
    f.write('--skipfoot %d '%info['skipfoot'][foot_idx])
    f.write('--blocks %d '%info['blocks'])
    f.write('--bin %d '%info['bins'])
    f.write('--temp %d \n\n'%info['temp'])
    f.close()

def get_sigma(file,cv):
    kernels = plumed.read_as_pandas(file)
    #print(kernels['sigma_%s'%cv])
    #print(list(kernels['sigma_%s'%cv])[-1])
    sigma = list(kernels['sigma_%s'%cv])[-1]
    return sigma

#if __name__ == '__main__':

path = "../"
info = {}
info["cv"]        = 's05'
info["mins"]      = -1.1
info["maxs"]      = 0.2
info["deltaFat"]  = -0.5
#info["sigma"]    = 0.02
cvname            = 's05'
info["sigma"]     = '%.4f'%get_sigma(os.path.join(path,'HILLS'),cvname)
info["skiprows"]  = 5e5
info["blocks"]    = 1
info["bins"]      = 200
info["temp"]      = 300 
info["tot_time"]  = 5e6
step              = 15
skipfoot  = []
for i in range(1,step+1):
    skipfoot.append((info["tot_time"]-info["skiprows"])-(info["tot_time"]-info["skiprows"])/step*i)
info["skipfoot"]      = skipfoot
print(info)

for i in range(step):
    name = str(10001+i)
    mkpath = os.path.join(path,'deltafes', name)
    mkdir(mkpath)
    subfile = os.path.join(mkpath,"q.sub")
    get_input(subfile,info,i,subjob='slurm')