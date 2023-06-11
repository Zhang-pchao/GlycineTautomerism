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

def get_1Dfes(file,cv,split):
    f = plumed.read_as_pandas(file)
    lcv = f[cv].to_numpy()
    lfes = f['file.free'].to_numpy()
    fesa = []
    fesb = []
    for i in range(len(lcv)):
        if lcv[i]<=split:
            fesa.append(lfes[i])
        else:
            fesb.append(lfes[i])
        
    return np.asarray(fesa),np.asarray(fesb)

def get_2Dfes(file,cv1,cv2,sp1,sp2):
    """
    sp1 = [cv1_low,cv1_high,cv2_low,cv2_high]
    sp2 = [cv1_low,cv1_high,cv2_low,cv2_high]
    """
    f = plumed.read_as_pandas(file)
    lcv1 = f[cv1].to_numpy()
    lcv2 = f[cv2].to_numpy()
    lfes = f['file.free'].to_numpy()
    fesa = []
    fesb = []
    for i in range(len(lcv1)):
        if sp1[0]<=lcv1[i]<=sp1[1] and sp1[2]<=lcv2[i]<=sp1[3]:
            fesa.append(lfes[i])
        if sp2[0]<=lcv1[i]<=sp2[1] and sp2[2]<=lcv2[i]<=sp2[3]:
            fesb.append(lfes[i])
        
    return np.asarray(fesa),np.asarray(fesb)

def delta_fes(fesa,fesb,kbt):
    fesA = -kbt*np.logaddexp.reduce(-1/kbt*fesa)
    fesB = -kbt*np.logaddexp.reduce(-1/kbt*fesb)
    deltaF = fesB-fesA
    return deltaF,fesB,fesA


path = "./"
Ffile= os.path.join(path,'fes-rew.dat')
Fout = open(os.path.join(path,'deltafes.dat'), 'w')
temp = 300
kbt  = temp*0.0083144621


fesa_1d, fesb_1d = get_1Dfes(Ffile,'s05',-0.5)
deltaFES_1d,fesB_1d,fesA_1d = delta_fes(fesa_1d,fesb_1d,kbt)
print("Delta 1D Free Energy = %4.4f - %4.4f = %4.4f kJ/mol\n"%(fesB_1d,fesA_1d,deltaFES_1d),file=Fout)

fesa_2d, fesb_2d = get_2Dfes(Ffile,'s05','d05',[-1.1,-0.5,-0.2,8.2],[-0.5,1.1,-0.2,8.2])
deltaFES_2d,fesB_2d,fesA_2d = delta_fes(fesa_2d,fesb_2d,kbt)
print("Delta 2D Free Energy = %4.4f - %4.4f = %4.4f kJ/mol\n"%(fesB_2d,fesA_2d,deltaFES_2d),file=Fout)

fesa_2d, fesb_2d = get_2Dfes(Ffile,'s05','d05',[-1.1,-0.5,1.5,8.0],[-0.5,0.2,1.5,6.0])
deltaFES_2d,fesB_2d,fesA_2d = delta_fes(fesa_2d,fesb_2d,kbt)
print("Delta 2D Free Energy = %4.4f - %4.4f = %4.4f kJ/mol\n"%(fesB_2d,fesA_2d,deltaFES_2d),file=Fout)

Fout.close()