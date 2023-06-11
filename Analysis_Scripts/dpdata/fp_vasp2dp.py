#!/usr/bin/env python
# coding: utf-8

import dpdata
import os

def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path) 

dirct = '/home/pengchao/cp2k/GlycineH2O/Glycine54H2O/SCAN_Dataset/'
index = '001'
save_path = os.path.join(dirct,'dataset',index)

all_poscar = None
all_outcar = None
for i in range(10007,10009):
    myposcar = os.path.join(dirct,index,str(i),'POSCAR')
    myoutcar = os.path.join(dirct,index,str(i),'OUTCAR')     
    if all_poscar is None:
        all_poscar = dpdata.System(myposcar, fmt = 'vasp/poscar')
        all_outcar = dpdata.LabeledSystem(myoutcar,  fmt = 'vasp/outcar')
    else:
        d_poscar = dpdata.System(myposcar, fmt = 'vasp/poscar')
        d_outcar = dpdata.LabeledSystem(myoutcar,  fmt = 'vasp/outcar')
        all_poscar.append(d_poscar)
        all_outcar.append(d_outcar)

mkdir(save_path)
all_outcar.to_deepmd_raw(save_path)
all_outcar.to_deepmd_npy(save_path)