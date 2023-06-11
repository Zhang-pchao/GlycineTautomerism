#!/usr/bin/env python
# coding: utf-8

import numpy as np
import plumed
import os
import sys
import MDAnalysis as mda

#reuse same plumed kernel, to avoid multiple warnings
PLUMED_KERNEL=plumed.Plumed()
from functools import partial
plumed.read_as_pandas = partial(plumed.read_as_pandas, kernel=PLUMED_KERNEL)

geo_path = "/data/HOME_BACKUP/pengchao/glycine/dpmd/geofile"
trj_path = "../"
trj_name = os.path.join(trj_path,"glycine_10.lammpstrj")
geo_name = os.path.join(geo_path,"CanonicalGlycine54H2O_opt_atomic.data")
cvs_file = plumed.read_as_pandas(os.path.join(trj_path,'COLVAR_tmp'))
cvs_type = ['s05','d05']
cvs_step = 1
trj_step = 1e1

flag   = "can2zwi"
#mysd  = [[0,0],[0,0.7],[0,1.4],[0,2.1],[0,2.8],[0,3.5]]
mysd   = []
sdsize = [0.125,0.5]
sdstep = [0.2,0.5]
sdstar = [0.0,0.0]
sdstop = [0.0,3.0]
mystep = 10
sstep  = (sdstop[0]-sdstar[0])/mystep
dstep  = (sdstop[1]-sdstar[1])/mystep
for i in range(mystep+1):
    mysd.append([sdstar[0]+i*sstep,sdstar[1]+i*dstep])
print('The range of cvs:\n',mysd)

u = mda.Universe(geo_name,trj_name,atom_style='id type x y z',format='LAMMPSDUMP')

step = int(trj_step/cvs_step)

cvs = cvs_file[::step]
sss = cvs[cvs_type[0]].to_numpy()
ddd = cvs[cvs_type[1]].to_numpy()

def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)

def lmp_index(trj_file):
    searchfile = open(trj_file, "r")
    lines = searchfile.readlines()
    searchfile.close()
    
    xyz_index = []
    i = 0
    for line in lines:
        if "ITEM: ATOMS id type x y z" in line:
            xyz_index.append(i) 
        i += 1
    
    start_index = []
    j = 0
    for line in lines:
        if "ITEM: TIMESTEP" in line:
            start_index.append(j)
        j += 1
            
    return xyz_index,start_index

def write_lmp_start(u,f,atoms,frame_idx,trj_step):
    """
    ITEM: TIMESTEP
    0
    ITEM: NUMBER OF ATOMS
    172
    ITEM: BOX BOUNDS pp pp pp
    0.0000000000000000e+00 1.2028000000000000e+01
    0.0000000000000000e+00 1.2028000000000000e+01
    0.0000000000000000e+00 1.2028000000000000e+01
    ITEM: ATOMS id type x y z
    ...
    """
    u.trajectory[frame_idx]
    f.write("ITEM: TIMESTEP\n")
    f.write("%d\n"%int(frame_idx*trj_step))
    f.write("ITEM: NUMBER OF ATOMS\n")
    f.write("%d\n"%atoms)
    f.write("ITEM: BOX BOUNDS pp pp pp\n")
    f.write("%16.8f%16.8f\n"%(0,u.dimensions[0]))
    f.write("%16.8f%16.8f\n"%(0,u.dimensions[1]))
    f.write("%16.8f%16.8f\n"%(0,u.dimensions[2]))
    f.write("ITEM: ATOMS id type x y z\n")

xyz_index,start_index = lmp_index(trj_name)

save_path = os.path.join(trj_path,flag)
mkdir(save_path)

atom_num = len(u.atoms)
for k in range(len(mysd)):
    shi = mysd[k][0]+sdsize[0]
    slo = mysd[k][0]-sdsize[0]
    dhi = mysd[k][1]+sdsize[1]
    dlo = mysd[k][1]-sdsize[1]    
    f = open(os.path.join(save_path, "%d_%s_%.2f_%.2f_%s_%.2f_%.2f_%s.lammpstrj"%(
        k+10001,cvs_type[0],slo,shi,cvs_type[1],dlo,dhi,flag)), 'w+')
    
    for i in range(len(sss)):
        if slo <= sss[i]<= shi and dlo <= ddd[i]<= dhi:
            u.trajectory[i]
            write_lmp_start(u,f,atom_num,i,trj_step)
            for j in range(atom_num):
                xyz = u.atoms[j].position
                f.write("%5d%5d%16.6f%16.6f%16.6f\n"%(j+1,int(u.atoms[j].type),xyz[0],xyz[1],xyz[2]))
    f.close()     