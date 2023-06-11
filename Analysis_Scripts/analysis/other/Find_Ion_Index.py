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

###NOTE###
#check element order: should be O H H O H H ... (in .lammpstrj file)

###CHANGE BELOW###
geo_path = "/data/HOME_BACKUP/pengchao/glycine/dpmd/geofile"
trj_path = '../'
save_path= './'
trj_name = os.path.join(trj_path,"glycine_10.lammpstrj")
geo_name = os.path.join(geo_path,"127w_1oh.data")
cvs_file = plumed.read_as_pandas(os.path.join(trj_path,'COLVAR'))
ion_type = 'im' # cv: 'im' for OH minus, 'ip' for H3O positive
tcs_type = 'tc' # cv: 'tc' for total change
o_type   = '2'  # original O type in .lammpstrj file
ionOtype = 9    # new type for O in Ion 
cvs_step = 1    # COLVAR file dump
trj_step = 1e1  # .lammpstrj file dump
###CHANGE ABOVE###

u = mda.Universe(geo_name,trj_name, atom_style='id type x y z',format='LAMMPSDUMP')
step = int(trj_step/cvs_step)
#print(cvs_file[::step])

cvs = cvs_file[::step]
tcs = cvs[tcs_type].to_numpy()
ion = cvs[ion_type].to_numpy()

def isInt(x):
    if x%1 == 0:
        return True
    else:
        return False

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
    385
    ITEM: BOX BOUNDS pp pp pp
    0.0000000000000000e+00 1.5526350000000001e+01
    0.0000000000000000e+00 1.5526350000000001e+01
    0.0000000000000000e+00 1.5526350000000001e+01
    ITEM: ATOMS id type x y z
    1 2 4.51246 8.45677 1.20054
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

ion_idx = []
for i in range(len(tcs)):
    u.trajectory[i]
    if 0.999 <= tcs[i]<= 1.001:
        if isInt(ion[i]):
            ion_idx.append(int(ion[i]*3+1)) # NOTE: *3 for O H H (.lammpstrj), +1 for serial in .data file
        else:
            ion_idx.append(ion_idx[-1])
    else:
        ion_idx.append(ion_idx[-1])
        
xyz_index,start_index = lmp_index(trj_name)

atom_num = len(u.atoms)
f = open(os.path.join(save_path, "add_Ion_O_type_9.lammpstrj"), 'w+')
for i in range(len(tcs)):
    write_lmp_start(u,f,atom_num+1,i,trj_step) #add one atom
    u.trajectory[i]
    for j in range(atom_num):
        xyz = u.atoms[j].position
        f.write("%5d%5d%16.6f%16.6f%16.6f\n"%(j+1,int(u.atoms[j].type),xyz[0],xyz[1],xyz[2]))
    for j in range(atom_num):           
        if j == ion_idx[i]-1:
            if u.atoms[ion_idx[i]-1].type == o_type:# O='2', -1 for index in mda
                xyz = u.atoms[j].position
                f.write("%5d%5d%16.6f%16.6f%16.6f\n"%(atom_num+1,ionOtype,xyz[0],xyz[1],xyz[2]))
            else:
                print('The Ion index %d (.lammpstrj) is not O type in frame %d'%(ion_idx[i],i*trj_step))
f.close()           

f = open(os.path.join(save_path, "only_Ion_O_type_9.lammpstrj"), 'w+')
for i in range(len(tcs)):
    write_lmp_start(u,f,1,i,trj_step)
    u.trajectory[i]
    for j in range(atom_num):
        xyz = u.atoms[j].position
        if j == ion_idx[i]-1:
            if u.atoms[ion_idx[i]-1].type == o_type:# O='2', -1 for index in mda
                f.write("%5d%5d%16.6f%16.6f%16.6f\n"%(1,ionOtype,xyz[0],xyz[1],xyz[2]))
            else:
                print('The Ion index %d (.lammpstrj) is not O type in frame %d'%(ion_idx[i],i*trj_step))
f.close()            
