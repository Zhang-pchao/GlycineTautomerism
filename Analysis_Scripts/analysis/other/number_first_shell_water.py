#!/usr/bin/env python
# coding: utf-8

import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from matplotlib import cm
import os
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.rdf import InterRDF,InterRDF_s
from MDAnalysis.lib import distances
import argparse

### Parser stuff ###
parser = argparse.ArgumentParser(description='get the RDF of glycine in water')
# files
parser.add_argument('--trjdir','-td',dest='tdirectory',type=str,default='./',help='the directory of trajectory files')
parser.add_argument('--savdir','-sd',dest='sdirectory',type=str,default='./',help='the directory of saving files')
# some easy parsing
args=parser.parse_args()

def find_water_num(u,noo_atom,shell_atom,cutoff):
    shell_idx,shell_idx_all = [],[]
    for atom1 in noo_atom:
        for j, atom2 in enumerate(shell_atom):
            distance = distances.calc_bonds(atom1.position, atom2.position,box=u.dimensions)
            if distance < cutoff:
                shell_idx.append(atom2.index)
    [shell_idx_all.append(i) for i in shell_idx if i not in shell_idx_all]
    #print(shell_idx_all)
    return len(shell_idx_all)

def get_avg_std(x):
    avg = np.average(x, axis=0)
    std = np.std(x, axis=0, ddof=1)
    return avg,std


####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/glycine/dpmd/geofile"
#myroot      ='/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/020_1e-6_add_water/008_gly128H2O_sd_B35_compress_f'
#trj_path    = os.path.join(myroot,'0-30ns/trj_file_2/ani2can_2/10001_0_10_ns')
#save_path   = os.path.join(myroot,'0-30ns/HBPathTrj_2/ani2can_2/10001_0_10_ns')
trj_path    = args.tdirectory
save_path   = args.sdirectory
data_geo    = "Glycine128H2O_opt_atomic.data"
trj_skip    = 0
mda_step    = 1
cutoff      = 3.6 # find the number of water in the cutoff 3.6 A
NOOdict     = {'N': 385,'O1': 391,'O2': 393,'C1': 387,'C2': 388} #serial in lammps .data file
####################change above####################
data_geo    = os.path.join(geo_path,data_geo)

trjs = []
for filename in os.listdir(trj_path):
    if filename.endswith('.lammpstrj'):
        trjs.append(filename)

File = open(os.path.join(save_path, "WaterNum_step%d.txt"%(mda_step)), 'w+')
File.write("#%12s%12s%12s\n"%("Index","avg","std"))
avg,std = [],[]

Index = []
for k in range(len(trjs)):
#for k in [0,-2]:    
    trj_file    = os.path.join(trj_path,trjs[k])
    u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')
    len_trj  = len(u.trajectory)
    trj_info = [len_trj,trj_skip,mda_step]
    num_list = []
    for i in range(trj_info[1],trj_info[0],trj_info[2]):
        u.trajectory[i]
        shell_atom = u.select_atoms("type 2 and not (index %d or index %d or index %d or index %d or index %d)"%(
                    NOOdict['N']-1,NOOdict['O1']-1,NOOdict['O2']-1,NOOdict['C1']-1,NOOdict['C2']-1))
        noo_atom   = u.select_atoms("index %d or index %d or index %d or index %d or index %d"%(
                    NOOdict['N']-1,NOOdict['O1']-1,NOOdict['O2']-1,NOOdict['C1']-1,NOOdict['C2']-1))
        o_num_in_shell = find_water_num(u,noo_atom,shell_atom,cutoff)
        num_list.append(o_num_in_shell)
        #print(o_num_in_shell)    
    Index.append(k+1)
    #print(num_list)
    avg,std = get_avg_std(num_list)
    File.write(' %12d%12.4f%12.4f\n' %(k,avg,std))
File.close()
