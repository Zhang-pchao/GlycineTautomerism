#!/usr/bin/env python
# coding: utf-8

import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import os
import MDAnalysis as mda
from MDAnalysis.lib import distances
import argparse

### Parser stuff ###
parser = argparse.ArgumentParser(description='get the RDF of glycine in water')
# files
parser.add_argument('--trjdir','-td',dest='tdirectory',type=str,default='./',help='the directory of trajectory files')
parser.add_argument('--savdir','-sd',dest='sdirectory',type=str,default='./',help='the directory of saving files')
parser.add_argument('--fstart',dest='fs',type=str,default='R',help='the name of starting frame')
parser.add_argument('--festop',dest='fe',type=str,default='P',help='the name of ending   frame')
# some easy parsing
args=parser.parse_args()

def find_min_bond_index(u,n_idx,o_idx1,o_idx2):
    # select hydrogen and oxygen atoms
    h_atoms = u.select_atoms("index %d"%(int(n_idx)-1))
    o_atoms = u.select_atoms("index %d %d"%(int(o_idx1)-1,int(o_idx2)-1))

    # initialize minimum distance and closest oxygen index
    min_distance = np.inf
    closest_o_index = -1
    
    # loop over hydrogen and oxygen atoms
    for h_atom in h_atoms:
        for j, o_atom in enumerate(o_atoms):
            # calculate distance between hydrogen and oxygen atoms
            distance = distances.calc_bonds(h_atom.position, o_atom.position,box=u.dimensions)
            
            # update closest oxygen index and minimum distance if current oxygen is closer
            if distance < min_distance:
                closest_o_index = j
                min_distance = distance
    
    # print the index of the oxygen atom with the shortest O-H bond distance
    #print(f"min bond index and length: {o_atoms[closest_o_index].index}, {min_distance}")
    Oh_index = o_atoms[closest_o_index].index
    O_index = int(o_idx1-1+o_idx2-1-Oh_index)
    return Oh_index,O_index

####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/glycine/dpmd/geofile"
trj_path    = args.tdirectory
save_path   = args.sdirectory
data_geo    = "Glycine128H2O_opt_atomic.data"
trj_skip    = 0
mda_step    = 1
NOOdict     = {'N': 385,'O1': 391,'O2': 393,'C1': 387,'C2': 388} #serial in lammps .data file
pathindex   = [args.fs,args.fe] #'cationic','anionic'
####################change above####################
data_geo    = os.path.join(geo_path,data_geo)

trjs = []
for filename in os.listdir(trj_path):
    if filename.endswith('.lammpstrj'):
        trjs.append(filename)

trj_file    = os.path.join(trj_path,trjs[0])
u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')
u.trajectory[1]
c1_idx,c2_idx = find_min_bond_index(u,NOOdict['N'],NOOdict['C1'],NOOdict['C2'])

All_BAD = np.array([])
for k in range(len(trjs)):
    trj_file    = os.path.join(trj_path,trjs[k])
    u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')
    len_trj  = len(u.trajectory)
    #trj_info = [len_trj,trj_skip,mda_step]
    BAD = np.array([])
    for j in range(len_trj):       
        u.trajectory[j]
        o1_idx,o2_idx = find_min_bond_index(u,NOOdict['N'],NOOdict['O1'],NOOdict['O2'])       
        O1 = u.select_atoms("index %d"%o1_idx)[0].position
        O2 = u.select_atoms("index %d"%o2_idx)[0].position
        N0 = u.select_atoms("index %d"%(int(NOOdict['N'])-1))[0].position
        C1 = u.select_atoms("index %d"%c1_idx)[0].position
        C2 = u.select_atoms("index %d"%c2_idx)[0].position
        #print(O1,O2,N0,C1,C2)
        B_N0O1     = distances.calc_bonds(N0, O1, box=u.dimensions)
        B_N0O2     = distances.calc_bonds(N0, O2, box=u.dimensions)
        A_N0C1C2   = distances.calc_angles(N0, C1, C2, box=u.dimensions)
        A_C1C2O1   = distances.calc_angles(C1, C2, O1, box=u.dimensions)
        A_C1C2O2   = distances.calc_angles(C1, C2, O2, box=u.dimensions)
        A_O1C2O2   = distances.calc_angles(O1, C2, O2, box=u.dimensions)
        D_N0C1C2O1 = abs(distances.calc_dihedrals(N0, C1, C2, O1, box=u.dimensions))
        D_N0C1C2O2 = abs(distances.calc_dihedrals(N0, C1, C2, O2, box=u.dimensions))

        b = np.array([B_N0O1,B_N0O2])
        ad = np.array([A_N0C1C2,A_C1C2O1,A_C1C2O2,A_O1C2O2,D_N0C1C2O1,D_N0C1C2O2])#/np.pi*180
        bad = np.concatenate((b,ad),axis=None)
        BAD = np.concatenate((BAD,bad),axis=0)
    a_BAD = np.average(BAD.reshape((-1,bad.shape[0])),axis=0)
    All_BAD = np.concatenate((All_BAD,a_BAD),axis=0)

All_BAD = All_BAD.reshape((len(trjs),-1))

BADfile = open(os.path.join(save_path, "BAD_%s_%s_step%d.txt"%(pathindex[0],pathindex[1],mda_step)), 'w+')
BADfile.write("#%12s%12s%12s%12s%12s%12s%12s%12s%12s\n"%("Index","B_N0O1", "B_N0O2",    
                                                         "A_N0C1C2","A_C1C2O1",  
                                                         "A_C1C2O2","A_O1C2O2",  
                                                         "D_N0C1C2O1","D_N0C1C2O2"))
Index=[]
for k in range(len(trjs)):
    Index.append(k+1)
    BADfile.write(' %12d%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f\n' %(k+1,
                    All_BAD[k][0],All_BAD[k][1],All_BAD[k][2],All_BAD[k][3],
                    All_BAD[k][4],All_BAD[k][5],All_BAD[k][6],All_BAD[k][7] ))
BADfile.close()
colors = [(31 ,59 ,115), 
          (47 ,146,148), 
          (80 ,178,141),
          (167,214,85 ),
          (255,224,62 ),
          (255,169,85 ),
          (225,110,65 ),          
          (189, 48,47 )]
colors=list(tuple(i/255 for i in color) for color in colors)
colors.reverse()
#set bigger font sizes
SMALL_SIZE = 13
MEDIUM_SIZE = 14
BIG_SIZE = 17
plt.rc('font', size=SMALL_SIZE)        # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)   # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE-2)  # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)   # fontsize of the figure title

def BAD_Fig(x,y,name):
    fig = plt.figure(figsize=(8,6), dpi=150, facecolor='white')
    ax = fig.add_subplot(1, 1, 1)
    for i in range(y.shape[0]):
        ax.plot(x, y[i], lw=2, c=colors[i],label="%s"%name[i])
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.legend()
    ax.set_ylabel("Bond/Angle/Dihedral")
    ax.set_xlabel("Index")
    fig.savefig(os.path.join(save_path, "BAD.png"), dpi=600, bbox_inches='tight')

BAD_Fig(Index,All_BAD.T,
             ["B_N0O1", "B_N0O2", "A_N0C1C2","A_C1C2O1","A_C1C2O2","A_O1C2O2", "D_N0C1C2O1","D_N0C1C2O2"])
