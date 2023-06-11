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
parser.add_argument('--lmptrj',dest='tj',type=str,default='10001',help='lammps trj file')
# some easy parsing
args=parser.parse_args()

def find_min_bond_index(u,n_idx,o_idx1,o_idx2):
    # select hydrogen and oxygen atoms
    h_atoms = u.select_atoms(n_idx)
    o_atoms = u.select_atoms("index %d %d"%(int(o_idx1)-1,int(o_idx2)-1))

    # initialize minimum distance and closest oxygen index
    min_distance = np.inf
    closest_o_index = -1
    
    # loop over hydrogen and oxygen atoms
    for h_atom in h_atoms:
        #print(h_atom.index)
        for j, o_atom in enumerate(o_atoms):
            # calculate distance between hydrogen and oxygen atoms
            distance = distances.calc_bonds(h_atom.position, o_atom.position,box=u.dimensions)
            
            # update closest oxygen index and minimum distance if current oxygen is closer
            if distance < min_distance:
                closest_o_index = j
                closest_h_index = h_atom.index
                min_distance = distance
    
    # print the index of the oxygen atom with the shortest O-H bond distance
    #print(f"min bond index and length: {o_atoms[closest_o_index].index}, {min_distance}")
    Oh_index = o_atoms[closest_o_index].index
    O_index = int(o_idx1-1+o_idx2-1-Oh_index)
    return Oh_index,O_index,min_distance,closest_h_index

def get_avg_std(x):
    avg = np.average(x, axis=0)
    std = np.std(x, axis=0, ddof=1)
    return avg,std

####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/glycine/dpmd/geofile"
trj_path    = args.tdirectory
save_path   = args.sdirectory
data_geo    = "Glycine128H2O_opt_atomic.data"
trj_skip    = 0
mda_step    = 1
NOOdict     = {'N': 385,'O1': 391,'O2': 393,'C1': 387,'C2': 388} #serial in lammps .data file
pathindex   = [args.fs,args.fe] #'cationic','anionic'
oh_cutoff   = 1.2
trjs        = ['%s.lammpstrj'%args.tj] 
####################change above####################
data_geo    = os.path.join(geo_path,data_geo)

trj_file    = os.path.join(trj_path,trjs[0])
u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')
u.trajectory[1]
c1_idx,c2_idx,d_nc1,_idx = find_min_bond_index(u,"index %d"%(int(NOOdict['N'])-1),NOOdict['C1'],NOOdict['C2'])

#All_BA = np.array([])
#for k in range(len(trjs)):
for k in [0]:    
    trj_file    = os.path.join(trj_path,trjs[k])
    u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')
    len_trj  = len(u.trajectory)
    #trj_info = [len_trj,trj_skip,mda_step]
    BA = np.array([])
    for j in range(len_trj):       
        u.trajectory[j]
        o1_idx,o2_idx,d_ho1,h_idx = find_min_bond_index(u,"type 1",NOOdict['O1'],NOOdict['O2']) 
        if d_ho1 <= oh_cutoff:
            #o1_idx,o2_idx,d_no1,_idx = find_min_bond_index(u,"index %d"%(int(NOOdict['N'])-1),NOOdict['O1'],NOOdict['O2'])               
            H0 = u.select_atoms("index %d"%h_idx)[0].position
            O1 = u.select_atoms("index %d"%o1_idx)[0].position
            O2 = u.select_atoms("index %d"%o2_idx)[0].position
            N0 = u.select_atoms("index %d"%(int(NOOdict['N'])-1))[0].position
            C1 = u.select_atoms("index %d"%c1_idx)[0].position
            C2 = u.select_atoms("index %d"%c2_idx)[0].position
            #print(O1,O2,N0,C1,C2)
            B_N0O1     = distances.calc_bonds(N0, O1, box=u.dimensions)
            A_C1O1H0   = distances.calc_angles(C1, O1, H0, box=u.dimensions)
            #D_N0C1C2O1 = abs(distances.calc_dihedrals(N0, C1, C2, O1, box=u.dimensions))
            ba = np.array([B_N0O1,A_C1O1H0])#/np.pi*180
            BA = np.concatenate((BA,ba),axis=0)
    #a_BAD = np.average(BAD.reshape((-1,bad.shape[0])),axis=0)
    #All_BAD = np.concatenate((All_BAD,a_BAD),axis=0)

BA = BA.reshape((-1,2))
#print(BA,BA.shape)

BAfile = open(os.path.join(save_path, "BA_%s_%s_step%d.txt"%(pathindex[0],pathindex[1],mda_step)), 'w+')
BAfile.write("#%12s%12s%12s\n"%("Index","B_N0O1", "A_C1O1H0"))
Index=[]
for k in range(BA.shape[0]):
    Index.append(k+1)
    BAfile.write(' %12d%12.4f%12.4f\n' %(k+1,BA[k][0],BA[k][1]))
BAfile.close()

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
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)   # fontsize of the figure title

def BA_Fig(data):
    fig = plt.figure(figsize=(7,6), dpi=150, facecolor='white')
    ax = fig.add_subplot(1, 1, 1)
    x = data[:,1]*180/np.pi  # first column
    y = data[:,0]  # second column
    h = ax.hist2d(x, y, bins=20, cmap='Blues')
    cbar = fig.colorbar(h[3], ax=ax)  
    #ax.legend()
    ax.set_ylabel("N-O Bond")
    ax.set_xlabel("C-O-H Angle")
    #ax.set_xlim(0,180)
    fig.savefig(os.path.join(save_path, "BA.png"), dpi=600, bbox_inches='tight')

BA_Fig(BA)