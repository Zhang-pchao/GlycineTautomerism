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

def find_min_bond_index(u,n_idx,h_idx):
    # select hydrogen and oxygen atoms
    h_atoms = u.select_atoms(h_idx)
    o_atoms = u.select_atoms(n_idx)

    # initialize minimum distance and closest oxygen index
    min_distance = np.inf
    closest_o_index = -1
    
    # loop over hydrogen and oxygen atoms
    d_list,h_list=[.0,.0,.0],[0,0,0]
    for h_atom in h_atoms:
        #print(h_atom.index)
        for j, o_atom in enumerate(o_atoms):
            # calculate distance between hydrogen and oxygen atoms
            distance = distances.calc_bonds(h_atom.position, o_atom.position,box=u.dimensions)
            
            # update closest oxygen index and minimum distance if current oxygen is closer
            if distance < min_distance:
                #closest_o_index = j
                closest_h_index = h_atom.index
                min_distance = distance
                d_list[2]=d_list[1]
                d_list[1]=d_list[0]
                d_list[0]=min_distance
                h_list[2]=h_list[1]
                h_list[1]=h_list[0]
                h_list[0]=closest_h_index                
    
    return d_list,h_list

def find_min_bond_index_old(u,n_idx,o_idx1,o_idx2):
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

trjs = []
for filename in os.listdir(trj_path):
    if filename.endswith('.lammpstrj'):
        trjs.append(filename)
#print(trjs)

trj_file    = os.path.join(trj_path,trjs[0])
u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')
u.trajectory[1]
c1_idx,c2_idx,d_nc1,_idx = find_min_bond_index_old(u,"index %d"%(int(NOOdict['N'])-1),NOOdict['C1'],NOOdict['C2'])

#for k in range(len(trjs)):
for k in [0]:    
    trj_file    = os.path.join(trj_path,trjs[k])
    u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')
    len_trj  = len(u.trajectory)
    #trj_info = [len_trj,trj_skip,mda_step]
    BA,AD = np.array([]),np.array([])
    for j in range(len_trj):       
        u.trajectory[j]
        d_list,h_list = find_min_bond_index(u,"index %d"%(int(NOOdict['N'])-1),"type 1") 
        if d_list[0] <= oh_cutoff and d_list[1] <= oh_cutoff:
            #print(d_list[0],d_list[1],d_list[2])
            if d_list[2] > oh_cutoff:#NH2
                #o1_idx,o2_idx,d_no1,_idx = find_min_bond_index(u,"index %d"%(int(NOOdict['N'])-1),NOOdict['O1'],NOOdict['O2'])               
                H1 = u.select_atoms("index %d"%h_list[0])[0].position
                H2 = u.select_atoms("index %d"%h_list[1])[0].position
                H3 = u.select_atoms("index %d"%h_list[2])[0].position
                N0 = u.select_atoms("index %d"%(int(NOOdict['N'])-1))[0].position
                C1 = u.select_atoms("index %d"%c1_idx)[0].position
                C2 = u.select_atoms("index %d"%c2_idx)[0].position
                #print(O1,O2,N0,C1,C2)
                B_C2H1     = distances.calc_bonds(C2, H1, box=u.dimensions)
                B_C2H2     = distances.calc_bonds(C2, H2, box=u.dimensions)
                A_C2N0H1   = distances.calc_angles(C2, N0, H1, box=u.dimensions)
                A_C2N0H2   = distances.calc_angles(C2, N0, H2, box=u.dimensions)
                #D_N0C1C2O1 = abs(distances.calc_dihedrals(N0, C1, C2, O1, box=u.dimensions))
                ba = np.array([np.average([B_C2H1,B_C2H2]),
                               np.average([A_C2N0H1,A_C2N0H2])])#/np.pi*180
                BA = np.concatenate((BA,ba),axis=0)
            else:#NH3
                H1 = u.select_atoms("index %d"%h_list[0])[0].position
                H2 = u.select_atoms("index %d"%h_list[1])[0].position
                H3 = u.select_atoms("index %d"%h_list[2])[0].position
                N0 = u.select_atoms("index %d"%(int(NOOdict['N'])-1))[0].position
                C1 = u.select_atoms("index %d"%c1_idx)[0].position
                C2 = u.select_atoms("index %d"%c2_idx)[0].position
                A_H1N0H2 = distances.calc_angles(H1, N0, H2, box=u.dimensions)
                A_H2N0H3 = distances.calc_angles(H2, N0, H3, box=u.dimensions)
                A_H3N0H1 = distances.calc_angles(H3, N0, H1, box=u.dimensions)
                D_C2C1N0H1 = abs(distances.calc_dihedrals(C2, C1, N0, H1, box=u.dimensions))
                D_C2C1N0H2 = abs(distances.calc_dihedrals(C2, C1, N0, H2, box=u.dimensions))
                D_C2C1N0H3 = abs(distances.calc_dihedrals(C2, C1, N0, H3, box=u.dimensions))
                ad = np.array([np.average([A_H1N0H2,A_H2N0H3,A_H3N0H1]),min(D_C2C1N0H1,D_C2C1N0H2,D_C2C1N0H3)])
                AD = np.concatenate((AD,ad),axis=0)
               
BA = BA.reshape((-1,2))
#print(BA,BA.shape)
AD = AD.reshape((-1,2))
#print(BAD,BAD.shape)

BAfile = open(os.path.join(save_path, "BA_CH_CNH.txt"), 'w+')
BAfile.write("#%12s%12s%12s\n"%("Index","B_C2H", "A_C2N0H"))
Index=[]
for k in range(BA.shape[0]):
    Index.append(k+1)
    BAfile.write(' %12d%12.4f%12.4f\n' %(k+1,BA[k][0],BA[k][1]))
BAfile.close()

ADfile = open(os.path.join(save_path, "AD_HNH_CCNH.txt"), 'w+')
ADfile.write("#%12s%12s%12s\n"%("Index","A_HN0H", "D_C2C1N0H"))
Index=[]
for k in range(BAD.shape[0]):
    Index.append(k+1)
    ADfile.write(' %12d%12.4f%12.4f\n' %(k+1,AD[k][0],AD[k][1]))
ADfile.close()

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

def BAD_Fig(data,name):
    fig = plt.figure(figsize=(7.4,6), dpi=150, facecolor='white')
    ax = fig.add_subplot(1, 1, 1)
    #ax.legend()
    #ax.set_xlim(0,180)
    if name == "BA":
        x = data[:,1]*180/np.pi  # first column
        y = data[:,0]  # second column
        h = ax.hist2d(x, y, bins=20, cmap='Blues',density=True)
        cbar = fig.colorbar(h[3], ax=ax)  
        ax.set_ylabel("C-H Bond")
        ax.set_xlabel("C-N-H Angle")
        fig.savefig(os.path.join(save_path, "BA_CH_CNH.png"), dpi=600, bbox_inches='tight')
    if name == "AD":
        x = data[:,1]*180/np.pi  # first column
        y = data[:,0]*180/np.pi  # second column
        h = ax.hist2d(x, y, bins=20, cmap='Blues',density=True)
        cbar = fig.colorbar(h[3], ax=ax)          
        ax.set_ylabel("H-N-H Angle")
        ax.set_xlabel("C-C-N-H Dihedral")
        fig.savefig(os.path.join(save_path, "AD_HNH_CCNH.png"), dpi=600, bbox_inches='tight')

BAD_Fig(BA,"BA")
BAD_Fig(AD,"AD")
