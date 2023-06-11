#!/usr/bin/env python
# coding: utf-8

import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import os
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.lib import distances

### Parser stuff ###
parser = argparse.ArgumentParser(description='get the RDF of glycine in water')
# files
parser.add_argument('--trjdir','-td',dest='tdirectory',type=str,default='./',help='the directory of trajectory files')
parser.add_argument('--savdir','-sd',dest='sdirectory',type=str,default='./',help='the directory of saving files')
parser.add_argument('--fstart',dest='fs',type=str,default='R',help='the name of starting frame')
parser.add_argument('--festop',dest='fe',type=str,default='P',help='the name of ending   frame')
# some easy parsing
args=parser.parse_args()

def set_hbonds(d_idx,a_idx,bond=3.5,angle=140):
    hbonds = HydrogenBondAnalysis(
        universe           = u,
        donors_sel         = d_idx, # O="type 2"
        hydrogens_sel      = "type 1", # H
        acceptors_sel      = a_idx,
        d_a_cutoff         = bond,      # <3.5
        d_h_a_angle_cutoff = angle,      # >140
        update_selections  = False
    )
    return hbonds

def get_hb_num(i,trj_info,donor=True):
    if donor:
        hbonds = set_hbonds("index %d"%(NOOdict[i]-1),"type 2") # serial-1=index in mda
    else:#accepter
        hbonds = set_hbonds("type 2","index %d"%(NOOdict[i]-1))
    hbonds.run(
        start   = trj_info[1],                   # skip 10% trajectorys, default:None
        stop    = None,
        step    = trj_info[2],
        verbose = True
    )
    counts_hbs = hbonds.count_by_time()
    
    return counts_hbs

def get_hb_num_4OO(Oidx,trj_info,donor=True):
    jj=0
    for j in range(trj_info[1],trj_info[0],trj_info[2]):
        #print(j)
        if donor:
            hbonds = set_hbonds("index %d"%Oidx[j],"type 2")
        else:#accepter
            hbonds = set_hbonds("type 2","index %d"%Oidx[j])
        
        hbonds.run(
            start   = j,                  
            stop    = j+1,
            step    = trj_info[2],
            verbose = False
        )
        if jj == 0:
            counts_hbs = hbonds.count_by_time()
        else:
            counts_hbs_tmp = hbonds.count_by_time()
            counts_hbs = np.concatenate((counts_hbs, counts_hbs_tmp), axis=0)
        jj+=1
        
    return counts_hbs

def find_Oh_O_index(u,o_idx1,o_idx2):
    # select hydrogen and oxygen atoms
    h_atoms = u.select_atoms("type 1")
    o_atoms = u.select_atoms("index %d %d"%(o_idx1-1,o_idx2-1))

    # initialize minimum distance and closest oxygen index
    min_distance = np.inf
    closest_o_index = -1
    
    # loop over hydrogen and oxygen atoms
    for h_atom in h_atoms:
        for j, o_atom in enumerate(o_atoms):
            # calculate distance between hydrogen and oxygen atoms
            #distance = np.linalg.norm(h_atom.position - o_atom.position)
            distance = distances.calc_bonds(h_atom.position, o_atom.position,box=u.dimensions)
            
            # update closest oxygen index and minimum distance if current oxygen is closer
            if distance < min_distance:
                closest_o_index = j
                min_distance = distance
    
    # print the index of the oxygen atom with the shortest O-H bond distance
    print(f"Oxygen atom index : {o_atoms[closest_o_index].index}, {min_distance}")
    Oh_index = o_atoms[closest_o_index].index
    O_index = int(o_idx1-1+o_idx2-1-Oh_index)
    return Oh_index,O_index

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
NOOdict     = {'N': 385,'O1': 391,'O2': 393} #serial in lammps .data file
pathindex   = [args.fs,args.fe] #'cationic','anionic'
####################change above####################
data_geo    = os.path.join(geo_path,data_geo)

trjs = []
for filename in os.listdir(trj_path):
    if filename.endswith('.lammpstrj'):
        trjs.append(filename)
#print(trjs)

HBfile = open(os.path.join(save_path, "HB_%s_%s_step%d.txt"%(pathindex[0],pathindex[1],mda_step)), 'w+')
HBfile.write("#%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n"%("Index",
                                        "N_d_avg","N_d_std","Oh_d_avg","Oh_d_std","O_d_avg","O_d_std",
                                        "N_a_avg","N_a_std","Oh_a_avg","Oh_a_std","O_a_avg","O_a_std"))

D__N_avg,D__N_std,D_O1_avg,D_O1_std,D_O2_avg,D_O2_std = [],[],[],[],[],[]
A__N_avg,A__N_std,A_O1_avg,A_O1_std,A_O2_avg,A_O2_std = [],[],[],[],[],[]
Index = []
for k in range(len(trjs)):
#for k in [3]:    
    trj_file    = os.path.join(trj_path,trjs[k])
    u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')
    len_trj  = len(u.trajectory)
    trj_info = [len_trj,trj_skip,mda_step]
    
    Oh_index_list,O_index_list=[],[]
    for i in range(trj_info[1],trj_info[0],trj_info[2]):
        u.trajectory[i]
        Oh_index,O_index = find_Oh_O_index(u,NOOdict['O1'],NOOdict['O2'])
        Oh_index_list.append(Oh_index)
        O_index_list.append(O_index)
    print(Oh_index_list,O_index_list)    
    # NOO as donors
    d__N = get_hb_num("N", trj_info,donor=True)
    d__N_avg,d__N_std = get_avg_std(d__N)
    d_O1 = get_hb_num_4OO(Oh_index_list,trj_info,donor=True)
    d_O2 = get_hb_num_4OO(O_index_list,trj_info,donor=True)
    d_O1_avg,d_O1_std = get_avg_std(d_O1)
    d_O2_avg,d_O2_std = get_avg_std(d_O2)
    # NOO as accepters
    a__N = get_hb_num("N", trj_info,donor=False)
    a__N_avg,a__N_std = get_avg_std(a__N)
    a_O1 = get_hb_num_4OO(Oh_index_list,trj_info,donor=False)
    a_O2 = get_hb_num_4OO(O_index_list,trj_info,donor=False)
    a_O1_avg,a_O1_std = get_avg_std(a_O1) 
    a_O2_avg,a_O2_std = get_avg_std(a_O2) 
    D__N_avg.append(d__N_avg)
    D__N_std.append(d__N_std)
    D_O1_avg.append(d_O1_avg)
    D_O1_std.append(d_O1_std)
    D_O2_avg.append(d_O2_avg)
    D_O2_std.append(d_O2_std)    
    A__N_avg.append(a__N_avg)
    A__N_std.append(a__N_std)
    A_O1_avg.append(a_O1_avg)
    A_O1_std.append(a_O1_std)
    A_O2_avg.append(a_O2_avg)
    A_O2_std.append(a_O2_std)    
    Index.append(k+1)
    
    HBfile.write(' %12d%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f\n' %(k,
                    d__N_avg,d__N_std,d_O1_avg,d_O1_std,d_O2_avg,d_O2_std,
                    a__N_avg,a__N_std,a_O1_avg,a_O1_std,a_O2_avg,a_O2_std))
HBfile.close()

colors = [(31 ,59 ,115), 
#          (47 ,146,148), 
#          (80 ,178,141),
          (167,214,85 ),
#          (255,224,62 ),
#          (255,169,85 ),
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

def HBnumber_Fig(x,y,s,name,pathidex):
    fig = plt.figure(figsize=(8,3), dpi=150, facecolor='white')
    ax = fig.add_subplot(1, 1, 1)
    for i in range(len(y)):
        ax.plot(x, y[i], lw=3, c=colors[i],label="%s of glycine"%name[i])
        ax.fill_between(x,y[i]-s[i],y[i]+s[i],facecolor=colors[i],alpha = 0.6)
    ax.legend()
    #ax.axes.xaxis.set_ticklabels([])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    #ax.set_xticks([x[0],(x[0]+x[-1])/2,x[-1]],[pathindex[0],'transitional',pathindex[1]])
    ax.set_ylim((-0.5, 3.5))
    ax.set_ylabel("No. of %ss"%name[-1])
    fig.savefig(os.path.join(save_path, "HB_%s_%s_%s_step%d.png"%(pathindex[0],pathindex[1],name[3],mda_step)), dpi=600, bbox_inches='tight')

HBnumber_Fig(Index,[np.asarray(D__N_avg),np.asarray(D_O1_avg),np.asarray(D_O2_avg)],
                   [np.asarray(D__N_std),np.asarray(D_O1_std),np.asarray(D_O2_std)],
             ['N','Oh','O','donor'],pathindex)
HBnumber_Fig(Index,[np.asarray(A__N_avg),np.asarray(A_O1_avg),np.asarray(A_O2_avg)],
                   [np.asarray(A__N_std),np.asarray(A_O1_std),np.asarray(A_O2_std)],
             ['N','Oh','O','acceptor'],pathindex)
