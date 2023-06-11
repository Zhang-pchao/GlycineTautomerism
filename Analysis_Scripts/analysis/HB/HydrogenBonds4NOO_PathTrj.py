#!/usr/bin/env python
# coding: utf-8

import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import os
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

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

def get_avg_std(x):
    avg = np.average(x, axis=0)
    std = np.std(x, axis=0, ddof=1)
    return avg,std

####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/glycine/dpmd/geofile"
trj_path    = "../"
data_geo    = "CanonicalGlycine54H2O_opt_atomic.data"
save_path   = "./"
trj_skip    = 1
mda_step    = 1
NOOdict     = {'N': 163,'O1': 169,'O2': 171} #serial in lammps .data file
pathindex   = ['anionic_high','anionic_low'] #'cationic','anionic'
####################change above####################
data_geo    = os.path.join(geo_path,data_geo)
trjs = []
for filename in os.listdir(trj_path):
    if filename.endswith('.lammpstrj'):
        trjs.append(filename)
print(trjs)

HBfile = open(os.path.join(save_path, "HB_%s_%s_step%d.txt"%(pathindex[0],pathindex[1],mda_step)), 'w+')
HBfile.write("#%12s%12s%12s%12s%12s%12s%12s%12s%12s\n"%("Index",
                                        "N_d_avg","N_d_std","O_d_avg","O_d_std",
                                        "N_a_avg","N_a_std","O_a_avg","O_a_std"))
                                        
D__N_avg,D__N_std,D__O_avg,D__O_std = [],[],[],[]
A__N_avg,A__N_std,A__O_avg,A__O_std = [],[],[],[]
Index = []
for k in range(len(trjs)):
    trj_file    = os.path.join(trj_path,trjs[k])
    u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')
    len_trj  = len(u.trajectory)
    trj_info = [len_trj,trj_skip,mda_step]
    
    # NOO as donors
    d__N = get_hb_num("N", trj_info,donor=True)
    d__N_avg,d__N_std = get_avg_std(d__N)
    d_O1 = get_hb_num("O1",trj_info,donor=True)
    d_O2 = get_hb_num("O2",trj_info,donor=True)
    d__O = np.concatenate((d_O1, d_O2), axis=0)
    d__O_avg,d__O_std = get_avg_std(d__O)
    # NOO as accepters
    a__N = get_hb_num("N", trj_info,donor=False)
    a__N_avg,a__N_std = get_avg_std(a__N)
    a_O1 = get_hb_num("O1",trj_info,donor=False)
    a_O2 = get_hb_num("O2",trj_info,donor=False)
    a__O = np.concatenate((a_O1, a_O2), axis=0)
    a__O_avg,a__O_std = get_avg_std(a__O) 
    
    D__N_avg.append(d__N_avg)
    D__N_std.append(d__N_std)
    D__O_avg.append(d__O_avg)
    D__O_std.append(d__O_std)
    A__N_avg.append(a__N_avg)
    A__N_std.append(a__N_std)
    A__O_avg.append(a__O_avg)
    A__O_std.append(a__O_std)
    Index.append(k+1)
    
    HBfile.write(' %12d%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f\n' %(k,
                    d__N_avg,d__N_std,d__O_avg,d__O_std,
                    a__N_avg,a__N_std,a__O_avg,a__O_std))
HBfile.close()

colors = [(31 ,59 ,115), 
#          (47 ,146,148), 
#          (80 ,178,141),
#          (167,214,85 ),
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
    ax.axes.xaxis.set_ticklabels([])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.set_xticks([x[0],(x[0]+x[-1])/2,x[-1]],[pathindex[0],'transitional',pathindex[1]])
    ax.set_ylim((-0.5, 3.5))
    ax.set_ylabel("No. of %ss"%name[-1])
    fig.savefig(os.path.join(save_path, "HB_%s_%s_%s.png"%(pathindex[0],pathindex[1],name[-1])), dpi=600, bbox_inches='tight')

HBnumber_Fig(Index,[np.asarray(D__N_avg),np.asarray(D__O_avg)],[np.asarray(D__N_std),np.asarray(D__O_std)],['N','O','donor'],pathindex)
HBnumber_Fig(Index,[np.asarray(A__N_avg),np.asarray(A__O_avg)],[np.asarray(A__N_std),np.asarray(A__O_std)],['N','O','acceptor'],pathindex)
