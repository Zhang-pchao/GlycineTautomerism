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

####################change below####################
geo_path    = "../"
trj_path    = "../"
save_path   = "./"
data_geo    = "1000000.data"
trj_name    = "glycine_1b.lammpstrj"
trj_skip    = 1
mda_step    = 1
NOOdict     = {'N': 3058,'O1': 3054,'O2': 3053} #serial in lammps .data file
Ion_flag    = "Zwitterionic"
####################change above####################

HBfile = open(os.path.join(save_path, "HBnumber_step{0}_{1}.txt".format(mda_step,trj_name.split('.')[0])), 'w+')
data_geo    = os.path.join(geo_path,data_geo)
trj_file    = os.path.join(trj_path,trj_name)
u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP')
#print(u.atoms[162].index)
len_trj     = len(u.trajectory)
trj_info = [len_trj,trj_skip,mda_step]

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

#trj_info = [len_trj,trj_skip,mda_step]
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
    hbs_num = hbonds.results.hbonds.shape[0]
    trj_num = (trj_info[0]-trj_info[1])//trj_info[2]
    hbs_avg = hbs_num/trj_num
    
    return hbs_avg

# NOO as donors
d__N = get_hb_num("N", trj_info,donor=True)
d_O1 = get_hb_num("O1",trj_info,donor=True)
d_O2 = get_hb_num("O2",trj_info,donor=True)
# NOO as accepters
a__N = get_hb_num("N", trj_info,donor=False)
a_O1 = get_hb_num("O1",trj_info,donor=False)
a_O2 = get_hb_num("O2",trj_info,donor=False)

HBfile.write("No. of Donors\n")
HBfile.write('# %8s%8s%8s\n' %('N','O1','O2'))
HBfile.write('  %8.2f%8.2f%8.2f\n\n' %(d__N,d_O1,d_O2))
HBfile.write("No. of Acceptors\n")
HBfile.write('# %8s%8s%8s\n' %('N','O1','O2'))
HBfile.write('  %8.2f%8.2f%8.2f\n\n' %(a__N,a_O1,a_O2))
HBfile.close()

# [frame, donor_index, hydrogen_index, acceptor_index, DA_distance, DHA_angle]
#print(hbonds.results.hbonds.shape)
#print(hbonds.results.hbonds[-1])

#for frame, donor_ix, *_ in hbonds.results.hbonds:
#    u.trajectory[frame.astype(int)]
#    donor = u.atoms[donor_ix.astype(int)]
#    zpos = donor.position[2]
#    hist, *_ = np.histogram(zpos, bins=bin_edges)
#    # multiply by two as each hydrogen bond involves two water molecules  
#    idx = int((frame-trj_skip)/mda_step)
#    counts_list[idx] += hist * 2

# HBnumber_Fig
fig = plt.figure(figsize=(8,4), dpi=150, facecolor='white')
ax = fig.add_subplot(1, 1, 1)
ax.bar([0.8,2.8,4.8], [d__N,d_O1,d_O2], width=0.3,alpha=0.9,label='Donor',    color='orange')
ax.bar([1.2,3.2,5.2], [a__N,a_O1,a_O2], width=0.3,alpha=0.9,label='Acceptor', color='royalblue')
ax.legend(fontsize = 14)
ax.axes.xaxis.set_ticklabels([])
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.text(0.95,-0.2,'N', color = 'black',fontsize=14)
ax.text(2.9, -0.2,'O1',color = 'black',fontsize=14)
ax.text(4.9, -0.2,'O2',color = 'black',fontsize=14)
ax.set_title("One %s glycine in water"%Ion_flag)
ax.tick_params(labelsize=12)
ax.set_ylim((0, 3))
ax.set_ylabel("No. of Donors or Acceptors",fontsize = 14)
#fig.show()
fig.savefig(os.path.join(save_path, "HBnumber_step{0}_{1}.png".format(mda_step,trj_name.split('.')[0])), dpi=600, bbox_inches='tight')
