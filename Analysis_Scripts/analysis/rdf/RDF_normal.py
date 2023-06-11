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
from MDAnalysis.analysis.rdf import InterRDF,InterRDF_s

colors = [(31 ,59 ,115), 
          (47 ,146,148), 
          (80 ,178,141),
#          (167,214,85 ),
          (255,224,62 ),
          (255,169,85 ),
          (189, 48,47 )]
colors=list(tuple(i/255 for i in color) for color in colors)
colors.reverse()

#set bigger font sizes
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIG_SIZE = 17
plt.rc('font', size=SMALL_SIZE)        # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)   # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE-1)  # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)   # fontsize of the figure title

####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/glycine/dpmd/geofile"
trj_path    = "../"
save_path   = "./"
data_geo    = "127w_1h3o.data"
trj_name    = "glycine_10.lammpstrj"
trj_skip    = 1
mda_step    = 10
type_map    = {'H': 1,'O': 2,'N': 3,'C': 4}
#NOOdict     = {'N': 163,'O1': 169,'O2': 171} #serial in lammps .data file
####################change above####################

data_geo    = os.path.join(geo_path,data_geo)
trj_file    = os.path.join(trj_path,trj_name)
u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP',dt=1e-3) #dt ps

def get_rdf(c1,c2,u,type_map,trj_skip,mda_step,save_path):
    s1 = u.select_atoms('type %s'%type_map[c1]) # H
    s2 = u.select_atoms('type %s'%type_map[c2]) # O
    half_box = u.dimensions[0]/2
    rdf = InterRDF(s1,s2,nbins=int(100*half_box), range=(0.01, half_box))
    rdf.run(start=trj_skip, step=mda_step) # stop=1000
    y_rdf = rdf.rdf
    #y_cdf = rdf.get_cdf()
    x_rdf = rdf.bins
    
    file = open(os.path.join(save_path, "RDF_{0}_{1}.txt".format(c1,c2)), 'w+')
    file.write("#RDF %4s-%4s\n"%(c1,c2))
    file.write('#%8s%8s\n' %('bin','rdf'))
    for i in range(len(x_rdf)):
        file.write(' %8.4f%8.4f\n' %(x_rdf[i],y_rdf[i]))
    file.close()
    return x_rdf,y_rdf

def get_fig(c1,c2,x_rdf,y_rdf,colors,save_path):
    fig = plt.figure(figsize=(8,3), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    ax.grid(True)
    ax.plot(x_rdf, y_rdf, alpha=0.8, lw=3, color=colors[0], label='%2s –%2s'%(c1,c2))
    ax.legend(loc="upper right")
    ax.set_xlim(0,6)
    ax.set_ylim(0,6)
    ax.set_xlabel("r "+r"$\ \rm (\AA)$")
    ax.set_ylabel("g(r)")
    fig.show()
    fig.savefig(os.path.join(save_path, "RDF_{0}_{1}.png".format(c1,c2)), dpi=600)

x_rdf,y_rdf = get_rdf('H','O',u,type_map,trj_skip,mda_step,save_path)
get_fig('H','O',x_rdf,y_rdf,colors,save_path)
x_rdf,y_rdf = get_rdf('H','H',u,type_map,trj_skip,mda_step,save_path)
get_fig('H','H',x_rdf,y_rdf,colors,save_path)
x_rdf,y_rdf = get_rdf('O','O',u,type_map,trj_skip,mda_step,save_path)
get_fig('O','O',x_rdf,y_rdf,colors,save_path)
