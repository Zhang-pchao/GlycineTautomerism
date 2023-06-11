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
#          (47 ,146,148), 
          (80 ,178,141),
#          (167,214,85 ),
          (255,224,62 ),
#          (255,169,85 ),
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
geo_path    = "/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/018_water/004_127H2O_H3O/FindIonIdex"
trj_path    = "/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/018_water/004_127H2O_H3O/FindIonIdex"
save_path   = trj_path
data_geo    = "addIon.data"
trj_name    = "add_Ion_O_type_9.lammpstrj"
trj_skip    = 1
mda_step    = 100
type_map    = {'H': 1,'O': 2,'N': 3,'C': 4,'IonO':9}
Ion_flag    = 'h3o'
####################change above####################
waterbox    = {'64': 6.23,'128': 7.84,'256': 9.88,'512': 12.45,'1024':15.69}
data_geo    = os.path.join(geo_path,data_geo)
trj_file    = os.path.join(trj_path,trj_name)
u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP') #dt ps,dt=1e-3

def get_rdf(c1,c2,u,type_map,trj_skip,mda_step,save_path):
    s1 = u.select_atoms('type %s'%type_map[c1[0]]) # H
    s2 = u.select_atoms('type %s'%type_map[c2]) # O
    half_box = u.dimensions[0]/2
    rdf = InterRDF(s1,s2,nbins=int(100*half_box), range=(0.01, half_box))
    rdf.run(start=trj_skip, step=mda_step) # stop=1000
    y_rdf = rdf.rdf
    #y_cdf = rdf.get_cdf()
    x_rdf = rdf.bins
    
    file = open(os.path.join(save_path, "RDF_{0}*_{1}.txt".format(c1[1],c2)), 'w+')
    file.write("#RDF %4s*-%s\n"%(c1[1],c2))
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
    ax.set_xlim(0,16)
    ax.set_ylim(0,3)
    ax.set_xlabel("r "+r"$\ \rm (\AA)$")
    ax.set_ylabel("g(r)")
    #fig.show()
    fig.savefig(os.path.join(save_path, "RDF_{0}_{1}.png".format(c1,c2)), dpi=600, bbox_inches='tight')

def get_rdf4Ion(c1,c2,u,type_map,trj_skip,mda_step,save_path,cutoff=1.3):
    s1 = u.select_atoms('(around %f type %s) and (type %s)'%(cutoff,type_map[c1[0]],type_map[c1[1]])) # H* of Ion
    s2 = u.select_atoms('type %s'%type_map[c2])
    half_box = u.dimensions[0]/2
    rdf = InterRDF(s1,s2,nbins=int(100*half_box), range=(0.01, half_box))
    rdf.run(start=trj_skip, step=mda_step) # stop=1000
    y_rdf = rdf.rdf
    #y_cdf = rdf.get_cdf()
    x_rdf = rdf.bins
    
    file = open(os.path.join(save_path, "RDF_{0}*_{1}.txt".format(c1[1],c2)), 'w+')
    file.write("#RDF %4s*-%s\n"%(c1[1],c2))
    file.write('#%8s%8s\n' %('bin','rdf'))
    for i in range(len(x_rdf)):
        file.write(' %8.4f%8.4f\n' %(x_rdf[i],y_rdf[i]))
    file.close()
    return x_rdf,y_rdf

def get_allfig(c1,c2,x_rdf,y_rdf,waterbox,colors,Ion_flag,save_path):
    fig = plt.figure(figsize=(8,3), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    ax.grid(True)
    label0='$\mathregular{H_2}$'+'O'
    label1='$\mathregular{H_3}$'+'$\mathregular{O^+}$'
    label2='$\mathregular{OH^-}$'
    for i in range(len(x_rdf)):
        ax.plot(x_rdf[i], y_rdf[i], alpha=0.8, lw=3, color=colors[i], label='%2s –%2s'%(c1[i],c2[i]))
    ax.legend(loc="upper right")
    for i in waterbox.keys():
        ax.vlines(waterbox[i],0,3, colors="black",lw=1,ls='--')
        ax.text(waterbox[i]-0.5,0.075,'%4s %s'%(i,label0), rotation=90, color = 'black',fontsize=9)
    ax.set_xlim(0,16)
    ax.set_ylim(0,3)
    if Ion_flag == 'h3o':
        ax.set_title("One %s in water"%label1)
    else: #'oh'
        ax.set_title("One %s in water"%label2)
    ax.set_xlabel("r "+r"$\ \rm (\AA)$")
    ax.set_ylabel("g(r)")
    #fig.show()
    fig.savefig(os.path.join(save_path, "RDF_all.png"), dpi=600, bbox_inches='tight')

x_rdf1,y_rdf1 = get_rdf(['IonO','O'],'O',u,type_map,trj_skip,mda_step,save_path)
get_fig('O*','O',x_rdf1,y_rdf1,colors,save_path)
x_rdf2,y_rdf2 = get_rdf(['IonO','O'],'H',u,type_map,trj_skip,mda_step,save_path)
get_fig('O*','H',x_rdf2,y_rdf2,colors,save_path)
x_rdf3,y_rdf3 = get_rdf4Ion(['IonO','H'],'O',u,type_map,trj_skip,mda_step,save_path)
get_fig('H*','O',x_rdf3,y_rdf3,colors,save_path)
x_rdf4,y_rdf4 = get_rdf4Ion(['IonO','H'],'H',u,type_map,trj_skip,mda_step,save_path)
get_fig('H*','H',x_rdf4,y_rdf4,colors,save_path)

get_allfig(['O*','O*','H*','H*'],['O','H','O','H'],
           [x_rdf1,x_rdf2,x_rdf3,x_rdf4],
           [y_rdf1,y_rdf2,y_rdf3,y_rdf4],
           waterbox,colors,Ion_flag,save_path)