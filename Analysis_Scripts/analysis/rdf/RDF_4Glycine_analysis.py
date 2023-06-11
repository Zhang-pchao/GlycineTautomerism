#!/usr/bin/env python
# coding: utf-8

import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import os
import re
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.rdf import InterRDF,InterRDF_s

colors = [(31 ,59 ,115), 
#          (47 ,146,148), 
#          (80 ,178,141),
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
rdfpath1    = "/rdf/zwitterion/1016_0.8ns"
rdfpath2    = "/rdf/canonical/1016_0.8ns"
save_path   = "/rdf/rdf_analysis/zwi_can"
txt_name    = ["RDF_N*_O.txt","RDF_N*_H.txt","RDF_O*_O.txt","RDF_O*_H.txt"]
waterbox    = {'64': 6.23,'128': 7.84,'256': 9.88,'512': 12.45,'1024':15.69}
all_flag    = ['Zwitterionic','Canonical','Anionic','Cationic']
Ion_flag    = [all_flag[0],all_flag[1]]

def get_allfig(x,y,ele,waterbox,colors,Ion_flag,save_path):
    fig = plt.figure(figsize=(8,3), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    ax.grid(True)
    label0='$\mathregular{H_2}$'+'O'
    label1='$\mathregular{H_3}$'+'$\mathregular{O^+}$'
    label2='$\mathregular{OH^-}$'
    for i in range(len(x)):
        ax.plot(x[i], y[i], alpha=0.8, lw=2, color=colors[i], 
                label='%12s [%2s –%2s]'%(Ion_flag[i],ele[0],ele[1]))
    deltax = x[0]-x[1]
    if deltax[0]==0.0 and deltax[-1]==0.0:
        ax.plot(x[0], y[0]-y[1], alpha=0.8, lw=2, color=colors[-1], 
                label='%15s [%2s –%2s]'%('Delta',ele[0],ele[1]))
    ax.legend(loc="upper right")
    for i in waterbox.keys():
        ax.vlines(waterbox[i],-1,3, colors="black",lw=1,ls='--')
        ax.text(waterbox[i]-0.5,-0.9,'%4s'%(i), rotation=90, color = 'black',fontsize=9)
    ax.set_xlim(0,16)
    ax.set_ylim(-1,3)
    #ax.set_title("One %s glycine in water"%Ion_flag)
    ax.set_xlabel("r "+r"$\ \rm (\AA)$")
    ax.set_ylabel("g(r)")
    fig.savefig(os.path.join(save_path, "RDF_delta_%s_%s_%s_%s.png"%(ele[0],ele[1],Ion_flag[0],Ion_flag[1])), 
                dpi=600, bbox_inches='tight')

for i in txt_name:
    rdftxt1     = os.path.join(rdfpath1,i)
    rdftxt2     = os.path.join(rdfpath2,i)
    rdf_x1      = np.loadtxt(rdftxt1,usecols=0)
    rdf_y1      = np.loadtxt(rdftxt1,usecols=1)
    rdf_x2      = np.loadtxt(rdftxt2,usecols=0)
    rdf_y2      = np.loadtxt(rdftxt2,usecols=1)
    get_allfig([rdf_x1,rdf_x2],[rdf_y1,rdf_y2],
           [i.split('_')[1],i.split('_')[2].split('.')[0]],
           waterbox,colors,Ion_flag,save_path)