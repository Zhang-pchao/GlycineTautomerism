#!/usr/bin/env python
# coding: utf-8

import plumed
import matplotlib.pyplot as plt
import os
import MDAnalysis
import numpy as np
import math
import sys
import argparse
#reuse same plumed kernel, to avoid multiple warnings
PLUMED_KERNEL=plumed.Plumed()
from functools import partial
plumed.read_as_pandas = partial(plumed.read_as_pandas, kernel=PLUMED_KERNEL)

#plot-related stuff
import matplotlib.pyplot as plt
from matplotlib import cm
from IPython.display import clear_output

### Parser stuff ###
parser = argparse.ArgumentParser(description='get the FES estimate used by OPES')
# files
parser.add_argument('--dir','-d',dest='directory',type=str,default='./',help='the directory of all files')
# some easy parsing
args=parser.parse_args()

#set bigger font sizes
SMALL_SIZE = 11
MEDIUM_SIZE = 12
BIG_SIZE = 15
plt.rc('font', size=SMALL_SIZE)        # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)   # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)   # fontsize of the figure title

def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)

def cmap2colorlist(cmap_name, color_numb):
    colormap = cm.get_cmap(cmap_name, color_numb+4)
    idx = np.arange(color_numb) + 3
    colorlist = colormap(idx)
    return colorlist

def fes_1D(data,detf,detx,xlist,ylist,err,ylabel,Xlabel,Ylabel,path='./',step=1,Legend=True):
    fig = plt.figure(figsize=(6,6), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    colors = cmap2colorlist('GnBu', len(data))
    for i in range(len(data)):
        ax.plot(data[i][xlist[0]][::step],data[i][ylist[0]][::step],
                   color=colors[i],label="%s"%(ylabel[i]))
                   #color=colors[i],label="%s Delta free energy = %5.2f kJ/mol"%(ylabel[i],detf[i]))
        ax.fill_between(data[i][xlist[0]][::step],
                        data[i][ylist[0]][::step]-data[i][err[0]][::step],
                        data[i][ylist[0]][::step]+data[i][err[0]][::step],
                        facecolor = colors[i],alpha = 0.3)
        #ax.vlines(detx[i],min(data[i][ylist[0]][::step]),max(data[i][ylist[0]][::step]),colors="black",lw=1,ls='--')
    ax.set_xlabel(Xlabel)
    ax.set_ylabel(Ylabel)
    if Legend:
        ax.legend()
    fig.savefig(os.path.join(path,"%s_%s.png"%(Xlabel,ylabel)), dpi=150, bbox_inches='tight')

def read_deltaf(file):
    with open(file) as f: 
        lines = f.readlines()
        for line in lines:
            if "DeltaF" in line:
                sp = line.split()[3]
    return float(sp)
    
def read_deltax(file):
    with open(file) as f: 
        lines = f.readlines()
        for line in lines:
            if "deltaFat" in line:
                sp = line.split()[1]
    return float(sp)    

root     = args.directory
dirlist  = ["005_Gly54H2O_sd_B45_compress_f/0-20ns","007_gly128H2O_sd_B35_compress_f/0-10.5ns"]
savepath ="./"
list1    = ["c100","d100","s100"]
list2    = ["c100","d100","s100"]
list3    = ["Charge_lambda100","Distance_lambda100","Solvation_lambda100"]
for i in range(len(list1)):
    path1  = os.path.join(root,dirlist[0],"fes_reweight",list1[i],"fes-rew.dat")
    path2  = os.path.join(root,dirlist[1],"fes_reweight",list1[i],"fes-rew.dat")
    data1  = plumed.read_as_pandas(path1)
    data2  = plumed.read_as_pandas(path2)
    detf1  = read_deltaf(path1)
    detf2  = read_deltaf(path2)
    detx1  = read_deltax(os.path.join(root,dirlist[0],"fes_reweight",list1[i],"q.pbs"))
    detx2  = read_deltax(os.path.join(root,dirlist[1],"fes_reweight",list1[i],"q.pbs"))
    fes_1D([data1,data2],[detf1,detf2],[detx1,detx2],[list2[i]],["file.free"],["uncertainty"],["54H2O","128H2O"],"CV: %s"%list3[i],"Free energy [kJ/mol]",savepath,1,True)
