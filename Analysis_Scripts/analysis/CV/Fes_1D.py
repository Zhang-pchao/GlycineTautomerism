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
parser.add_argument('--dir','-d',dest='directory',type=str,default='../',help='the directory of all files')
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

def fes_1D(data,xlist,ylist,err,ylabel,Xlabel,Ylabel,path='./',step=1,Legend=True):
    data=data
    fig = plt.figure(figsize=(6,6), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    colors = cmap2colorlist('GnBu', len(xlist))
    for i in range(len(xlist)):
        #ax.scatter(data[xlist[i]][::step],data[ylist[i]][::step],
        #           s=3,alpha=0.7,color=colors[i],label=ylabel[i])
        ax.plot(data[xlist[i]][::step],data[ylist[i]][::step],
                   color=colors[i],label=ylabel[i])
        ax.fill_between(data[xlist[i]][::step],
                        data[ylist[i]][::step]-data[err[i]][::step],
                        data[ylist[i]][::step]+data[err[i]][::step],
                        facecolor = colors[i],alpha = 0.3)
    ax.set_xlabel(Xlabel)
    ax.set_ylabel(Ylabel)
    if Legend:
        ax.legend()
    fig.savefig(os.path.join(path,"%s_%s.png"%(Xlabel,ylabel)), dpi=150, bbox_inches='tight')


root= args.directory
savepath=os.path.join(root,"fes1D")
mkdir(savepath)
#c_data=plumed.read_as_pandas(os.path.join(root,"fes_reweight","TotalCharge_deltaG",  "fes-rew.dat"))
d_data=plumed.read_as_pandas(os.path.join(root,"fes_reweight","IonDistance_deltaG",  "fes-rew.dat"))


#fes_1D(s_data,["cc"],  ["file.free"],["FreeEnergy"],"CV: TotalCharge","Free energy [kJ/mol]",savepath,1,False)
fes_1D(d_data,["logd"],["file.free"],["uncertainty"],["FreeEnergy"],"CV:F(IonDistance)","Free energy [kJ/mol]",savepath,1,False)
