#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import plumed
import sys
import argparse
#reuse same plumed kernel, to avoid multiple warnings
PLUMED_KERNEL=plumed.Plumed()
from functools import partial
plumed.read_as_pandas = partial(plumed.read_as_pandas, kernel=PLUMED_KERNEL)

### Parser stuff ###
parser = argparse.ArgumentParser(description='get the FES estimate used by OPES')
# files
parser.add_argument('--dir','-d',dest='directory',type=str,default='../',help='the directory of all files')
# some easy parsing
args=parser.parse_args()


#plot-related stuff
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from IPython.display import clear_output
import matplotlib.colors 
import matplotlib.ticker
from matplotlib.pyplot import MultipleLocator
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
    
def cv_fig(data,xlist,ylist,ylabel,Xlabel,Ylabel,path='./',step=1,Legend=True):
    data=data
    fig = plt.figure(figsize=(6,6), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    colors = cmap2colorlist('GnBu', len(xlist))
    for i in range(len(xlist)):
        ax.scatter(data[xlist[i]][::step],data[ylist[i]][::step],
                   s=1,alpha=0.7,color=colors[i],label=ylabel[i])
    ax.set_xlabel(Xlabel)
    ax.set_ylabel(Ylabel)
    if Legend:
        ax.legend(markerscale=5)
    fig.savefig(os.path.join(path,"%s_%s.png"%(Xlabel,Ylabel)), dpi=150, bbox_inches='tight')

def cv_colorbar_fig(data,xlist,ylist,ylabel,Xlabel,Ylabel,bias,path='./',step=1,Legend=True):
    data=data
    fig = plt.figure(figsize=(8,6), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    colors = cmap2colorlist('GnBu', len(xlist))
    
    mybias = data[bias[0]][::step]
    min_t = min(mybias)
    max_t = max(mybias)
    color = [plt.get_cmap("viridis", 100)(int((i-min_t)/(max_t-min_t)*100)) for i in mybias]
     
    # colorbar
    im = ax.scatter(data[xlist[0]][::step],data[ylist[0]][::step], s=2, c=color, marker='.')    
    fmt = matplotlib.ticker.FuncFormatter(lambda x,pos:'%.1f' %(x*(max_t-min_t)+min_t))
    fig.colorbar(im, shrink=1.0,format = fmt,label='Bias (kJ/mol)')    

    ax.set_xlabel(Xlabel)
    ax.set_ylabel(Ylabel)    
    if Legend:
        ax.legend(markerscale=5)
    fig.savefig(os.path.join(path,"%s_%s.png"%(Xlabel,Ylabel)), dpi=150, bbox_inches='tight')
    
root= args.directory
data=plumed.read_as_pandas(os.path.join(root,"COLVAR"))
savepath=os.path.join(root,'cvfig')
mkdir(savepath)

cv_colorbar_fig(data,["time"],["ss"], ["ss"],  "Time (ps)","CV:Solvation with bias", ["opes.bias"], savepath,10,False)

cv_fig(data,["ss"],["dd"],["dd"],"CV:Solvation","CV:Iondistance",savepath,1,False)
cv_fig(data,["ss"],["cc"],["cc"],"CV:Solvation","CV:Totalcharge",savepath,1,False)
cv_fig(data,["time"],["ss"], ["ss"],  "Time (ps)","CV:Solvation",  savepath,1,True)
cv_fig(data,["time"],["dd"], ["dd"],  "Time (ps)","CV:Iondistance",savepath,1,True)
cv_fig(data,["time"],["cc"], ["Q^2_H2O"],    "Time (ps)","No. of total charge",       savepath,1,True)
cv_fig(data,["time","time"],["cp","cm"], ["Q_H3O+","Q_OH-"],    "Time (ps)","No. of charge",       savepath,1,True)
