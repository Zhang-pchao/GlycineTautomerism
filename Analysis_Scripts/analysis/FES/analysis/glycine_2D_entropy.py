#!/usr/bin/env python
# coding: utf-8

import numpy as np
import sys
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.patches as patches
import matplotlib.patheffects as path_effects

import os
import plumed
import argparse

#reuse same plumed kernel, to avoid multiple warnings
PLUMED_KERNEL=plumed.Plumed()
from functools import partial
plumed.read_as_pandas = partial(plumed.read_as_pandas, kernel=PLUMED_KERNEL)

#some other plotting functions
from matplotlib.pyplot import MultipleLocator
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from IPython.display import clear_output

#set bigger font sizes
SMALL_SIZE = 16
MEDIUM_SIZE = 17
BIG_SIZE = 20
plt.rc('font', size=SMALL_SIZE)        # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)   # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)   # fontsize of the figure title

def load_fes(fes_folder):
    filename = os.path.join(fes_folder,'fes-rew.dat')
    X = plumed.read_as_pandas(filename, usecols=[2]).to_numpy()
    nbins = int(np.sqrt(len(X)))
    fes = []
    return X.reshape(nbins,nbins), fes

def get_xy_min_max(filename,cvs):
    with open(filename) as f: 
        lines = f.readlines()
        for line in lines:
            if "SET min_%s"%cvs[0] in line:
                minss = float(line.split()[3])
            if "SET max_%s"%cvs[0] in line:
                maxss = float(line.split()[3])    
            if "SET min_%s"%cvs[1] in line:
                mindd = float(line.split()[3])
            if "SET max_%s"%cvs[1] in line:
                maxdd = float(line.split()[3])
    return (minss,maxss,mindd,maxdd)

def plot_ala2_fes(ax, fes, myextent,max_fes=None):
    if max_fes is None:
        max_fes = np.amax(fes)
    im = ax.imshow(fes, vmax=max_fes, origin='lower', extent=myextent)
    cb = fig.colorbar(im, ax=ax, label='Relative entropy(kJ/mol)')
    #ax.set_yticks([-2,0,2])
    ax.set_aspect('auto')
    ax.set_xlabel('CV: SOLVATION')
    ax.set_ylabel('CV: IONDISTANCE')
    #ax.set_title(title)
    return cb

###change below
folder_f     = "../fes2D"
folder_p     = "../potential2D"
folder_e     = "../entropy2D"
max_fes    = 90 #kJ/mol
min_ent    = 4170 #kJ/mol
pathlist   = ['s05_d05_bin20']
cvnames    = ['s05','d05']
flag       = 'nocharge'
Temp       = 300.0
set_bar_ticks = True # modify color bar ticks 
###change above

colors = [(255,255,255), #white
#         (31 ,59 ,115),
          (40 ,111,134), 
          (53 ,152,146),
          (73 ,171,142 ),
          (114,192,118 ),
          (167,214,85 ),
          (219,220,71 ),
          (255,180,80 ),
          (228,120,69 ),          
          (188, 47,46 )]  # R -> G -> B          
colornum=10
colors=list(tuple(i/255 for i in color) for color in colors)
#colors.reverse()
cm = LinearSegmentedColormap.from_list('my_list', colors, N=colornum)

for k in pathlist:
    multiplot,axes = plt.subplots(1, 1)
    multiplot.set_size_inches((8, 8))
    
    x_f,y_f,z_f=np.loadtxt(os.path.join(folder_f,k,"fes-rew.dat"),unpack=True)[:3]
    x_p,y_p,z_p=np.loadtxt(os.path.join(folder_p,k,"fes-rew.dat"),unpack=True)[:3]
    #z=z-z.min()
    #z_f[z_f>max_fes]=max_fes
    #z_p[z_p==-1]=z_p.max()
    #z=z/Temp #Temprature
    #z[z>max_fes]=max_fes
    #z[z<min_fes]=min_fes
    
    if x_f[int(len(x_f)/2)] == x_p[int(len(x_p)/2)] and y_f[int(len(y_f)/2)] == y_p[int(len(y_p)/2)]:
        x=x_f
        y=y_f
        z=z_p-z_f
        min_tmp=10000
        for i in range(len(x_f)):
            if z_p[i]!=-1 and z_f[i]< max_fes:
                if z[i]<min_tmp:
                    min_tmp=z[i]
        for i in range(len(x_f)):
            if z_p[i]==-1 or z_f[i]>= max_fes:
                z[i]=-15
            else:
                z[i]-=min_tmp
                    
                    
    yi = np.linspace(min(y), max(y), 500)
    xi = np.linspace(min(x), max(x), 500)
    X,Y=np.meshgrid(xi,yi)

    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='nearest')
    
    extent_opes = get_xy_min_max(os.path.join(folder_f,k,"fes-rew.dat"),cvnames)
    axes.set_xlim([extent_opes[0],extent_opes[1]])
    axes.set_ylim([extent_opes[2],extent_opes[3]])

    axes.contour(xi, yi, zi, colornum, colors='k', linewidths=0.8)
    mlt=axes.contourf(xi, yi, zi, colornum,cmap=cm,)
    axes.xaxis.set_major_locator(MultipleLocator(0.5))

    axes.set_xlabel(r'$s_p$'+': glycine protonation',fontsize='20')
    if flag == 'charge':
        axes.set_ylabel(r'$s_c$'+': ion charge',fontsize='20')
        axes.yaxis.set_major_locator(MultipleLocator(0.5))
    else:
        axes.set_ylabel(r'$s_d$'+': ion distance'+" ("+r"$ \rm \AA)$",fontsize='20')
        axes.yaxis.set_major_locator(MultipleLocator(2))

    cbar_ax = multiplot.add_axes([0.92, 0.11, 0.018, 0.77])
    if set_bar_ticks:
        mycbar = multiplot.colorbar(mlt, cax=cbar_ax, ticks=[x*15 for x in [-1,0,1,2,3,4,5,6,7]],
                            label=r'$T$'+r'$\Delta$'+r'$S$'+': relative entropy at 300 K (kJ/mol)')
        mycbar.ax.set_yticklabels([' ']+['%s'%str(x*15) for x in [0,1,2,3,4,5,6,7]])                   
    else:
        multiplot.colorbar(mlt, cax=cbar_ax, label=r'$T$'+r'$\Delta$'+r'$S$'+': relative entropy at 300 K (kJ/mol)')
    multiplot.savefig(os.path.join(os.path.join(folder_e,k),"%s_%s_%s.png"%(k,"fes_reweight","2")),
                dpi=600, bbox_inches='tight')