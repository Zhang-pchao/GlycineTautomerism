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
    cb = fig.colorbar(im, ax=ax, label='Free energy [kJ/mol]')
    #ax.set_yticks([-2,0,2])
    ax.set_aspect('auto')
    ax.set_xlabel('CV: SOLVATION')
    ax.set_ylabel('CV: IONDISTANCE')
    #ax.set_title(title)
    return cb

###change below
folder     = "./"
max_fes    = 96 #kJ/mol
pathlist   = ['s05_d05_bin500']
cvnames    = ['s05','d05']
flag       = 'nocharge'
set_bar_ticks = True # modify color bar ticks 
###change above

for k in pathlist:
    X1, fes_opes1, = load_fes(fes_folder=os.path.join(folder,k))
    extent_opes = get_xy_min_max(os.path.join(folder,k,"fes-rew.dat"),cvnames)
    fig, (ax1) = plt.subplots(1, 1)
    fig.set_size_inches((8, 6))
    plot_ala2_fes(ax1, X1, extent_opes,max_fes=max_fes)
    plt.savefig(os.path.join(os.path.join(folder,k),"%s_%s_%s.png"%(k,"fes_reweight","1")),
                dpi=150, bbox_inches='tight')

colors = [(31 ,59 ,115),
          (40 ,111,134), 
          (53 ,152,146),
          (73 ,171,142 ),
          (114,192,118 ),
#          (167,214,85 ),
#          (198,217,83 ),         
          (219,220,71 ),
#          (231,206,74 ),
#          (243,193,77 ),          
          (255,180,80 ),
          (240,160,75 ),          
          (228,120,69 ),
          (200, 80,57 ),          
          (188, 47,46 )]  # R -> G -> B          
colornum=12
colors=list(tuple(i/255 for i in color) for color in colors)
colors.reverse()
cm = LinearSegmentedColormap.from_list('my_list', colors, N=colornum)

for k in pathlist:
    multiplot,axes = plt.subplots(1, 1)
    multiplot.set_size_inches((8, 8))
    ffile = open(os.path.join(folder,k,"fes-analysis.dat"),'w')
    fskip1 = open(os.path.join(folder,k,"fes-analskip1.dat"),'w')
    fskip2 = open(os.path.join(folder,k,"fes-analskip2.dat"),'w')
    
    x,y,z=np.loadtxt(os.path.join(folder,k,"fes-rew.dat"),unpack=True)[:3]
    z=z-z.min()
    z[z>max_fes]=max_fes
    #z[z<11]=11

    yi = np.linspace(min(y), max(y), 500)
    xi = np.linspace(min(x), max(x), 500)
    X,Y=np.meshgrid(xi,yi)
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='nearest')
    ffile.write('#%16s%16s%16s\n'%('x','y','fes'))
    fskip1.write('#%16s%16s%16s\n'%('x','y','fes'))
    fskip2.write('#%16s%16s%16s\n'%('x','y','fes'))
    xx = list(xi)
    yy = list(yi)
    zz = list(zi)
    for i in range(len(yy)):
        for j in range(len(xx)):
            ffile.write(' %16.8f%16.8f%16.8f\n'%(xx[j],yy[i],zz[i][j]))
            #for ani2zwi2can
            if -1<xx[j]<-0.5 and 0<yy[i]<2:
                fskip1.write(' %16.8f%16.8f%16.8f\n'%(xx[j],yy[i],max_fes))
            else:
                fskip1.write(' %16.8f%16.8f%16.8f\n'%(xx[j],yy[i],zz[i][j]))
            #for can2zwi2cat    
            if -1<xx[j]<-0.5 and 0<yy[i]<2:
                fskip2.write(' %16.8f%16.8f%16.8f\n'%(xx[j],yy[i],max_fes))
            elif -0.5<=xx[j]<0.4 and 1<yy[i]<1.2:
                fskip2.write(' %16.8f%16.8f%16.8f\n'%(xx[j],yy[i],max_fes))                
            else:
                fskip2.write(' %16.8f%16.8f%16.8f\n'%(xx[j],yy[i],zz[i][j]))               
    ffile.close()
    fskip1.close()
    fskip2.close()
    
    extent_opes = get_xy_min_max(os.path.join(folder,k,"fes-rew.dat"),cvnames)
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
        mycbar = multiplot.colorbar(mlt, cax=cbar_ax, ticks=[xx*8 for xx in [0,1,2,3,4,5,6,7,8,9,10,11,12]],
                            label=r'$\Delta$'+r'$F$'+': relative free energy (kJ/mol)')
        mycbar.ax.set_yticklabels(["%d"%int(xx*8) for xx in [0,1,2,3,4,5,6,7,8,9,10,11]]+['>96'])                   
    else:       
        multiplot.colorbar(mlt, cax=cbar_ax, label='$Delta$'+r'$F$'+': relative free energy (kJ/mol)')  
    multiplot.savefig(os.path.join(os.path.join(folder,k),"%s_%s_%s.png"%(k,"fes_reweight","2")),
                dpi=600, bbox_inches='tight')