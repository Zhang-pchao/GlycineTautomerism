#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.patches as patches
import matplotlib.patheffects as path_effects
import matplotlib.image as mpimg
from matplotlib.pyplot import MultipleLocator
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
import matplotlib.ticker

import os
import sys
module_dir='/home/pengchao/jupy/mymodules/colorbar'
sys.path.append(os.path.abspath(module_dir))
import fessa

import plumed
import argparse

#reuse same plumed kernel, to avoid multiple warnings
PLUMED_KERNEL=plumed.Plumed()
from functools import partial
plumed.read_as_pandas = partial(plumed.read_as_pandas, kernel=PLUMED_KERNEL)
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
plt.rc('legend', fontsize=SMALL_SIZE-3)  # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)   # fontsize of the figure title

path0='/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/020_1e-6_add_water/008_gly128H2O_sd_B35_compress_f/0-30ns'
#path0='/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/020_1e-6_add_water/008_gly128H2O_sd_B35_compress_f/22.9-30ns'
folder_cv   = os.path.join(path0,"COLVAR_tmp")
save_path=os.path.join(path0,'final_figs')

cvdata = plumed.read_as_pandas(folder_cv)
cvstep = 10

fig = plt.figure(figsize=(18,8), dpi=150, facecolor='white')
grid = plt.GridSpec(46, 37, wspace=0.5, hspace=0.5)
ax11  = fig.add_subplot(grid[ 0:22,  0:36])
ax12  = fig.add_subplot(grid[24:46,  0:36])
axcc  = fig.add_subplot(grid[ 0:46, 36:37])
#ax1c  = fig.add_subplot(grid[ 0:22, 36:37])
#ax2c  = fig.add_subplot(grid[24:46, 36:37])
lll=['a','b']
AX=[ax11,ax12]
for i in range(2):
    AX[i].text(0.01, 0.98, '(%s)'%lll[i], transform=AX[i].transAxes, fontsize=20,verticalalignment='top')

choose1  = 500000 # skip 500ps
choose2  = -1
mylambda = 5
mybias = cvdata["opes.bias"][choose1:choose2:cvstep]
min_t  = min(mybias)
max_t  = max(mybias)
# Map your bias values to a range from 0-1 using a Normalize object
norm = matplotlib.colors.Normalize(vmin=min_t, vmax=max_t)
color = [matplotlib.cm.viridis(norm(b)) for b in mybias]
# Define your Formatter function
fmt = matplotlib.ticker.FuncFormatter(lambda x,pos:'%.1f' %(x*(max_t-min_t)+min_t))
# Create the scatter plot and colorbar
im = ax11.scatter(cvdata['time'][choose1:choose2:cvstep]/1e3, cvdata['s05'][choose1:choose2:cvstep],
                  s=1, c=color, marker='.')
fig.colorbar(im, cax=axcc, format=fmt, label='Bias (kJ/mol)',cmap="viridis")

ax11.set_xlim([-0.5,30.5])
ax11.set_ylim([-1.1,1.1])
#ax11.set_xlabel('Simulation time (ns)',fontsize='20')
ax11.set_ylabel(r'$CV_p$'+r'$^{\lambda=%d}$'%mylambda,fontsize='20')
ax11.xaxis.set_major_locator(MultipleLocator(5))
ax11.xaxis.set_minor_locator(MultipleLocator(5/5))
ax11.yaxis.set_major_locator(MultipleLocator(0.5))
ax11.yaxis.set_minor_locator(MultipleLocator(0.5/5))

ax12.set_xlim([-0.5,30.5])
ax12.set_ylim([-1,12.5])
im = ax12.scatter(cvdata['time'][choose1:choose2:cvstep]/1e3,cvdata['d05'][choose1:choose2:cvstep],
                  s=1, c=color, marker='.')    
#fig.colorbar(im, cax=axcc,format = fmt,label='Bias (kJ/mol)',cmap="viridis") 
ax12.set_xlabel('Simulation time (ns)',fontsize='20')
ax12.set_ylabel(r'$CV_d$'+r'$^{\lambda=%d}$'%mylambda,fontsize='20')
ax12.xaxis.set_major_locator(MultipleLocator(5))
ax12.xaxis.set_minor_locator(MultipleLocator(5/5))
ax12.yaxis.set_major_locator(MultipleLocator(2.5))
ax12.yaxis.set_minor_locator(MultipleLocator(2.5/5))

fig.savefig(os.path.join(save_path, "cv_bias_SI.png"), dpi=600, bbox_inches='tight')