#!/usr/bin/env python
# coding: utf-8

import os
import sys
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.pyplot import MultipleLocator
import numpy as np
import glob
import scipy.signal

def cmap2colorlist(cmap_name, color_numb):
    colormap = cm.get_cmap(cmap_name, color_numb+4)
    idx = np.arange(color_numb) + 3
    colorlist = colormap(idx)
    return colorlist

def get_data(f):
    x = np.loadtxt(f,usecols=0)
    y = np.loadtxt(f,usecols=1)
    return x,y

###change below###
path = '/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/020_1e-6_add_water/rdf'
rdfs = ['1N-allH.dat','1N-allO.dat','2O-allH.dat','2O-allO.dat']
grs  = [["N","H"],["N","O"],["O","H"],["O","O"]]
dirs = ['008_1canonical_360H2O','007_1canonical_256H2O','006_1canonical_128H2O','005_1canonical_58H2O']
nums = [360,256,128,54]
idx  = 0
###change above###

X = []
Y = []
for i in range(len(dirs)):
    x,y = get_data(os.path.join(path,dirs[i],rdfs[idx]))
    X.append(x)
    Y.append(y)

#colors = cmap2colorlist('GnBu', len(Y))
colors = [(31 ,59 ,115), 
          (47 ,146,148), 
#          (80 ,178,141),
#          (167,214,85 ),
#          (255,224,62 ),
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

fig = plt.figure(figsize=(8,3), dpi=150, facecolor='white')
ax = fig.add_subplot(111)
ax.grid(True)
label_h2o = '$\mathregular{H_{2}O}$'

for i in range(len(Y)):
    ax.plot(X[i], Y[i], alpha=0.8, color=colors[i], lw=3, label="%4d%4s"%(nums[i],label_h2o))

ax.set_title("%s(glycine) - %s(all)"%(grs[idx][0],grs[idx][1]))
ax.legend(loc="upper right")
ax.set_xlim(0,12)
ax.set_ylim(0,2)
ax.set_xlabel("r "+r"$\ \rm (\AA)$", fontsize = label_fontsize)
ax.set_ylabel("g(r)", fontsize = label_fontsize)
savepath = path
fig.savefig(os.path.join(savepath, "RDF_%s(glycine)_%s(all).png"%(grs[idx][0],grs[idx][1])), dpi=600, bbox_inches='tight')
