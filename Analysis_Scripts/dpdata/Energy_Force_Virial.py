#!/usr/bin/env python
# coding: utf-8

import plumed
import matplotlib.pyplot as plt
import os
import MDAnalysis
import numpy as np
import math
import sys
import re
import glob
import shutil

#reuse same plumed kernel, to avoid multiple warnings
PLUMED_KERNEL=plumed.Plumed()
from functools import partial
plumed.read_as_pandas = partial(plumed.read_as_pandas, kernel=PLUMED_KERNEL)

import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from matplotlib.offsetbox import AnchoredText

atoms = 172
save_path = './'


def cmap2colorlist(cmap_name, color_numb):
    """gets a list of color from selected matplotlib cmap
    'cmap_name': str, name of cmap, for example 'plasma'.
                 See available cmaps in tutorial "Choosing Colormaps in Matplotlib"
    'color_numb': int, length of desired color list
    Returns to a numpy array. Its shape is (color_numb, 4),
    and each column (a (4,) vector) represents a color.
    example:
        clist = cmap2colorlist('plasma', 4)
        for ii in range(4):
            ax.plot(x, y[ii], color=clist[ii])
    """
    colormap = cm.get_cmap(cmap_name, color_numb+4)
    idx = np.arange(color_numb) + 3
    colorlist = colormap(idx)
    return colorlist

def min_max(ef_list1,ef_list2,bound):
    final_min = min(int(min(ef_list1))-bound,int(min(ef_list2))-bound)
    final_max = max(int(max(ef_list1))+bound,int(max(ef_list2))+bound)
    return final_min, final_max

def replace_txt(filename):
    f_path = 'detail.%s.out' %filename
    searchfile = open(f_path, "r")
    lines = searchfile.readlines()
    searchfile.close()
    f = open (f_path, "r+")
    if filename == 'e':
        replace_line = '#! FIELDS data_e pred_e\n'
    elif filename == 'f':
        replace_line = '#! FIELDS  data_fx data_fy data_fz pred_fx pred_fy pred_fz\n'
    else:
        replace_line = '#! FIELDS data_vxx data_vxy data_vxz data_vyx data_vyy data_vyz data_vzx data_vzy data_vzz pred_vxx pred_vxy pred_vxz pred_vyx pred_vyy pred_vyz pred_vzx pred_vzy pred_vzz\n'        
    open('detail.%s.out.2' %filename, 'w').write(re.sub(lines[0], replace_line, f.read()))

replace_txt('e')
replace_txt('f')
replace_txt('v')

F = ['data_fx', 'data_fy', 'data_fz', 'pred_fx', 'pred_fy', 'pred_fz']
V = ['data_vxx', 'data_vxy', 'data_vxz', 'data_vyx', 'data_vyy', 'data_vyz', 'data_vzx', 'data_vzy', 'data_vzz', 'pred_vxx', 'pred_vxy', 'pred_vxz', 'pred_vyx', 'pred_vyy', 'pred_vyz', 'pred_vzx', 'pred_vzy', 'pred_vzz']

colors = cmap2colorlist('GnBu', 4)
fig = plt.figure(figsize=(17,5), dpi=600, facecolor='white')    

#Energy
ax = fig.add_subplot(131)
data_e=plumed.read_as_pandas("detail.e.out.2")
ax.scatter(data_e["data_e"], data_e["pred_e"], alpha=1, label = "Energy", color=colors[0])
ax.set_xlabel("Energy calculated by DFT (eV)", fontsize = 14)
ax.set_ylabel("Energy predicted by DP (eV)", fontsize = 14)
ax.tick_params(axis = 'both', which = 'major', labelsize = 12)    
emin,emax = min_max(data_e["data_e"], data_e["pred_e"],0.5)    
e_locator = int((emax-emin)/4)
ax.xaxis.set_major_locator(MultipleLocator(e_locator))
ax.yaxis.set_major_locator(MultipleLocator(e_locator))   
ax.plot([emin,emax], [emin,emax], ls="--", c=".3")
ax.set(xlim=(emin,emax), ylim=(emin,emax))
at = AnchoredText("(a)", prop=dict(size=15), frameon=False, loc='upper left')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)    
    
#Force
ax = fig.add_subplot(132) 
data_f=plumed.read_as_pandas("detail.f.out.2")
fmin,fmax = 0.,0.
for i in range(3):
    ax.scatter(data_f[F[i]], data_f[F[i+3]], alpha=1, label = "Force", color=colors[1])
    tmpmin,tmpmax = min_max(data_f[F[i]], data_f[F[i+3]],0.8)
    if tmpmin < fmin:
        fmin = tmpmin
    if tmpmax > fmax:
        fmax = tmpmax           
ax.set_xlabel("Force calculated by DFT (eV/"+r"$ \rm \AA)$", fontsize = 14)
ax.set_ylabel("Force predicted by DP (eV/"+r"$ \rm \AA)$", fontsize = 14)
ax.tick_params(axis = 'both', which = 'major', labelsize = 12)  
f_locator  = int((fmax-fmin)/6)
ax.xaxis.set_major_locator(MultipleLocator(f_locator))
ax.yaxis.set_major_locator(MultipleLocator(f_locator))        
ax.plot([fmin,fmax], [fmin,fmax], ls="--", c=".3")
ax.set(xlim=(fmin,fmax), ylim=(fmin,fmax))
at = AnchoredText("(b)", prop=dict(size=15), frameon=False, loc='upper left')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)       
    
#Virial
ax = fig.add_subplot(133) 
data_v=plumed.read_as_pandas("detail.v.out.2")
vmin,vmax = 0.,0.
for i in range(9):
    ax.scatter(data_v[V[i]], data_v[V[i+9]], alpha=1, label = "Virial", color=colors[2])
    tmpmin,tmpmax = min_max(data_v[V[i]], data_v[V[i+9]],2)
    if tmpmin < vmin:
        vmin = tmpmin
    if tmpmax > vmax:
        vmax = tmpmax 
ax.set_xlabel("Virial calculated by DFT (eV)", fontsize = 14)
ax.set_ylabel("Virial predicted by DP (eV)", fontsize = 14)
ax.tick_params(axis = 'both', which = 'major', labelsize = 12)  
v_locator  = int((vmax-vmin)/6)
ax.xaxis.set_major_locator(MultipleLocator(v_locator))
ax.yaxis.set_major_locator(MultipleLocator(v_locator))        
ax.plot([vmin,vmax], [vmin,vmax], ls="--", c=".3")
ax.set(xlim=(vmin,vmax), ylim=(vmin,vmax))
at = AnchoredText("(c)", prop=dict(size=15), frameon=False, loc='upper left')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)  

fig.savefig(os.path.join(save_path, "E_F_V.png"), dpi = 600, bbox_inches='tight')

all_data_f = np.array([])
all_pred_f = np.array([])
for i in range(3):
    all_data_f = np.hstack((all_data_f,data_f[F[i]].to_numpy()))
    all_pred_f = np.hstack((all_pred_f,data_f[F[i+3]].to_numpy()))

all_data_v = np.array([])
all_pred_v = np.array([])
for i in range(9):
    all_data_v = np.hstack((all_data_v,data_v[V[i]].to_numpy()))
    all_pred_v = np.hstack((all_pred_v,data_v[V[i+9]].to_numpy()))

fig = plt.figure(figsize=(10,10), dpi=600, facecolor='white')    
grid = plt.GridSpec(17, 13, wspace=0.5, hspace=0.5)
label_fontsize = 14
###Energy###
ax = fig.add_subplot(grid[0:5,0:9])
ax.scatter(data_e["data_e"], abs(data_e["data_e"]-data_e["pred_e"])/atoms*1000, alpha=0.5, color=colors[0],s=10)
ax.set_xlabel('$\it{E}$'+'$\mathregular{_{DFT}}$' + " (eV)", fontsize = label_fontsize)
ax.set_ylabel('|$\it{E}$'+'$\mathregular{_{DP}}$'+ " – " + '$\it{E}$'+'$\mathregular{_{DFT}|}$' + " (meV/atom)", fontsize = label_fontsize)
ax.tick_params(axis = 'both', which = 'major', labelsize = label_fontsize-2)
ax.xaxis.set_major_locator(MultipleLocator(e_locator))
ax.set(ylim=(0,2))
at = AnchoredText("(a)", prop=dict(size=label_fontsize), frameon=False, loc='upper left')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)    

ax = fig.add_subplot(grid[0:5,9:13])
ax.hist(abs(data_e["data_e"]-data_e["pred_e"])/atoms*1000,100,orientation='horizontal',color=colors[0])
ax.tick_params(axis = 'both', which = 'major', labelsize = label_fontsize-2)
ax.set_xlabel("Count", fontsize = label_fontsize)
ax.set(ylim=(0,2))
ax.axes.yaxis.set_ticklabels([])
###Force###
ax = fig.add_subplot(grid[6:11,0:9])
ax.scatter(all_data_f, abs(all_data_f-all_pred_f), alpha=0.5, color=colors[1],s=10)
ax.set_xlabel('$\it{F}$'+'$\mathregular{_{DFT}}$'+" (eV/"+r"$ \rm \AA)$", fontsize = label_fontsize)
ax.set_ylabel('|$\it{F}$'+'$\mathregular{_{DP}}$'+ " – " +'$\it{F}$'+ '$\mathregular{_{DFT}|}$' +" (eV/"+r"$ \rm \AA)$", fontsize = label_fontsize)
ax.tick_params(axis = 'both', which = 'major', labelsize = label_fontsize-2)
ax.xaxis.set_major_locator(MultipleLocator(f_locator))
ax.set(ylim=(0,0.2))  
at = AnchoredText("(b)", prop=dict(size=label_fontsize), frameon=False, loc='upper left')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)       

ax = fig.add_subplot(grid[6:11,9:13])
ax.hist(abs(all_data_f-all_pred_f),5000,orientation='horizontal', color=colors[1])
ax.tick_params(axis = 'both', which = 'major', labelsize = label_fontsize-2)
ax.set_xlabel("Count", fontsize = label_fontsize)
ax.set(ylim=(0,0.2))
ax.axes.yaxis.set_ticklabels([])
###Virial###
ax = fig.add_subplot(grid[12:17,0:9])
ax.scatter(all_data_v, abs(all_data_v-all_pred_v)/atoms*1000, alpha=0.5, color=colors[2],s=10)
ax.set_xlabel('$\it{V}$'+'$\mathregular{_{DFT}}$' + " (eV)", fontsize = label_fontsize)
ax.set_ylabel('|$\it{V}$'+'$\mathregular{_{DP}}$'+ " – " + '$\it{V}$'+'$\mathregular{_{DFT}|}$' + " (meV/atom)", fontsize = label_fontsize)
ax.tick_params(axis = 'both', which = 'major', labelsize = label_fontsize-2)
ax.xaxis.set_major_locator(MultipleLocator(v_locator))
ax.set(ylim=(0,10))
at = AnchoredText("(c)", prop=dict(size=label_fontsize), frameon=False, loc='upper left')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)    

ax = fig.add_subplot(grid[12:17,9:13])
ax.hist(abs(all_data_v-all_pred_v)/atoms*1000,100,orientation='horizontal',color=colors[2])
ax.tick_params(axis = 'both', which = 'major', labelsize = label_fontsize-2)
ax.set_xlabel("Count", fontsize = label_fontsize)
ax.set(ylim=(0,10))
ax.axes.yaxis.set_ticklabels([])    

fig.savefig(os.path.join(save_path, "E_F_V_hist.png"), dpi = 600, bbox_inches='tight')

os.remove("detail.e.out.2")
os.remove("detail.f.out.2")
os.remove("detail.v.out.2")