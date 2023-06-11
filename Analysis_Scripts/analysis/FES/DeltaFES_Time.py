#!/usr/bin/env python
# coding: utf-8

import numpy as np
import plumed
import os
import sys

#reuse same plumed kernel, to avoid multiple warnings
PLUMED_KERNEL=plumed.Plumed()
from functools import partial
plumed.read_as_pandas = partial(plumed.read_as_pandas, kernel=PLUMED_KERNEL)

#plot-related stuff
import matplotlib.pyplot as plt
from matplotlib import cm
from IPython.display import clear_output
#set bigger font sizes
SMALL_SIZE = 13
MEDIUM_SIZE = 14
BIG_SIZE = 17
plt.rc('font', size=SMALL_SIZE)        # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)   # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)   # fontsize of the figure title

def cmap2colorlist(cmap_name, color_numb):
    colormap = cm.get_cmap(cmap_name, color_numb+4)
    idx = np.arange(color_numb) + 3
    colorlist = colormap(idx)
    return colorlist

def read_deltaf(file):
    with open(file) as f: 
        lines = f.readlines()
        for line in lines:
            if "DeltaF" in line:
                sp = line.split()[3]
    return float(sp) #kJ/mol
    
def read_time(file):
    with open(file) as f: 
        lines = f.readlines()
        for line in lines:
            if "total_sample_size" in line:
                sp = float(line.split()[3])/1e6
    return sp #ns

def read_cvname(file):
    with open(file) as f: 
        lines = f.readlines()
        for line in lines:
            if "FIELDS" in line:
                sp = line.split()[2]
    return sp

def read_cvrange(file):
    with open(file) as f: 
        lines = f.readlines()
        for line in lines:
            if "python" in line:
                sp = line.split()
                for i in range(len(sp)):
                    if sp[i] == '--min':
                        mins = float(sp[i+1])
                    if sp[i] == '--max':
                        maxs = float(sp[i+1])
                    if sp[i] == '--deltaFat':
                        delta = float(sp[i+1])                           
    return [mins,maxs,delta]

def cv_fig(x,y,cvname,choosecv,path):
    fig = plt.figure(figsize=(6,6), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    colors = cmap2colorlist('GnBu', 4)
    ax.plot(x,y,lw=2,color=colors[0],label='%s[%4.2f~%4.2f] - %s[%4.2f~%4.2f]'%(cvname,choosecv[1],choosecv[2],
                                                                                cvname,choosecv[0],choosecv[2]))
    ax.legend()
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('$\Delta F$ (kJ/mol)')
    fig.savefig(os.path.join(path,"DeltaFES_Time.png"), dpi=150, bbox_inches='tight')


#if __name__ == '__main__':
root='../'
data=plumed.read_as_pandas(os.path.join(root,"COLVAR_tmp"))
savepath=os.path.join(root,'deltafes')

step = 15
T = []
F = []
for i in range(step):
    T.append(read_time(os.path.join(root,'deltafes',str(10001+i),'fes-rew.dat')))
    F.append(read_deltaf(os.path.join(root,'deltafes',str(10001+i),'fes-rew.dat')))

cvname = read_cvname(os.path.join(root,'deltafes',str(10001),'fes-rew.dat'))
choosecv = read_cvrange(os.path.join(root,'deltafes',str(10001),'q.sub'))
cv_fig(T,F,cvname,choosecv,savepath)