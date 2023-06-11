#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import os
import numpy as np
import math
import sys
import re
import glob
import shutil

import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from matplotlib.offsetbox import AnchoredText

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

def min_max(ef_list1,ef_list2,bound):
    final_min = min(int(min(ef_list1))-bound,int(min(ef_list2))-bound)
    final_max = max(int(max(ef_list1))+bound,int(max(ef_list2))+bound)
    return final_min, final_max

def check_npy(loadData):
    print(type(loadData))
    print(loadData.dtype)
    print(loadData.ndim)
    print(loadData.shape)
    print(loadData.size)
    print(loadData)    

def plotfig(data_e,pred_e,data_fx,pred_fx,path):
    colors = cmap2colorlist('GnBu', 4)
    fig = plt.figure(figsize=(16,7.5), dpi=600, facecolor='white')    

    ax = fig.add_subplot(121)
    ax.scatter(data_e, pred_e, alpha=1, label = "The relative energy", color=colors[0])
    ax.legend(prop={'size': 15},loc="lower right")
    ax.set_xlabel("Energy from cp2k (eV)", fontsize = 20)
    ax.set_ylabel("Energy from abacus+deepks (eV)", fontsize = 20)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 16)    
    emin,emax = min_max(data_e,pred_e,1)    
    e_locator = int((emax-emin)/4)
    ax.xaxis.set_major_locator(MultipleLocator(e_locator))
    ax.yaxis.set_major_locator(MultipleLocator(e_locator))   
    ax.plot([emin,emax], [emin,emax], ls="--", c=".3")
    ax.set(xlim=(emin,emax), ylim=(emin,emax))
    at = AnchoredText("(a)", prop=dict(size=15), frameon=True, loc='upper left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)
    
    ax = fig.add_subplot(122) 
    #all of the xyz force. fx is for label
    ax.scatter(data_fx, pred_fx, alpha=1, label = "The atomic force", color=colors[1])
    ax.legend(prop={'size': 15},loc="lower right")
    ax.set_xlabel("Force from cp2k (eV/"+r"$ \rm \AA)$", fontsize = 20)
    ax.set_ylabel("Force from abacus+deepks (eV/"+r"$ \rm \AA)$", fontsize = 20)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 16)  
    fxmin,fxmax = min_max(data_fx,pred_fx,0.8)
    fx_locator  = int((fxmax-fxmin)/6)
    ax.xaxis.set_major_locator(MultipleLocator(fx_locator))
    ax.yaxis.set_major_locator(MultipleLocator(fx_locator))        
    ax.plot([fxmin,fxmax], [fxmin,fxmax], ls="--", c=".3")
    ax.set(xlim=(fxmin,fxmax), ylim=(fxmin,fxmax))
    at = AnchoredText("(b)", prop=dict(size=15), frameon=True, loc='upper left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at) 
    
    fig.savefig(os.path.join(path, "Energy_Force.png"), dpi = 600, bbox_inches='tight')

def plotfig2(data_e,pred_e,data_f,pred_f,path,atoms):
    fig = plt.figure(figsize=(10,6.5), dpi=600, facecolor='white')    
    grid = plt.GridSpec(11, 13, wspace=0.5, hspace=0.5)
    label_fontsize = 14
    colors = cmap2colorlist('GnBu', 4)
    ###Energy###
    ax = fig.add_subplot(grid[0:5,0:9])
    ax.scatter(data_e, np.absolute(data_e-pred_e)/atoms*1000, alpha=0.5, color=colors[0],s=10)
    ax.set_xlabel('$\it{E}$'+'$\mathregular{_{DFT}}$' + " (eV)", fontsize = label_fontsize)
    ax.set_ylabel('|$\it{E}$'+'$\mathregular{_{DPKS}}$'+ " – " + '$\it{E}$'+'$\mathregular{_{DFT}|}$' + " (meV/atom)", fontsize = label_fontsize)
    ax.tick_params(axis = 'both', which = 'major', labelsize = label_fontsize-2)
    emin,emax = min_max(data_e,pred_e,1)
    e_locator = int((emax-emin)/4)
    ax.xaxis.set_major_locator(MultipleLocator(e_locator))
    ax.set(ylim=(0,2))
    at = AnchoredText("(a)", prop=dict(size=label_fontsize), frameon=False, loc='upper left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)    
    
    ax = fig.add_subplot(grid[0:5,9:13])
    ax.hist(np.absolute(data_e-pred_e)/atoms*1000,20,orientation='horizontal',color=colors[0])
    ax.tick_params(axis = 'both', which = 'major', labelsize = label_fontsize-2)
    ax.set_xlabel("Count", fontsize = label_fontsize)
    ax.set(ylim=(0,2))
    ax.axes.yaxis.set_ticklabels([])
    
    ###Force###
    ax = fig.add_subplot(grid[6:11,0:9])
    ax.scatter(data_f, np.absolute(data_f-pred_f), alpha=0.5, color=colors[1],s=10)
    ax.set_xlabel('$\it{F}$'+'$\mathregular{_{DFT}}$'+" (eV/"+r"$ \rm \AA)$", fontsize = label_fontsize)
    ax.set_ylabel('|$\it{F}$'+'$\mathregular{_{DPKS}}$'+ " – " +'$\it{F}$'+ '$\mathregular{_{DFT}|}$' +" (eV/"+r"$ \rm \AA)$", fontsize = label_fontsize)
    ax.tick_params(axis = 'both', which = 'major', labelsize = label_fontsize-2)
    fmin,fmax = min_max(data_f,pred_f,0.8)
    f_locator  = int((fmax-fmin)/6)
    ax.xaxis.set_major_locator(MultipleLocator(f_locator))
    ax.set(ylim=(0,0.2))  
    at = AnchoredText("(b)", prop=dict(size=label_fontsize), frameon=False, loc='upper left')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)       
    
    ax = fig.add_subplot(grid[6:11,9:13])
    ax.hist(np.absolute(data_f-pred_f),1000,orientation='horizontal', color=colors[1])
    ax.tick_params(axis = 'both', which = 'major', labelsize = label_fontsize-2)
    ax.set_xlabel("Count", fontsize = label_fontsize)
    ax.set(ylim=(0,0.2))
    ax.axes.yaxis.set_ticklabels([])
    
    fig.savefig(os.path.join(path, "E_F_hist.png"), dpi = 600, bbox_inches='tight')

def get_error(y_label,y_predict):
    mae = np.sum(np.abs(y_label-y_predict))/len(y_label)
    mse= np.sum((y_label-y_predict)**2)/len(y_label)
    rmse=np.sqrt(mse)
    return mae, rmse

calc_path = "/home/pengchao/cp2k/GlycineH2O/Glycine54H2O/M062X_Dataset/dp_dataset_newmap/"
dpks_path = "/data/HOME_BACKUP/pengchao/glycine/deepks/glycine/pbe_2/test/dp_dataset/"
frame_idx = "008"
atom_nums = 172
errfile   = os.path.join("/data/HOME_BACKUP/pengchao/glycine/deepks/glycine/pbe_2/test",
                         'e_f_test_fig',frame_idx,'mas_rmse.data')
savepath  = os.path.join("/data/HOME_BACKUP/pengchao/glycine/deepks/glycine/pbe_2/test",
                         'e_f_test_fig',frame_idx)
mkdir(savepath)

calc_e = np.load(os.path.join(calc_path,frame_idx,'set.000','energy.npy'))
dpks_e = np.load(os.path.join(dpks_path,frame_idx,'set.000','energy.npy'))
calc_f = np.load(os.path.join(calc_path,frame_idx,'set.000','force.npy'))
dpks_f = np.load(os.path.join(dpks_path,frame_idx,'set.000','force.npy'))

data_e_abs = np.reshape(calc_e,calc_e.size)
pred_e_abs = np.reshape(dpks_e,dpks_e.size)
data_f = np.reshape(calc_f,calc_f.size) 
pred_f = np.reshape(dpks_f,dpks_f.size)

min_e  = np.min(data_e_abs)
data_e = data_e_abs-min_e
pred_e = pred_e_abs-min_e

#check_npy(data_e)

plotfig(data_e,pred_e,data_f,pred_f,savepath)

e_mae, e_rmse = get_error(data_e,pred_e)
f_mae, f_rmse = get_error(data_f,pred_f)

write_err = open(errfile, 'w')
write_err.write("dir index = %s\n" %frame_idx)
write_err.write("energy mae = %8.2f (eV),           rmse = %8.2f (eV)\n"   %(e_mae, e_rmse))
write_err.write("energy mae = %8.2f (meV/atom),     rmse = %8.2f (meV/atom)\n"   %(e_mae*1000/atom_nums, e_rmse*1000/atom_nums))
write_err.write("force  mae = %8.2f (meV/Angstrom), rmse = %8.2f (meV/Angstrom)\n" %(f_mae*1000, f_rmse*1000))
write_err.close()

plotfig2(data_e,pred_e,data_f,pred_f,savepath,atom_nums)