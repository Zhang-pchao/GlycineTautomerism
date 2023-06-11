#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import plumed
import sys
#sys.path.append(root_path + 'scripts') #dirty way to import our script
#from calcFES import calcFES, calcDeltaF, calcESS

#reuse same plumed kernel, to avoid multiple warnings
PLUMED_KERNEL=plumed.Plumed()
from functools import partial
plumed.read_as_pandas = partial(plumed.read_as_pandas, kernel=PLUMED_KERNEL)

#plot-related stuff
import matplotlib.pyplot as plt
#some other plotting functions
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from IPython.display import clear_output
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

try:
# https://github.com/luigibonati/fessa-color-palette/blob/master/fessa.py
    import fessa
    plt.set_cmap('fessa')
except:
    pass # no big deal

def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)

def plot_colvar(ax, method, cv, fmt, colvars):
    colvar = colvars[method]
    ax.plot(colvar['time'][::10], colvar[cv][::10], fmt, color=colors[method], label=cv)
    ax.set_xlim(colvar['time'].min(), colvar['time'].max())
    ax.set_xlabel('time [ps]')
    if cv == 'n1':
        ax.set_ylabel('CV: SOLVATION')
    elif cv == 'd1':
        ax.set_ylabel('CV: IONDISTANCE')
    else:
        ax.set_ylabel(cv+" (kJ/mol)")
    ax.set_title(labels[method])

def plot_two_cvs(method, cv1, cv2, path, fmt, colvars):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches((13, 6))
    plot_colvar(ax1, method, cv1, fmt, colvars=colvars)
    plot_colvar(ax2, method, cv2, fmt, colvars=colvars)
    #plt.show()
    plt.savefig(os.path.join(path,"%s_%s_%s.png"%(method, cv1, cv2)), dpi=150, bbox_inches='tight')
    
def two_cvs(method, cv1, cv2, path, bins, colvars):
    fig = plt.figure(figsize=(13,6), dpi=150, facecolor='white')
    ax = fig.add_subplot(121)
    colvar = colvars[method]
    ax.scatter(colvar[cv1][::], colvar[cv2][::], s=2, color=colors[method])
    ax.set_xlabel('CV: SOLVATION')
    ax.set_ylabel('CV: IONDISTANCE')
    
    ax = fig.add_subplot(122)
    h = ax.hist2d(colvar[cv1][::], colvar[cv2][::], bins=bins,
                  range=[[min(colvar[cv1]),max(colvar[cv1])],[min(colvar[cv2]),max(colvar[cv2])]])
    fig.colorbar(h[3], shrink=1.0)   
    ax.set_xlabel('CV: SOLVATION')
    ax.set_ylabel('CV: IONDISTANCE')
    
    fig.savefig(os.path.join(path,"%s_%s_%s.png"%(method, cv1, cv2)), dpi=150, bbox_inches='tight')

def two_cvs_time(method, cv1, cv2, path, colvars):
    step1 = 2
    step2 = 5
    mdtime = colvars[method]["time"]
    choose1 = int(len(mdtime)/step1/step2)
    choose2 = int(len(mdtime)/step1)
    fig = plt.figure(figsize=(20,10), dpi=150, facecolor='white')
    for j in range(step1):
        for i in range(step2):
            mdtime = colvars[method]["time"][choose1*i+choose2*j:choose1*(i+1)+choose2*j]
            min_t = min(mdtime)
            max_t = max(mdtime)
            color = [plt.get_cmap("viridis", 100)(int(float(i-min_t)/(max_t-min_t)*100)) for i in mdtime]
            ax = fig.add_subplot(step1,step2,i+1+j*step2)
            ax.set_title("%d ps"%max_t)
            if j == step1-1:
                ax.set_xlabel('CV: SOLVATION')
            if i == 0:
                ax.set_ylabel('CV: IONDISTANCE')
            im = ax.scatter(colvars[method]["n1"][choose1*i+choose2*j:choose1*(i+1)+choose2*j],
                            colvars[method]["d1"][choose1*i+choose2*j:choose1*(i+1)+choose2*j], 
                            c=color, marker='.')    
            fig.colorbar(im, shrink=1.0)
    fig.savefig(os.path.join(path,"%s_%s_%s_%s.png"%(method, cv1, cv2, 'time')), dpi=150, bbox_inches='tight')

def plot_ala2_fes(ax, fes, myextent,title=None, max_fes=None):
    if max_fes is None:
        max_fes = np.amax(fes)
    im = ax.imshow(fes, vmax=max_fes, origin='lower', extent=myextent)
    cb = fig.colorbar(im, ax=ax, label='free energy [kJ/mol]')
    #ax.set_yticks([-2,0,2])
    ax.set_aspect('auto')
    ax.set_xlabel('CV: SOLVATION')
    ax.set_ylabel('CV: IONDISTANCE')
    ax.set_title(title)
    return cb

def load_fes(method, tot, fes_folder):
    filename = fes_folder%method+f'fes-rew.dat'
    X = plumed.read_as_pandas(filename, usecols=[2]).to_numpy()
    fes = []
    for i in range(tot):
        filename = fes_folder%method+f'fes-rew_{1+i}.dat'
        fes_i = plumed.read_as_pandas(filename, usecols=[2]).to_numpy()
        nbins = int(np.sqrt(len(fes_i)))
        fes.append(fes_i.reshape(nbins, nbins))
        print(f'loading FES files... {(i+1)/tot:.0%}  ', end='\r')
    f = open(filename, 'r')
    line = f.readline()
    line = f.readline()
    line = f.readline()
    line = f.readline()
    min_n1 = float(line.split()[3])
    line = f.readline()
    max_n1 = float(line.split()[3])
    line = f.readline()
    line = f.readline()
    line = f.readline()
    min_d1 = float(line.split()[3])
    line = f.readline()
    max_d1 = float(line.split()[3])
    return X.reshape(nbins,nbins), fes, (min_n1,max_n1,min_d1,max_d1)

#if __name__ == '__main__':

methods = ['OPES', 'OPES_explore', 'OPES + OPES_explore']
labels = {'OPES' : 'OPES', 'OPES_explore' : 'OPES_explore', 'OPES + OPES_explore' : 'OPES + OPES_explore'}
colors = {'OPES' : "#E69138", 'OPES_explore' : "#674EA7", 'OPES + OPES_explore' : "#6AA84F"}

###change below
root_path  = '/home/pengchao/cp2k/GlycineH2O/Glycine54H2O/GFN1-xTB/'
folder_    = root_path + '%s/Glycine53H2O_OH_NH3CH2COOH/ProtonCV_7/'
file_index = 'OPES_EXP'
method     = methods[1]
sigma      = [0.01,0.02,0.04]#[0.02,0.04,0.06]
tot        = 5
max_fes    = 150 #kJ/mol
###change above

save_path  = os.path.join(folder_%file_index,'analysis')
mkdir(save_path)
colvars    = {}
colvars[method] = plumed.read_as_pandas(folder_%file_index+'COLVAR')

###change below
#plot_two_cvs(method, 'd1', 'opes.bias', save_path , '.', colvars)
plot_two_cvs(method, 'd1', 'opese.bias',save_path, '.', colvars)
#plot_two_cvs(method, 'n1', 'opes.bias', save_path , '.', colvars)
plot_two_cvs(method, 'n1', 'opese.bias',save_path, '.', colvars)
###change above

two_cvs(method,'n1','d1', save_path,500,colvars)
two_cvs_time(method,'n1','d1', save_path,colvars)

X1, fes_opes1, extent_opes = load_fes(file_index,tot,fes_folder=folder_+'fes_reweight/sigma_%.2f_%.2f/'%(sigma[0],sigma[0]))
X2, fes_opes2, extent_opes = load_fes(file_index,tot,fes_folder=folder_+'fes_reweight/sigma_%.2f_%.2f/'%(sigma[1],sigma[1]))
X3, fes_opes3, extent_opes = load_fes(file_index,tot,fes_folder=folder_+'fes_reweight/sigma_%.2f_%.2f/'%(sigma[2],sigma[2]))

fig, (ax1,ax2,ax3) = plt.subplots(1, 3)
fig.set_size_inches((18, 6))
plot_ala2_fes(ax1, X1, extent_opes, "sigma=%.2f,%.2f"%(sigma[0],sigma[0]), max_fes=max_fes)
plot_ala2_fes(ax2, X2, extent_opes, "sigma=%.2f,%.2f"%(sigma[1],sigma[1]), max_fes=max_fes)
plot_ala2_fes(ax3, X3, extent_opes, "sigma=%.2f,%.2f"%(sigma[2],sigma[2]), max_fes=max_fes)
plt.savefig(os.path.join(save_path,"%s_%s_%s.png"%(method, "fes_reweight", "all")),
            dpi=150, bbox_inches='tight')

for i in range(tot):
    fig, (ax1,ax2,ax3) = plt.subplots(1, 3)
    fig.set_size_inches((18, 6))
    plot_ala2_fes(ax1, fes_opes1[i], extent_opes, "sigma=%.2f,%.2f"%(sigma[0],sigma[0]), max_fes=max_fes)
    plot_ala2_fes(ax2, fes_opes2[i], extent_opes, "sigma=%.2f,%.2f"%(sigma[1],sigma[1]), max_fes=max_fes)
    plot_ala2_fes(ax3, fes_opes2[i], extent_opes, "sigma=%.2f,%.2f"%(sigma[2],sigma[2]), max_fes=max_fes)
    plt.savefig(os.path.join(save_path,"%s_%s_%s.png"%(method, "fes_reweight", str(i+1))),
                dpi=150, bbox_inches='tight')
