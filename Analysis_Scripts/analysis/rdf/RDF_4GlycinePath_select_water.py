#!/usr/bin/env python
# coding: utf-8

import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from matplotlib import cm
import os
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.rdf import InterRDF,InterRDF_s

### Parser stuff ###
parser = argparse.ArgumentParser(description='get the RDF of glycine in water')
# files
parser.add_argument('--trjdir','-td',dest='tdirectory',type=str,default='./',help='the directory of trajectory files')
parser.add_argument('--savdir','-sd',dest='sdirectory',type=str,default='./',help='the directory of saving files')
parser.add_argument('--fstart',dest='fs',type=str,default='anion',help='the name of starting frame')
parser.add_argument('--festop',dest='fe',type=str,default='anion',help='the name of ending   frame')
# some easy parsing
args=parser.parse_args()

#set bigger font sizes
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIG_SIZE = 17
plt.rc('font', size=SMALL_SIZE)        # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)   # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE-4)  # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)   # fontsize of the figure title

####################change below####################
geo_path    = "/data/HOME_BACKUP/pengchao/glycine/dpmd/geofile"
trj_path    = args.tdirectory
save_path   = args.sdirectory
data_geo    = "Glycine128H2O_opt_atomic.data"
trj_skip    = 1
mda_step    = 1
type_map    = {'H': 1,'O': 2,'N': 3,'C': 4,'IonO':9}
NOOdict     = {'N': 385,'O1': 391,'O2': 393} #serial in lammps .data file
Ion_flag    = [args.fs,args.fe]
eles        = [['N','H'],['O','H'],['N','O'],['O','O']]
####################change above####################
waterbox    = {'64': 6.23,'128': 7.84,'256': 9.88,'512': 12.45,'1024':15.69}
data_geo    = os.path.join(geo_path,data_geo)
#trj_file    = os.path.join(trj_path,trj_name)

trjs = []
for filename in os.listdir(trj_path):
    if filename.endswith('.lammpstrj'):
        trjs.append(filename)
print(trjs)

def get_rdf(c1,c2,u,NOOdict,type_map,trj_skip,mda_step,save_path,save_name):
    s1_ = 'index'
    for i in c1:
        s1_ += ' %d'%(NOOdict[i]-1)
    s1 = u.select_atoms(s1_) # lmp serail - 1 = mda index
    s2 = u.select_atoms('type %s'%type_map[c2])
    half_box = u.dimensions[0]/2
    rdf = InterRDF(s1,s2,nbins=int(100*half_box), range=(0.01, half_box))
    rdf.run(start=trj_skip, step=mda_step) # stop=1000
    y_rdf = rdf.rdf
    #y_cdf = rdf.get_cdf()
    x_rdf = rdf.bins
    final_save = os.path.join(save_path,'rdf_data')
    mkdir(final_save)
    file = open(os.path.join(final_save, "RDF_{0}_{1}*_{2}.txt".format(save_name,c1[0][0],c2)), 'w+')
    file.write("#RDF %4s*-%s\n"%(c1[0][0],c2))
    file.write('#%8s%8s\n' %('bin','rdf'))
    for i in range(len(x_rdf)):
        file.write(' %8.4f%8.4f\n' %(x_rdf[i],y_rdf[i]))
    file.close()
    return x_rdf,y_rdf

def get_rdf_N_O(c1,c2,u,NOOdict,type_map,trj_skip,mda_step,save_path,save_name):
    s1_ = 'index'
    for i in c1:
        s1_ += ' %d'%(NOOdict[i]-1)
    s1 = u.select_atoms(s1_) # lmp serail - 1 = mda index
    s2 = u.select_atoms('type %s and not (index %d %d)'%(type_map[c2],NOOdict['O1']-1,NOOdict['O2']-1))
    half_box = u.dimensions[0]/2
    rdf = InterRDF(s1,s2,nbins=int(100*half_box), range=(0.01, half_box))
    rdf.run(start=trj_skip, step=mda_step) # stop=1000
    y_rdf = rdf.rdf
    #y_cdf = rdf.get_cdf()
    x_rdf = rdf.bins
    final_save = os.path.join(save_path,'rdf_data')
    mkdir(final_save)
    file = open(os.path.join(final_save, "RDF_{0}_{1}g_{2}w.txt".format(save_name,c1[0][0],c2)), 'w+')
    file.write("#RDF %4s*-%s\n"%(c1[0][0],c2))
    file.write('#%8s%8s\n' %('bin','rdf'))
    for i in range(len(x_rdf)):
        file.write(' %8.4f%8.4f\n' %(x_rdf[i],y_rdf[i]))
    file.close()
    return x_rdf,y_rdf

def get_rdf_NorO_H(c1,c2,u,NOOdict,type_map,trj_skip,mda_step,save_path,save_name):
    s1_ = 'index'
    for i in c1:
        s1_ += ' %d'%(NOOdict[i]-1)
    s1 = u.select_atoms(s1_) # lmp serail - 1 = mda index
    s2 = u.select_atoms('type %s and not ((around 1.3 type %s) or (around 1.1 type %s) or (around 1.1 index %d %d))'%(
             type_map[c2],type_map['C'],type_map['N'],NOOdict['O1']-1,NOOdict['O2']-1),updating=True)
    half_box = u.dimensions[0]/2    
    rdf = InterRDF(s1,s2,nbins=int(100*half_box), range=(0.01, half_box))
    rdf.run(start=trj_skip, step=mda_step) # stop=1000
    y_rdf = rdf.rdf
    x_rdf = rdf.bins
    final_save = os.path.join(save_path,'rdf_data')
    mkdir(final_save)
    file = open(os.path.join(final_save, "RDF_{0}_{1}g_{2}w.txt".format(save_name,c1[0][0],c2)), 'w+')
    file.write("#RDF %4s*-%s\n"%(c1[0][0],c2))
    file.write('#%8s%8s\n' %('bin','rdf'))
    for i in range(len(x_rdf)):
        file.write(' %8.4f%8.4f\n' %(x_rdf[i],y_rdf[i]))
    file.close()
    return x_rdf,y_rdf

def get_rdf_O_O(c1,c2,u,NOOdict,type_map,trj_skip,mda_step,save_path,save_name):
    half_box = u.dimensions[0]/2
    
    s1 = u.select_atoms('index %d'%int(NOOdict[c1[0]]-1)) # lmp serail - 1 = mda index
    s2 = u.select_atoms('type %s and not (index %d)'%(type_map[c2],NOOdict[c1[1]]-1))
    rdf = InterRDF(s1,s2,nbins=int(100*half_box), range=(0.01, half_box))
    rdf.run(start=trj_skip, step=mda_step) # stop=1000
    y1_rdf = rdf.rdf
    #y_cdf = rdf.get_cdf()
    x_rdf = rdf.bins
    
    s1 = u.select_atoms('index %d'%int(NOOdict[c1[1]]-1)) # lmp serail - 1 = mda index
    s2 = u.select_atoms('type %s and not (index %d)'%(type_map[c2],NOOdict[c1[0]]-1))
    rdf = InterRDF(s1,s2,nbins=int(100*half_box), range=(0.01, half_box))
    rdf.run(start=trj_skip, step=mda_step) # stop=1000
    y2_rdf = rdf.rdf
    #y_cdf = rdf.get_cdf()
    x_rdf = rdf.bins
    
    y_rdf = np.average(np.concatenate((y1_rdf.reshape((1, -1)), y2_rdf.reshape((1, -1))), axis=0), axis=0)
    #print(y1_rdf.shape, y2_rdf.shape, y_rdf.shape)
    final_save = os.path.join(save_path,'rdf_data')
    mkdir(final_save)
    file = open(os.path.join(final_save, "RDF_{0}_{1}g_{2}w.txt".format(save_name,c1[0][0],c2)), 'w+')
    file.write("#RDF %4s*-%s\n"%(c1[0][0],c2))
    file.write('#%8s%8s\n' %('bin','rdf'))
    for i in range(len(x_rdf)):
        file.write(' %8.4f%8.4f\n' %(x_rdf[i],y_rdf[i]))
    file.close()
    return x_rdf,y_rdf

def cmap2colorlist(cmap_name, color_numb):
    colormap = cm.get_cmap(cmap_name, color_numb+4)
    idx = np.arange(color_numb) + 3
    colorlist = colormap(idx)
    return colorlist

def mkdir(path):
    path=path.strip() 
    path=path.rstrip("\\") 
    isExists=os.path.exists(path) 
    if not isExists:
        os.makedirs(path) 

def get_fig(c1,c2,x_rdf,y_rdf,colors,save_path,save_name):
    fig = plt.figure(figsize=(8,3), dpi=150, facecolor='white')
    ax = fig.add_subplot(111)
    ax.grid(True)
    ax.plot(x_rdf, y_rdf, alpha=0.8, lw=3, color=colors[0], label='%2s –%2s'%(c1,c2))
    ax.legend(loc="upper right")
    #ax.set_xlim(0,16)
    ax.set_ylim(0,3)
    ax.set_xlabel("r "+r"$\ \rm (\AA)$")
    ax.set_ylabel("g(r)")
    #fig.show()
    #fig.savefig(os.path.join(save_path, "RDF_{0}_{1}_{2}.png".format(save_name,c1,c2)), dpi=600, bbox_inches='tight')

def get_allfig(eles,path_idx,x_rdf,y_rdf,colors,save_path):
    fig = plt.figure(figsize=(8,4), dpi=150, facecolor='white')
    grid = plt.GridSpec(10, 10, wspace=0.5, hspace=0.5)
    ax2 = fig.add_subplot(grid[0:1,1:9])
    cbar = fig.colorbar(cm.ScalarMappable(cmap=cm.GnBu),cax=ax2,# ticklocation='top',
                       orientation='horizontal', ticks=[0,0.5,1])
    cbar.ax.set_xticklabels([path_idx[0],'transitional',path_idx[1]])
    
    ax = fig.add_subplot(grid[2:10,0:10])
    ax.grid(True)
    for i in range(len(x_rdf)):
        ax.plot(x_rdf[i], y_rdf[i], alpha=0.8, lw=2, color=colors[i])
        
    ax.set_xlim(0,8)
    ax.set_ylim(0,3)
    #ax.set_title("%s* - %sw"%(path_idx[0],path_idx[1]))
    ax.set_xlabel("r "+r"$\ \rm (\AA)$")
    ax.set_ylabel("g(r) of [$\mathregular{%s^G}$ – $\mathregular{%s^W}$]"%(eles[0],eles[1]))
    
    fig.savefig(os.path.join(save_path, "RDF_%sg_%sw_all.png"%(eles[0],eles[1])), dpi=600, bbox_inches='tight')

colors = cmap2colorlist('GnBu', len(trjs))
for j in range(len(eles)):
    X,Y,save_idx = [],[],[]
    if eles[j][0] == 'O':
        ele_list = ['O1','O2']
    else: #N
        ele_list = ['N']
    for i in range(len(trjs)): 
        save_name = trjs[i].split('_')[0]
        trj_file = os.path.join(trj_path,trjs[i])
        u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP') #dt ps,dt=1e-3
        if eles[j][0]   == 'N' and eles[j][1] == 'O':
            x,y = get_rdf_N_O(ele_list,'%s'%eles[j][1],u,NOOdict,type_map,trj_skip,mda_step,save_path,save_name)
        elif eles[j][0] == 'N' and eles[j][1] == 'H':
            x,y = get_rdf_NorO_H(ele_list,'%s'%eles[j][1],u,NOOdict,type_map,trj_skip,mda_step,save_path,save_name)
        elif eles[j][0] == 'O' and eles[j][1] == 'O':
            x,y = get_rdf_O_O(ele_list,'%s'%eles[j][1],u,NOOdict,type_map,trj_skip,mda_step,save_path,save_name)
        else: # 'O','H'
            x,y = get_rdf_NorO_H(ele_list,'%s'%eles[j][1],u,NOOdict,type_map,trj_skip,mda_step,save_path,save_name)    
        get_fig('%s*'%eles[j][0],'%s'%eles[j][1],x,y,colors,save_path,save_name)
        X.append(x)
        Y.append(y)
        save_idx.append(save_name)
    get_allfig(eles[j],Ion_flag,X,Y,colors,save_path)
