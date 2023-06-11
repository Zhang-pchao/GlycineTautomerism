#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os, sys
import glob
import shutil

def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path) 

def get_dirlist(dirct):
    dirList = []
    files=os.listdir(dirct)
    for f in files:
        if os.path.isdir(dirct + '/' + f):
            dirList.append(f)
    return dirList

def read_atom_num(file):
    with open(file) as f: 
        lines = f.readlines()[0]
        sp = lines.split()[0]
    return int(sp)

def read_xyz_old(file,typedict,atom_num):
    total = []
    with open(file) as f: 
        lines = f.readlines()[2:]
        for i in range(atom_num):
            sp = lines[i].split()
            for k in typedict.keys():
                if sp[0] == k:
                    total.append(float(typedict[k]))
            for j in range(3):
                total.append(float(sp[j+1]))
    return total

def read_xyz(file,typemap,typedict,atom_num):
    total = []
    with open(file) as f: 
        lines = f.readlines()[2:]
        for k in typemap:
            for i in range(atom_num):
                sp = lines[i].split()
                if sp[0] == k:
                    total.append(float(typedict[k]))
                    for j in range(3):
                        total.append(float(sp[j+1]))
    return total

def read_force_old(file,atom_num):
    total = []
    with open(file) as f: 
        lines = f.readlines()[:]
        idx = 0
        for line in lines:
            idx += 1
            if "ATOMIC FORCES in" in line:
                break
        for i in range(atom_num):
            sp = lines[idx+2+i].split()
            for j in range(3):
                total.append(float(sp[j+3]))
    return total

def read_force(file,typemap,typedict,atom_num):
    total = []
    with open(file) as f: 
        lines = f.readlines()[:]
        idx = 0
        for line in lines:
            idx += 1
            if "ATOMIC FORCES in" in line:
                break
        for k in typemap:
            for i in range(atom_num):
                sp = lines[idx+2+i].split()
                if sp[2] == k:
                    for j in range(3):
                        total.append(float(sp[j+3]))        
    return total

def read_energy(file):
    total = []
    with open(file) as f: 
        lines = f.readlines()[::-1] # reverse list for reading final scf results 
        for line in lines:
            if "Total energy:" in line:
                sp = line.split()
                total.append(float(sp[2]))
                break
    return total

def check_scf(file):
    total = 0
    with open(file) as f: 
        lines = f.readlines()[:]
        for line in lines:
            if "SCF run converged in" in line:
                total = 1
                break
    return total    

def check_npy(loadData):
    print(type(loadData))
    print(loadData.dtype)
    print(loadData.ndim)
    print(loadData.shape)
    print(loadData.size)
    print(loadData)    

# conversion unit here, modify if you need
au2eV = 27.211386245988
au2A = 0.529177249
# input data path here, string, this directory should contains
data_path = "/data/HOME_BACKUP/pengchao/glycine/deepks/M062X_Dataset/58H2O_ion"
typemap   = ["H","O"]
typedict  = {"H": 1.0, "O": 8.0}
all_dir   = ["001_tst15","001_trn75"]
for i in all_dir:
    scf_num   = 0
    dirct     = os.path.join(data_path,i)
    dir_list  = get_dirlist(dirct)
    
    save_path1 = os.path.join(data_path,  "ks_dataset",i)
    mkdir(save_path1)
    coord_path = os.path.join(save_path1,  "atom.npy")
    force_path = os.path.join(save_path1,  "force.npy")
    energy_path= os.path.join(save_path1,  "energy.npy")
    
    choose_frames = dir_list[:]
    total_xyz     = []
    total_force   = []
    total_energy  = []

    for j in choose_frames:
        #print("system %s frame %s"%(i,j))
        cp2k_path = os.path.join(dirct,j)
        cp2k_xyz  = os.path.join(cp2k_path,"POSCAR1.xyz")
        cp2k_inp  = os.path.join(cp2k_path,"m062x/m062x.inp")
        cp2k_out  = os.path.join(cp2k_path,"m062x/m062x.out")
        atom_num = read_atom_num(cp2k_xyz)
        #cell = read_cell(cp2k_inp)
        total_xyz.append(read_xyz(cp2k_xyz,typemap,typedict,atom_num))
        total_force.append(read_force(cp2k_out,typemap,typedict,atom_num))
        total_energy.append(read_energy(cp2k_out))
        scf_num += check_scf(cp2k_out)
    
    print("%s scf_num = %d/%d" %(i,scf_num,len(choose_frames)))
    total_xyz_np = np.asarray(total_xyz,dtype=float).reshape(len(choose_frames),atom_num,4)
    np.save(coord_path, total_xyz_np)
    
    total_force_np = np.asarray(total_force,dtype=float).reshape(len(choose_frames),atom_num,3)
    np.save(force_path, total_force_np)
    
    total_energy_np = np.asarray(total_energy,dtype=float).reshape(len(choose_frames),1)
    np.save(energy_path, total_energy_np)   

check_npy(loadData = np.load(os.path.join(save_path1,'atom.npy')))
check_npy(loadData = np.load(os.path.join(save_path1,'energy.npy')))
check_npy(loadData = np.load(os.path.join(save_path1,'force.npy')))
