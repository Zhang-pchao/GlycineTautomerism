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
        lines = f.readlines()
        for line in lines:
            if "TOTAL ATOM NUMBER =" in line:
                sp = line.split()[4]
    return int(sp)

def read_cell(file):
    cell = []
    with open(file) as f: 
        lines = f.readlines()[:]
        idx = 0
        for line in lines:
            idx += 1
            if "Lattice vectors:" in line:
                break
        for i in range(3):
            sp = lines[idx+i].split()
            cell.append([float(sp[0]),float(sp[1]),float(sp[2])])
    return cell

def read_xyz(file,typemap,atom_num):
    total = []
    with open(file) as f:
        idx = 0
        lines = f.readlines()
        for line in lines:
            idx += 1
            if "CARTESIAN COORDINATES (" in line:
                break
                
        for i in range(atom_num):
            sp = lines[idx+i+1].split()
            for j in range(3):
                total.append(float(sp[j+1]))
    return total

def read_force(file,typemap,atom_num):
    total = []
    with open(file) as f: 
        lines = f.readlines()[:]
        idx = 0
        for line in lines:
            idx += 1
            if "TOTAL-FORCE (eV/Angstrom)" in line:
                break
        for i in range(atom_num):
            sp = lines[idx+4+i].split()
            for j in range(3):
                total.append(float(sp[j+1]))        
    return total

def read_energy(file):
    total = []
    with open(file) as f: 
        lines = f.readlines()[:]
        for line in lines:
            if "!FINAL_ETOT_IS" in line:
                sp = line.split()
                total.append(float(sp[1]))
    return total

def type_raw(file, typemap, type_path, map_path):
    w_map = open(map_path, 'w')
    for i in typemap:
        w_map.write(i + '\n')
    w_map.close()
    
    type_flag = 0
    type_num  = []
    w_type = open(type_path, 'w')
    with open(file) as f: 
        lines = f.readlines()
        for k in range(len(typemap)):
            for line in lines:
                if "atom label for species %d = %s"%(k+1,typemap[k]) in line:
                    type_flag += 1                    
        if type_flag == len(typemap):
            for line in lines:
                if "number of atom for this type" in line:
                    type_num.append(int(line.split()[7]))
            for i in range(len(typemap)):
                for j in range(type_num[i]):
                    w_type.write(str(i) + '\n')
        else:
            print("check the atom type order!")
    w_type.close()

def check_scf(file):
    total = 0
    with open(file) as f: 
        lines = f.readlines()[:]
        for line in lines:
            if "charge density convergence is achieved" in line:
                total = 1
                break
    return total    

def choose_my_scf(file):
    with open(file) as f: 
        lines = f.readlines()[::-1]
        for line in lines:
            if "Density error is" in line:
                scf_err = float(line.split()[3])
                break
    return scf_err   

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
scf_error = 1e-6
# input data path here, string, this directory should contains
data_path = "/data/HOME_BACKUP/pengchao/glycine/deepks/water_ion/pbe/test/"
typemap   = ["H","O"]
all_dir   = ["002_tst140"]
for i in all_dir:
    dirct     = os.path.join(data_path,i)
    dir_list  = get_dirlist(dirct)
    scf_num   = 0
    my_scf_num= 0
    save_path1 = os.path.join(data_path,  "dp_dataset_scf_1e-6",i)
    save_path2 = os.path.join(save_path1, "set.000")
    mkdir(save_path2)
    type_path  = os.path.join(save_path1,  "type.raw")
    map_path   = os.path.join(save_path1,  "type_map.raw")
    box_path   = os.path.join(save_path2,  "box.npy")
    coord_path = os.path.join(save_path2,  "coord.npy")
    force_path = os.path.join(save_path2,  "force.npy")
    energy_path= os.path.join(save_path2,  "energy.npy")
    
    choose_frames = dir_list[:]
    total_xyz     = []
    total_box     = []
    total_force   = []
    total_energy  = []
    scf_frames    = []

    for j in choose_frames:
        abacus_path = os.path.join(dirct,j)
        abacus_out  = os.path.join(abacus_path,"running_scf.log")
        atom_num = read_atom_num(abacus_out)
        cell = read_cell(abacus_out)

        choose_scf = check_scf(abacus_out)
        my_scf_err = choose_my_scf(abacus_out)
        scf_num   += check_scf(abacus_out)
        #if choose_scf == 1: #err = 1e-7
        if my_scf_err <= scf_error:
            my_scf_num += 1
            total_xyz.append(read_xyz(abacus_out,typemap,atom_num))
            total_force.append(read_force(abacus_out,typemap,atom_num))
            total_energy.append(read_energy(abacus_out))
            total_box.append(cell)
            scf_frames.append(j)
        
    print("%s scf_num [below %1.8f] = %d/%d" %(i,1e-7,scf_num,len(choose_frames))) 
    print("%s scf_num [below %1.8f] = %d/%d" %(i,scf_error,my_scf_num,len(choose_frames))) 
    print("scf_frames = ",scf_frames)
    type_raw(abacus_out, typemap, type_path, map_path)
    
    total_xyz_np = np.asarray(total_xyz,dtype=float).reshape(len(scf_frames),atom_num*3)
    np.save(coord_path, total_xyz_np)
    
    total_force_np = np.asarray(total_force,dtype=float).reshape(len(scf_frames),atom_num*3)
    np.save(force_path, total_force_np)
    
    total_energy_np = np.asarray(total_energy,dtype=float).reshape(len(scf_frames),1)
    np.save(energy_path, total_energy_np)   
    
    total_box_np = np.asarray(total_box,dtype=float).reshape(len(scf_frames),9)
    np.save(box_path, total_box_np)

check_npy(loadData = np.load(os.path.join(save_path2,'coord.npy')))
check_npy(loadData = np.load(os.path.join(save_path2,'box.npy')))
check_npy(loadData = np.load(os.path.join(save_path2,'force.npy')))
check_npy(loadData = np.load(os.path.join(save_path2,'energy.npy')))