#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import ase
import scipy
from scipy import constants
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_bonds

root_path = '/data/HOME_BACKUP/pengchao/glycine/cp2k/Wannier/026_Gly54H2O_sd_B35_cmps_c2.4_g/ani2can_2'

# units
e = constants.e # C
# conversion factor 1 Debye = 3.34e-30 Cm
conv = 1 / (3.34e-30 * 1e10)
# lambda parameter of the switching function
l = 500
# radial cutoff
lmax = 0.7
# edge of the box and create cubic box
L = 12.028

fixed_atoms = {'GlyN':'162',
               'GlyO': ['170','168'],
               'GlyCn':'164', # C colse to N
               'GlyCo':'165', # C colse to O
               'GlyHc':['166','167'] # H close to C
              }            

path_idxs = [] # glycine tautomerization free energy path
for file in os.listdir(root_path):
    if os.path.isdir(os.path.join(root_path, file)):
        path_idxs.append(file)

def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)

def get_distance(u,atom_indices,L=12.028,pbc=True):
    atom_positions = u.atoms.positions[atom_indices]
    u.dimensions = np.array([L, L, L, 90, 90, 90])
    diff = atom_positions[0] - atom_positions[1]
    if pbc:
        diff = diff - np.round(diff / u.dimensions[:3]) * u.dimensions[:3]
    bond = np.sqrt(np.sum(diff ** 2))
    return diff,bond

def get_voronoi(hatoms,oatoms,u,dim=5):
    OHcouple = np.zeros([len(oatoms),dim])*np.nan
    for j in range(len(hatoms)):
        sum = 0
        w = np.zeros(len(oatoms))
        for i in range(len(oatoms)):
            d,dmod = get_distance(u,[int(oatoms[i].index),int(hatoms[j].index)])
            sum = sum + np.exp(-l*dmod)
        for i in range(len(oatoms)):
            d,dmod = get_distance(u,[int(oatoms[i].index),int(hatoms[j].index)])
            w[i] = np.exp(-l*dmod) / sum
        if np.size(np.where(w>0.8)) == 0: # voronoi cv failed
            break
        else:
            ohcouple = np.where(w>0.8)[0][0]
            for i in range(dim):
                if np.isnan(OHcouple[ohcouple][i]):
                    OHcouple[ohcouple][i] = j
                    break
    if np.size(np.where(w>0.8)) == 0: # voronoi cv failed
        return np.zeros([3,3])*np.nan,False
    else:
        return OHcouple,True

def get_coord(OHcouple_idx,hatoms_xyz):
    if not np.isnan(OHcouple_idx):
        r = hatoms_xyz[int(OHcouple_idx)]
    else:
        r = np.zeros(3)
    return r

def get_mu(dipole):
    mu = np.sqrt(dipole[0]**2 + dipole[1]**2 + dipole[2]**2)
    return mu

def get_avg_std(x):
    avg = np.average(x, axis=0)
    std = np.std(x, axis=0, ddof=1)
    return avg,std

def get_gly_mu(o_n_c_atoms,dipoles):
    gly_dipole = np.zeros(3)
    all_dipole = np.zeros(3)
    rem_dipole = np.zeros(3)
    gln_dipole = np.zeros(3)
    glo_dipole = np.zeros(3)
    for i in range(len(o_n_c_atoms)):
        all_dipole += dipoles[i]
        if o_n_c_atoms[i].index in [162,164,165,168,170]: # index N C C O O of glcine
            gly_dipole += dipoles[i]
        else:
            rem_dipole += dipoles[i]  
            
        if o_n_c_atoms[i].index in [162,164]: # index N C of glcine
            gln_dipole += dipoles[i]
        elif o_n_c_atoms[i].index in [165,168,170]: # index C O of glcine
            glo_dipole += dipoles[i]        
    
    all_mu = get_mu(all_dipole)        
    gly_mu = get_mu(gly_dipole)
    rem_mu = get_mu(rem_dipole)
    gln_mu = get_mu(gln_dipole)
    glo_mu = get_mu(glo_dipole)
    return [gly_mu,gln_mu,glo_mu,rem_mu,all_mu]

def get_gly_solvation_mu(o_n_c_atoms,dipoles,cutoff):
    solvation_list = []
    solvation_list2 = []
    gly_dipole = np.zeros(3)
    for i in range(len(o_n_c_atoms)):
        for j in [162,164,165,168,170]:
            d,dmod = get_distance(u,[int(o_n_c_atoms[i].index),j],pbc=True)
            if dmod < cutoff:
                solvation_list.append(int(o_n_c_atoms[i].index))
    [solvation_list2.append(i) for i in solvation_list if i not in solvation_list2]
    for i in range(len(o_n_c_atoms)):
        if o_n_c_atoms[i].index in solvation_list2:
            gly_dipole += dipoles[i]
    gly_mu = get_mu(gly_dipole)
    return gly_mu

def get_dipole(hatoms,watoms,oatoms,OHcouple,OWcouple,u,L,e,conv,modify_xyz=True,dim=3):
    hatoms_xyz = np.zeros([len(hatoms),dim])*np.nan
    watoms_xyz = np.zeros([len(watoms),dim])*np.nan
    for x in range(np.shape(OHcouple)[0]):
        for y in range(np.shape(OHcouple)[1]):
            if not np.isnan(OHcouple[x][y]):
                at = int(OHcouple[x][y])
            else:
                break
            if modify_xyz:
                d,dmod = get_distance(u,[int(hatoms[at].index),int(oatoms[x].index)],pbc=False)
                position = u.atoms.positions[int(hatoms[at].index)]
                for i in range(3):
                    if d[i] < -L*0.5:
                        if x > np.shape(OHcouple)[0]-5:
                            print('modify_xyz1')
                        hatoms_xyz[at][i] = position[i] + L
                    elif d[i] > L*0.5:
                        if x > np.shape(OHcouple)[0]-5:
                            print('modify_xyz2')
                        hatoms_xyz[at][i] = position[i] - L
                    else:
                        hatoms_xyz[at][i] = position[i]
            else:
                position = u.atoms.positions[int(hatoms[at].index)]
                for i in range(3):
                    hatoms_xyz[at][i] = position[i]
                    
        for z in range(np.shape(OWcouple)[1]):
            if not np.isnan(OWcouple[x][z]):
                atw = int(OWcouple[x][z])
            else:
                break
            if modify_xyz:    
                m,dmod = get_distance(u,[int(watoms[atw].index),int(oatoms[x].index)],pbc=False)
                position = u.atoms.positions[int(watoms[atw].index)]
                for j in range(3):
                    if m[j] < -L*0.5:
                        if x > np.shape(OHcouple)[0]-5:
                            print('modify_xyz3')
                        watoms_xyz[atw][j] = position[j] + L
                    elif m[j] > L*0.5:
                        if x > np.shape(OHcouple)[0]-5:
                            print('modify_xyz4')
                        watoms_xyz[atw][j] = position[j] - L
                    else:
                        watoms_xyz[atw][j] = position[j]
            else:
                position = u.atoms.positions[int(watoms[atw].index)]
                for j in range(3):
                    watoms_xyz[atw][j] = position[j]   
    
    mu = np.zeros(len(oatoms))
    dipoles = np.zeros([len(oatoms),dim])
    for x in range(np.shape(OHcouple)[0]):
        rW = np.zeros(3)
        for k in range(np.shape(OWcouple)[1]):
            rW += get_coord(OWcouple[x][k],watoms_xyz)
        rH = np.zeros(3)
        for k in range(np.shape(OHcouple)[1]):
            rH += get_coord(OHcouple[x][k],hatoms_xyz)
        
        o_type = oatoms[x].type
        if   o_type == 'O':
            rO = u.atoms.positions[int(oatoms[x].index)]*6
        elif o_type == 'N':
            rO = u.atoms.positions[int(oatoms[x].index)]*5
        elif o_type == 'C':
            rO = u.atoms.positions[int(oatoms[x].index)]*4        
            
        dipoles[x] = e*(rO + rH - 2*rW)*conv
        mu[x] = get_mu(dipoles[x])    
    return mu,dipoles,hatoms_xyz,watoms_xyz

def get_wannier_distance(noo,OWcouple,watoms,oatoms,u):
    wannier_d = np.zeros([len(noo),np.shape(OWcouple)[1]])*np.nan
    for i in range(len(noo)):
        for j in range(np.shape(OWcouple)[1]):
            if not np.isnan(OWcouple[noo[i]][j]):
                at = int(OWcouple[noo[i]][j])
                d,dmod = get_distance(u,[int(watoms[at].index),int(oatoms[noo[i]].index)],pbc=True)
                wannier_d[i][j]=dmod   
    return wannier_d

def write_wannier_distance(wannier_d,noo_name,path):
    for i in range(len(noo_name)):
        file = open(os.path.join(path, "%s.txt"%(noo_name[i])), 'w+')
        file.write('#%12s\n'%'distance')
        ele_num=np.shape(wannier_d)[1]
        for j in range(int(np.shape(wannier_d)[0]/ele_num)):
            for k in range(np.shape(wannier_d)[1]):
                if not np.isnan(wannier_d[ele_num*j+i][k]):
                    file.write(' %12.6f\n'%wannier_d[ele_num*j+i][k])
        file.close()

def shift_coord(u,shift_center,L):
    centered_coords = u.atoms.positions - shift_center
    centered_coords = np.where(centered_coords > L/2, centered_coords-L, centered_coords)
    centered_coords = np.where(centered_coords < -L/2, centered_coords+L, centered_coords)
    return centered_coords

u = mda.Universe("%s/%s/%s/m062x/%s"%(root_path,'10001','10001','M062X-HOMO_centers_s1-1_0.xyz'))
h_atoms = u.select_atoms("type H and not (index %s or index %s)"%(fixed_atoms['GlyHc'][0],fixed_atoms['GlyHc'][1]))
all_h_atoms = u.select_atoms("type H")
o_n_atoms = u.select_atoms("type O or type N")
o_n_c_atoms = u.select_atoms("type O or type N or type C")
x_atoms = u.select_atoms("not (type O or type N or type C or type H)")

EleFile = open(os.path.join(root_path, "Element_dipole.txt"), 'w+')
SolFile = open(os.path.join(root_path, "Solvate_dipole.txt"), 'w+')
nco_list = [162,164,165,168,170]
EleFile.write("#%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n"%(
'index','N_%d_a'%nco_list[0],'N_%d_s'%nco_list[0],
        'C_%d_a'%nco_list[1],'C_%d_s'%nco_list[1],
        'C_%d_a'%nco_list[2],'C_%d_s'%nco_list[2],
        'O_%s_a'%"withH"    ,'O_%s_s'%"withH",
        'O_%s_a'%"noH"      ,'O_%s_s'%"noH"))
cut_list = [1.6,2.0,2.4,2.8,3.2,3.6,4.0]
SolFile.write("#%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n"%(
'index','cut%.1f_a'%cut_list[0],'cut%.1f_s'%cut_list[0],
        'cut%.1f_a'%cut_list[1],'cut%.1f_s'%cut_list[1],
        'cut%.1f_a'%cut_list[2],'cut%.1f_s'%cut_list[2],
        'cut%.1f_a'%cut_list[3],'cut%.1f_s'%cut_list[3],
        'cut%.1f_a'%cut_list[4],'cut%.1f_s'%cut_list[4],
        'cut%.1f_a'%cut_list[5],'cut%.1f_s'%cut_list[5],
        'cut%.1f_a'%cut_list[6],'cut%.1f_s'%cut_list[6]))

noo = [-5,-4,-3,-2,-1]
noo_name = ['N','C_n','N_o','O_with_h','O_no_h']

for i in range(len(path_idxs)):
    fram_idxs = [] # frames for wannier center in one directory
    for file in os.listdir(os.path.join(root_path, path_idxs[i])):
        if os.path.isdir(os.path.join(root_path, path_idxs[i],file)):
            fram_idxs.append(file)

    flag = 0 # for sum the number of successful voronoi cv
    for j in range(len(fram_idxs)):
        trj_path = "%s/%s/%s/m062x/%s"%(root_path,path_idxs[i],fram_idxs[j],'M062X-HOMO_centers_s1-1_0.xyz')
        u = mda.Universe(trj_path)
        u.dimensions = np.array([L, L, L, 90, 90, 90])
        shift_center = u.atoms[164].position # shift C_164 as periodic box center
        u.atoms.positions = shift_coord(u,shift_center,L)
        ONCHcouple,HFlag = get_voronoi(all_h_atoms,o_n_c_atoms,u)
        ONCWcouple,WFlag = get_voronoi(x_atoms,o_n_c_atoms,u)
        if HFlag and WFlag: # if False means that skip this loop
            mu_ele,dp_ele,hatoms_xyz,watoms_xyz = get_dipole(all_h_atoms,x_atoms,o_n_c_atoms,
                                                             ONCHcouple,ONCWcouple,u,L,e,conv)

            mu_sol = []
            for k in range(len(cut_list)):
                mu_sol_tmp = get_gly_solvation_mu(o_n_c_atoms,dp_ele,cutoff=cut_list[k])
                mu_sol.append(mu_sol_tmp)
            mu_sol   = np.array(mu_sol)
            
            wannier_d = get_wannier_distance(noo,ONCWcouple,x_atoms,o_n_c_atoms,u)
            if mu_ele[-2] > mu_ele[-1]: # order: O with H, O without H
                tmp1mu = mu_ele[-1]
                tmp2mu = mu_ele[-2]
                mu_ele[-1] = tmp2mu
                mu_ele[-2] = tmp1mu
                wannier_d = wannier_d[[0,1,2,4,3], :] # only for N C C O O five atoms
                
            if flag == 0:
                mu_ele_all = mu_ele[-5:]
                mu_sol_all = mu_sol
                wannier_d_all = wannier_d
            else:    
                mu_ele_all = np.vstack((mu_ele_all, mu_ele[-5:]))
                mu_sol_all = np.vstack((mu_sol_all, mu_sol))
                wannier_d_all = np.vstack((wannier_d_all, wannier_d))
            flag += 1
            
    mu_ele_all_avg,mu_ele_all_std = get_avg_std(mu_ele_all)
    mu_sol_all_avg,mu_sol_all_std = get_avg_std(mu_sol_all) 
    save_path = os.path.join(root_path,'wannier_distance','%d'%(i+10001))
    mkdir(save_path)
    write_wannier_distance(wannier_d_all,noo_name,save_path)
    
    EleFile.write(' %12d'%(10001+i))
    for k in range(len(mu_ele_all_avg)):
        EleFile.write('%12.6f'%mu_ele_all_avg[k])
        EleFile.write('%12.6f'%mu_ele_all_std[k])
    EleFile.write('\n')
    
    SolFile.write(' %12d'%(10001+i))
    for k in range(len(mu_sol_all_avg)):
        SolFile.write('%12.6f'%mu_sol_all_avg[k])
        SolFile.write('%12.6f'%mu_sol_all_std[k])
    SolFile.write('\n')  
    
EleFile.close()
SolFile.close()
