#!/usr/bin/env python
# coding: utf-8

# In[13]:


import numpy as np
import os
import sys
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis


# In[14]:


import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.patches as patches
import matplotlib.patheffects as path_effects
import matplotlib.image as mpimg
from matplotlib.pyplot import MultipleLocator
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from IPython.display import clear_output


# In[15]:


import plumed
PLUMED_KERNEL=plumed.Plumed()
from functools import partial
plumed.read_as_pandas = partial(plumed.read_as_pandas, kernel=PLUMED_KERNEL)


# In[16]:


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


# In[35]:


geo_path = "/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/020_1e-6_add_water/008_gly128H2O_sd_B35_compress_f/0-30ns"
trj_path = '/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/020_1e-6_add_water/008_gly128H2O_sd_B35_compress_f/0-30ns'
geo_name = os.path.join(geo_path,"Glycine128H2O_opt_atomic.data")
trj_name = os.path.join(trj_path,"trj_file_2/anion","anion_00-30ns.lammpstrj")
#trj_name = os.path.join(trj_path,"trj_file_2/anion","cation_00-30ns.lammpstrj")
cvs_file = plumed.read_as_pandas(os.path.join(trj_path,"trj_file_2/anion",'COLVAR'))
ion_type = ['cp','ip'] # ['cm','im']
save_path = os.path.join(trj_path,"IonHBnumber",ion_type[1])


# In[18]:


u = mda.Universe(geo_name,trj_name, atom_style='id type x y z',format='LAMMPSDUMP')


# In[19]:


cvs_step = 1
trj_step = 1
step = int(trj_step/cvs_step)
#print(cvs_file[::step])


# In[20]:


cvs = cvs_file[::step]
tcs = cvs['tc'].to_numpy()
cpm = cvs[ion_type[0]].to_numpy()
ion = cvs[ion_type[1]].to_numpy()
d05 = cvs["d05"].to_numpy()


# In[21]:


print(tcs[:5])
print(ion[:5])


# In[22]:


def isInt(x):
    if x%1 == 0:
        return True
    else:
        return False


# In[23]:


def lmp_index(trj_file):
    searchfile = open(trj_file, "r")
    lines = searchfile.readlines()
    searchfile.close()
    
    xyz_index = []
    i = 0
    for line in lines:
        if "ITEM: ATOMS id type x y z" in line:
            xyz_index.append(i) 
        i += 1
    
    start_index = []
    j = 0
    for line in lines:
        if "ITEM: TIMESTEP" in line:
            start_index.append(j)
        j += 1
            
    return xyz_index,start_index


# In[24]:


def get_hbonds(u,index,scheme):
    if scheme == "donor":
        hbonds = HydrogenBondAnalysis(
            universe           = u,
            donors_sel         = "index {0}".format(index), # O
            hydrogens_sel      = "type 1", # H
            acceptors_sel      = "type 2", # O
            d_a_cutoff         = 3.5,      # <3.5
            d_h_a_angle_cutoff = 140,      # >140
            update_selections  = False
        )
    else:  # scheme == "acceptor"
        hbonds = HydrogenBondAnalysis(
            universe           = u,
            donors_sel         = "type 2", # O
            hydrogens_sel      = "type 1", # H
            acceptors_sel      = "index {0}".format(index), # O
            d_a_cutoff         = 3.5,      # <3.5
            d_h_a_angle_cutoff = 140,      # >140
            update_selections  = False
        )    
    return hbonds


# In[32]:


def get_ion_hb_dict(u,o_id,scheme,d05,frame,hb_dict,name):    
    hbonds = get_hbonds(u,o_id,scheme)
    hbonds.run(start=frame,stop=frame+1,step=None,verbose=False)
    u.trajectory[frame]
    o_d_or_a = hbonds.results.hbonds.shape[0]
    hb_dict[d05] = o_d_or_a
    #print(name,scheme,o_id,d05,o_d_or_a)    
    return hb_dict,o_d_or_a,hbonds.results.hbonds


# In[33]:


def get_neg_hb_dict(u,o_id,scheme,d05,frame,hb_dict,name):    
    hbonds = get_hbonds(u,o_id,scheme)
    hbonds.run(start=frame,stop=frame+1,step=None,verbose=False)
    u.trajectory[frame]
    o_d_or_a = hbonds.results.hbonds.shape[0]
    hb_dict[d05] = o_d_or_a
    #print(name,scheme,o_id,d05,o_d_or_a)    
    return hb_dict


# In[26]:


def write_dict(path,ion_dict,name):
    file = open(os.path.join(path, "IonHBnumber_{0}.txt".format(name)), 'w+')
    file.write('# %16s%16s\n' %('D(Angstrom)','HBnumber'))
    for i in ion_dict.keys():
        file.write('%16.8f%16d\n' %(i,ion_dict[i]))
    file.close()  


# In[36]:


def mkdir(path):
    path=path.strip()
    path=path.rstrip("\\")
    isExists=os.path.exists(path) 

    if not isExists:
        os.makedirs(path) 


# In[45]:


def get_hist2d(hb_dict,edge_start,edge_end,binss,name):
    fig = plt.figure(figsize=(7.9,6.5), dpi=150, facecolor='white')    
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(name)
    h = ax.hist2d(hb_dict.keys(),hb_dict.values(),bins=binss,
                  range=[[edge_start-1, edge_end+1],[-0.75, 7.5]],
                  cmap='Blues')
    fig.colorbar(h[3], shrink=1.0)
    ax.set_xlabel("CV d05 "+r"$\ \rm (\AA)$")
    ax.set_ylabel("Number of adjacent water molecules")
    fig.savefig(os.path.join(save_path, "HBNeg_%s.png"%(name)), dpi=600, bbox_inches='tight')


# In[27]:


ion_idx = []
for i in range(len(tcs)):
#for i in range(15):    
    u.trajectory[i]
    if 0.999<=tcs[i]<= 1.001 and 0.999<=np.abs(cpm[i])<= 1.001:
        if isInt(ion[i]):
            ion_idx.append(int(ion[i]*3+1)) # *3 for O H H, +1 for serial in .data file
        else:
            ion_idx.append(-10)
    else:
        ion_idx.append(-10)


# In[34]:


ion_donor_dict = {} # find ion as donor
ion_accpt_dict = {}
ion_d_a_d_dict = {} # ion as donor, find its acceptor, then find its acceptor's donor
ion_d_a_a_dict = {}
ion_a_d_d_dict = {}
ion_a_d_a_dict = {}
#for i in range(1,len(ion_idx),trj_step):
for i in [1,2,3]: 
    if ion_idx[i] != -10:
        ion_o_id = ion_idx[i]-1 # -1: index in mda
        ion_donor_dict,ion_o_donors,hbonds_results = get_ion_hb_dict(u,ion_o_id,"donor",d05[i],i,ion_donor_dict,"ion_donor_dict")        
        if ion_o_donors > 0:
            for frame, donor_index, hydrogen_index, acceptor_index, DA_distance, DHA_angle in hbonds_results:
                ion_d_a = acceptor_index.astype(int)
                ion_d_a_d_dict = get_neg_hb_dict(u,ion_d_a,"donor",   d05[i],i,ion_d_a_d_dict,"ion_d_a_d_dict")
                ion_d_a_a_dict = get_neg_hb_dict(u,ion_d_a,"acceptor",d05[i],i,ion_d_a_a_dict,"ion_d_a_a_dict")
                        
        ion_accpt_dict,ion_o_accpts,hbonds_results = get_ion_hb_dict(u,ion_o_id,"acceptor",d05[i],i,ion_accpt_dict,"ion_accpt_dict")
        if ion_o_accpts > 0:
            for frame, donor_index, hydrogen_index, acceptor_index, DA_distance, DHA_angle in hbonds_results:
                ion_a_d = donor_index.astype(int)
                ion_a_d_d_dict = get_neg_hb_dict(u,ion_a_d,"donor",   d05[i],i,ion_a_d_d_dict,"ion_a_d_d_dict")
                ion_a_d_a_dict = get_neg_hb_dict(u,ion_a_d,"acceptor",d05[i],i,ion_a_d_a_dict,"ion_a_d_a_dict")


# In[37]:


mkdir(save_path)
write_dict(save_path,ion_donor_dict,"ion_donor_dict")
write_dict(save_path,ion_d_a_d_dict,"ion_d_a_d_dict")
write_dict(save_path,ion_d_a_a_dict,"ion_d_a_a_dict")
write_dict(save_path,ion_accpt_dict,"ion_accpt_dict")
write_dict(save_path,ion_a_d_d_dict,"ion_a_d_d_dict")
write_dict(save_path,ion_a_d_a_dict,"ion_a_d_a_dict")


# In[ ]:


get_hist2d(ion_donor_dict,-0.5,12.5,[150,25],"ion_donor")
get_hist2d(ion_d_a_d_dict,-0.5,12.5,[150,25],"ion_d_a_d")
get_hist2d(ion_d_a_a_dict,-0.5,12.5,[150,25],"ion_d_a_a")
get_hist2d(ion_accpt_dict,-0.5,12.5,[150,25],"ion_accpt")
get_hist2d(ion_a_d_d_dict,-0.5,12.5,[150,25],"ion_a_d_d")
get_hist2d(ion_a_d_a_dict,-0.5,12.5,[150,25],"ion_a_d_a")


# In[ ]:




