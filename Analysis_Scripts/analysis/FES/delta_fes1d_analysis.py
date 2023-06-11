#!/usr/bin/env python
# coding: utf-8

# In[121]:


import numpy as np
import sys
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.patches as patches
import matplotlib.patheffects as path_effects


# In[122]:


import os
import plumed
import argparse

#reuse same plumed kernel, to avoid multiple warnings
PLUMED_KERNEL=plumed.Plumed()
from functools import partial
plumed.read_as_pandas = partial(plumed.read_as_pandas, kernel=PLUMED_KERNEL)

#some other plotting functions
from matplotlib.pyplot import MultipleLocator
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
from IPython.display import clear_output


# In[123]:


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


# In[124]:


colors = [(255,255,255),
#          (31,59,115),
          (34,75,121),
          (40,108,133), 
          (46,141,146),
          (57,156,146),
          (70,168,143),
          (85,180,138),
          (118,194,116),         
          (184,216,81),
          (216,220,72),
          (250,223,63),          
          (255,207,69),
          (255,186,78),          
          (252,164,83),
          (237,133,74),          
#          (222,102,64),
          (214,87,59)
          ]  # R -> G -> B 
          
colornum=16
colors=list(tuple(i/255 for i in color) for color in colors)
colors.reverse()
cm = LinearSegmentedColormap.from_list('my_list', colors, N=colornum)


# In[125]:


path0='/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/020_1e-6_add_water/008_gly128H2O_sd_B35_compress_f'
path1='/0-30ns/fes2D_2/s05_d05_bin500_block3/block3_ani2cat_2'
data1='ani2cat_2.dat.std'
path2='/0-30ns/fes2D_2/s05_d05_bin500_block3/block3_ani2can_2'
data2='ani2can_2.dat.std'
path3='/0-30ns/fes2D_2/s05_d05_bin500_block3/'
data3='deltafes.dat'


# In[126]:


fx1,fy1,fz1,fe1=np.loadtxt(os.path.join(path0+path1,data1),unpack=True)[1:5]
fx1=np.array(fx1)
fy1=np.array(fy1)
fz1=np.array(fz1)
fe1=np.array(fe1)
fi1=[]
for i in range(len(fx1)):
    fi1.append(i)
fi1=np.array(fi1)
fi1=fi1[60:925]-60
fz1=fz1[60:925]
fe1=fe1[60:925]
#print(fi1)


# In[127]:


fx2,fy2,fz2,fe2=np.loadtxt(os.path.join(path0+path2,data2),unpack=True)[1:5]
fx2=np.array(fx2)
fy2=np.array(fy2)
fz2=np.array(fz2)
fe2=np.array(fe2)
fi2=[]
for i in range(len(fx2)):
    fi2.append(i)
fi2=np.array(fi2)
fi2=fi2[60:]-60
fz2=fz2[60:]
fe2=fe2[60:]
#print(fi2)


# In[128]:


def get_min_max(my_arr,flag='min'):
    if flag=='min':
        # Find the index of the minimum value
        min_index = np.argmin(my_arr)
        return min_index
    else:#max        
        # Find the index of the maximum value
        max_index = np.argmax(my_arr)
        return max_index


# In[129]:


PATH1=[]
rang1=[[0,3],[150,200],[400,500],[400,500],[450,600],[700,800],[800,864]]
miax1=['min','max'    ,'min'    ,'max'    ,'min'    ,'max'    ,'min']
name1=['A'  ,'NA'     ,'N'      ,'ZN'     ,'Z'      ,'ZC'     ,'C']
j=0
for i in rang1:
    idx_tmp=get_min_max(fz1[i[0]:i[1]],miax1[j])
    #print(idx_tmp)
    PATH1.append(int(idx_tmp)+i[0])
    j+=1
print(PATH1)


# In[130]:


PATH2=[]
rang2=[[0,3],[50,150],[300,350],[455,460],[450,600],[705,708]]
miax2=['min','max','min','max','max','min']
name2=['A','ZA','Z','ZC','NC','N']
j=0
for i in rang2:
    idx_tmp=get_min_max(fz2[i[0]:i[1]],miax2[j])
    #print(idx_tmp)
    PATH2.append(int(idx_tmp)+i[0])
    j+=1
print(PATH2)


# In[131]:


#fy3,fe3=np.loadtxt(os.path.join(path0+path3,data3),unpack=True)[1:3]
fx3=[0,1,2,3]
fy3=[0.0000,34.2418,44.1679,48.6994]         
fe3=[0.0814,1.4674,2.8366,2.3175]   
fz3=['Z','N','A','C']


# In[132]:


def get_sdn_data(file):
    fx2,fy2,fe2=np.loadtxt(file,unpack=True)[0:3]
    fx2=np.array(fx2)
    fy2=np.array(fy2)
    fe2=np.array(fe2)
    return fx2,fy2,fe2


# In[133]:


#path0='/data/HOME_BACKUP/pengchao/glycine/dpmd/franklin/opes/020_1e-6_add_water/008_gly128H2O_sd_B35_compress_f'
path4='/0-30ns/fes_reweight_block3'
s05x,s05y,s05e=get_sdn_data(os.path.join(path0+path4,'s05/fes-rew.dat'))
s08x,s08y,s08e=get_sdn_data(os.path.join(path0+path4,'s08/fes-rew.dat'))
s100x,s100y,s100e=get_sdn_data(os.path.join(path0+path4,'s100/fes-rew.dat'))
d05x,d05y,d05e=get_sdn_data(os.path.join(path0+path4,'d05/fes-rew.dat'))
d08x,d08y,d08e=get_sdn_data(os.path.join(path0+path4,'d08/fes-rew.dat'))
d100x,d100y,d100e=get_sdn_data(os.path.join(path0+path4,'d100/fes-rew.dat'))
c05x,c05y,c05e=get_sdn_data(os.path.join(path0+path4,'c05/fes-rew.dat'))
c08x,c08y,c08e=get_sdn_data(os.path.join(path0+path4,'c08/fes-rew.dat'))
c100x,c100y,c100e=get_sdn_data(os.path.join(path0+path4,'c100/fes-rew.dat'))


# In[135]:


fig = plt.figure(figsize=(16,8), dpi=150, facecolor='white')
yrange=(-5,75)

ax2 = fig.add_subplot(3,2,2)
ax2.plot(fi1, fz1, lw=2, c=colors[0],label="Minimum free energy path")
ax2.fill_between(fi1,fz1-fe1,fz1+fe1,facecolor=colors[0],alpha = 0.6)
ax2.legend(loc='lower left')
ax2.axes.xaxis.set_ticklabels([])
PATH1tmp=[2, 198, 441, 463, 558, 789, 856]
PATH1tmp[2]-=20
PATH1tmp[3]+=20
#ax2.set_xticks([2, 198, 558, 789, 856],['A','NA','Z','ZC','C'])

# Set x-axis tick labels and adjust the position of some of them
x_tick_labels = name1
x_tick_positions = PATH1tmp
ax2.set_xticks(x_tick_positions)
ax2.set_xticklabels(x_tick_labels)
ax2.tick_params(axis='x', length=0)
    
for i in PATH1:
    ax2.axvline(x=i, color='gray', linestyle='--')
ax2.set_xlim((-60,910))
#ax2.set_xticks([])
#ax2.set_xlabel("Minimum free energy path")
ax2.set_ylim(yrange)
ax2.yaxis.set_major_locator(MultipleLocator(25))
ax2.yaxis.set_minor_locator(MultipleLocator(5))

ax4 = fig.add_subplot(3,2,4)
ax4.plot(fi2, fz2, lw=2, c=colors[-2],label="Another free energy path")#label="%s"%path_file
ax4.fill_between(fi2,fz2-fe2,fz2+fe2,facecolor=colors[-2],alpha = 0.6)
ax4.legend(loc='lower right')
ax4.axes.xaxis.set_ticklabels([])
ax4.set_xticks(PATH2,name2)
for i in PATH2:
    ax4.axvline(x=i, color='gray', linestyle='--')
ax4.tick_params(axis='x', length=0)
#ax4.set_xticks([])
#ax4.set_xlabel("Another free energy path")
ax4.set_ylabel("Free energy (kJ/mol)")
ax4.set_ylim(yrange)
ax4.set_xlim((-50,750))
ax4.yaxis.set_major_locator(MultipleLocator(25))
ax4.yaxis.set_minor_locator(MultipleLocator(5))

ax6 = fig.add_subplot(3,2,6)
ax6.scatter(fx3, fy3, marker='o', color=colors[9],label='liquid phase')
ax6.errorbar(fx3, fy3, yerr=fe3,linestyle="None", capsize=6,color=colors[9]) #
ax6.axes.xaxis.set_ticklabels([])
ax6.legend(loc='lower right')
ax6.set_xlim((-0.5, 3.5))
ax6.set_xticks(fx3,fz3)
ax6.set_ylim((-5,60))
ax6.yaxis.set_major_locator(MultipleLocator(20))
ax6.yaxis.set_minor_locator(MultipleLocator(4))

ax1 = fig.add_subplot(3,2,1)
ax1.plot(s100x, s100y, lw=2, c=colors[-5],label=r'$CV_p$'+r'$^{\lambda=%d}$'%100)
ax1.fill_between(s100x, s100y-s100e,s100y+s100e,facecolor=colors[-5],alpha = 0.6)
ax1.plot(s08x, s08y, lw=2, c=colors[7],label=r'$CV_p$'+r'$^{\lambda=%d}$'%8)
ax1.fill_between(s08x, s08y-s08e,s08y+s08e,facecolor=colors[7],alpha = 0.6)
ax1.plot(s05x, s05y, lw=2, c=colors[2],label=r'$CV_p$'+r'$^{\lambda=%d}$'%5)
ax1.fill_between(s05x, s05y-s05e,s05y+s05e,facecolor=colors[2],alpha = 0.6)
ax1.legend(loc='lower right',frameon=False)
ax1.set_ylim(yrange)
ax1.set_xlim((-1.3,1.3))
ax1.yaxis.set_major_locator(MultipleLocator(25))
ax1.yaxis.set_minor_locator(MultipleLocator(5))
ax1.xaxis.set_minor_locator(MultipleLocator(0.1))

ax3 = fig.add_subplot(3,2,3)
ax3.plot(d100x, d100y, lw=2, c=colors[-5],label=r'$CV_d$'+r'$^{\lambda=%d}$'%100)
ax3.fill_between(d100x, d100y-d100e,d100y+d100e,facecolor=colors[-5],alpha = 0.6)
ax3.plot(d08x, d08y, lw=2, c=colors[7],label=r'$CV_d$'+r'$^{\lambda=%d}$'%8)
ax3.fill_between(d05x, d05y-d05e,d05y+d05e,facecolor=colors[7],alpha = 0.6)
ax3.plot(d05x, d05y, lw=2, c=colors[2],label=r'$CV_d$'+r'$^{\lambda=%d}$'%5)
ax3.fill_between(d05x, d05y-d05e,d05y+d05e,facecolor=colors[2],alpha = 0.6)
ax3.legend(loc='lower right',frameon=False)
ax3.set_ylabel("Free energy (kJ/mol)")
ax3.set_ylim(yrange)
ax3.set_xlim((-1.5,12.5))
ax3.yaxis.set_major_locator(MultipleLocator(25))
ax3.yaxis.set_minor_locator(MultipleLocator(5))
ax3.xaxis.set_minor_locator(MultipleLocator(0.4))

ax5 = fig.add_subplot(3,2,5)
ax5.plot(c100x, c100y, lw=2, c=colors[-5],label=r'$CV_n$'+r'$^{\lambda=%d}$'%100)
ax5.fill_between(c100x, c100y-c100e,c100y+c100e,facecolor=colors[-5],alpha = 0.6)
ax5.plot(c08x, c08y, lw=2, c=colors[7],label=r'$CV_n$'+r'$^{\lambda=%d}$'%8)
ax5.fill_between(c08x, c08y-c08e,c08y+c08e,facecolor=colors[7],alpha = 0.6)
ax5.plot(c05x, c05y, lw=2, c=colors[2],label=r'$CV_n$'+r'$^{\lambda=%d}$'%5)
ax5.fill_between(c05x, c05y-c05e,c05y+c05e,facecolor=colors[2],alpha = 0.6)
ax5.legend(loc='lower right',frameon=False)
ax5.set_ylim((-5,90))
ax5.yaxis.set_major_locator(MultipleLocator(30))
ax5.yaxis.set_minor_locator(MultipleLocator(6))
ax5.xaxis.set_minor_locator(MultipleLocator(0.1))

lax=[ax1,ax2,ax3,ax4,ax5,ax6]
lll=['a','d','b','e','c','f']
for i in range(6):
    lax[i].text(0.01, 0.97, '(%s)'%lll[i], transform=lax[i].transAxes, fontsize=16,
              verticalalignment='top')
save_path=os.path.join(path0,'0-30ns/final_figs')
fig.savefig(os.path.join(save_path,"delta_fes_1d.png"), dpi=600, bbox_inches='tight')


# In[ ]:




