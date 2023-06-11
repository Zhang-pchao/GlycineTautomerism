#!/usr/bin/env python
# coding: utf-8

import pickle
import numpy as np
np.set_printoptions(linewidth=100)
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import os
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from matplotlib.offsetbox import AnchoredText

####################change below####################
geo_path    = "/home/pengchao/glycine/geofile"
trj_path    = "../"
save_path   = "./"
data_geo    = "1015H2O_OH_cation_pbc_atomic_2.data"
trj_name    = "glycine_10.lammpstrj"
trj_skip    = 1
mda_step    = 1
dt_trj      = 1e-3*mda_step*10 #ps
NOOdict     = {'N': 3056,'O1': 3052,'O2': 3053} #serial in lammps .data file
noolist     = ['N','O1','O2']
hb_cutoff   = 0.3
####################change above####################

HBfile = open(os.path.join(save_path, "HB_Time_step{0}_{1}.txt".format(mda_step,trj_name.split('.')[0])), 'w+')
data_geo    = os.path.join(geo_path,data_geo)

def get_avg_std(x):
    avg = np.average(x, axis=0)
    std = np.std(x, axis=0, ddof=1)
    return avg,std

def fit_biexponential(tau_timeseries, ac_timeseries):
    """Fit a biexponential function to a hydrogen bond time autocorrelation function
    Return the two time constants
    """
    from scipy.optimize import curve_fit

    def model(t, A, tau1, B, tau2):
        """Fit data to a biexponential function.
        """
        return A * np.exp(-t / tau1) + B * np.exp(-t / tau2)
    params, params_covariance = curve_fit(model, tau_timeseries, ac_timeseries, [1, 0.5, 1, 2])
    fit_t = np.linspace(tau_timeseries[0], tau_timeseries[-1], 1000)
    fit_ac = model(fit_t, *params)
    return params, fit_t, fit_ac

def set_hbonds(d_idx,a_idx,bond=3.5,angle=140):
    hbonds = HydrogenBondAnalysis(
        universe           = u,
        donors_sel         = d_idx, # O="type 2"
        hydrogens_sel      = "type 1", # H
        acceptors_sel      = a_idx,
        d_a_cutoff         = bond,      # <3.5
        d_h_a_angle_cutoff = angle,      # >140
        update_selections  = False
    )
    return hbonds

def get_hb_num(u,i,trj_info,da='donor',hb_cutoff=0.3,tau_max=100,window_step=1):
    if da == 'donor':
        hbonds = set_hbonds("index %d"%(NOOdict[i]-1),"type 2") # serial-1=index in mda
    else:#accepter
        hbonds = set_hbonds("type 2","index %d"%(NOOdict[i]-1))
    hbonds.run(start=trj_info[1],stop=None,step=trj_info[2],verbose=True)
    counts_hbs = hbonds.count_by_time()
    counts_hbs_avg,counts_hbs_std = get_avg_std(counts_hbs)
    if counts_hbs_avg < hb_cutoff: # no hb
        tau_frames,hbond_lifetime,tau_times = 0,0,0
        nohb = True
    else:
        tau_frames, hbond_lifetime = hbonds.lifetime(tau_max=tau_max, window_step=window_step)
        tau_times = tau_frames * u.trajectory.dt
        nohb = False

    return counts_hbs_avg,counts_hbs_std,tau_frames,tau_times,hbond_lifetime,nohb

def save_data_fit(tau_frames,hbond_lifetime,fit_t,fit_ac,atom,da_flag,save_path,taumax):
    HBfile = open(os.path.join(save_path, "HB_lifetime_%s_%s.txt"%(atom,da_flag)), 'w+')
    HBfile.write('#%8s as %8s, tau_max = %d [ps]\n' %(atom,da_flag,taumax))
    HBfile.write("#%16s%16s%16s%16s\n" %('tau_frames[ps]','hb_lifetime','fit_t[ps]','fit_ac'))
    for i in range(len(tau_frames)):
        HBfile.write(' %16.4f%16.4f%16.4f%16.4f\n' %(tau_frames[i],hbond_lifetime[i],fit_t[i],fit_ac[i]))
    HBfile.close()
    
    fig = plt.figure(figsize=(4,4), dpi=150, facecolor='white')
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(tau_frames, hbond_lifetime, label="data")
    ax.plot(fit_t, fit_ac, label="fit tau_max=%d ps"%taumax)
    ax.set_title(r"HB lifetime of %s as %s"%(atom,da_flag), weight="bold")
    ax.set_xlabel(r"$\tau\ \rm (ps)$")
    ax.set_ylabel(r"$C(\tau)$")
    ax.legend()
    fig.savefig(os.path.join(save_path, "HB_lifetime_%s_%s.png"%(atom,da_flag)), dpi=150, bbox_inches='tight')

trj_file    = os.path.join(trj_path,trj_name)
u = mda.Universe(data_geo,trj_file, atom_style='id type x y z',format='LAMMPSDUMP',dt=dt_trj) # dt[ps]

#print(u.atoms[162].index)
len_trj     = len(u.trajectory)
trj_info = [len_trj,trj_skip,mda_step]

daavg,dastd,datct = [],[],[]
HBfile.write("#%16s%16s%16s%16s%16s\n"%('atom','d_or_a','hb_avg','hb_std','hb_time[ps]'))
for da_flag in ['donor','acceptor']:
    for i in range(len(noolist)):
        taumax = 50
        navg,nstd,tau_frames,tau_times,hbond_lifetime,nohb=get_hb_num(u,noolist[i],trj_info,da_flag,hb_cutoff,taumax)       
        if nohb: # no hb
            time_constant=0
        else:           
            while (hbond_lifetime[-1]>0.03):
                taumax += 50
                if taumax > 900:
                    print("tau_max is not enough for fitting")
                    break
                navg,nstd,tau_frames,tau_times,hbond_lifetime,nohb=get_hb_num(u,noolist[i],trj_info,da_flag,hb_cutoff,taumax)
                    
            params,fit_t,fit_ac = fit_biexponential(tau_frames, hbond_lifetime)
            A, tau1, B, tau2 = params
            time_constant = A * tau1 + B * tau2
            save_data_fit(tau_frames,hbond_lifetime,fit_t,fit_ac,noolist[i],da_flag,save_path,taumax)
        daavg.append(navg)
        dastd.append(nstd)
        datct.append(time_constant)
        HBfile.write(' %16s%16s%16.4f%16.4f%16.4f\n' %(noolist[i],da_flag,navg,nstd,time_constant))
HBfile.close()

colors = [(31 ,59 ,115), 
#          (47 ,146,148), 
#          (80 ,178,141),
#          (167,214,85 ),
#          (255,224,62 ),
#          (255,169,85 ),
          (189, 48,47 )]
colors=list(tuple(i/255 for i in color) for color in colors)
colors.reverse()
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

# HBnumber_Fig
x=[1,2,3]
fig = plt.figure(figsize=(8,8), dpi=150, facecolor='white')

ax = fig.add_subplot(2, 2, 1)
ax.bar(x,daavg[:3],alpha=0.9,label='Donor',color=colors[0],width=0.4)#width=0.3
ax.axes.xaxis.set_ticklabels([])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.set_xticks(x,noolist)
ax.set_ylim((0, 3.5))
ax.set_ylabel("No. of H bonds")
ax.legend(loc='upper right')
at = AnchoredText("(a)", prop=dict(size=12), frameon=True, loc='upper left')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)

ax = fig.add_subplot(2, 2, 2)
ax.bar(x,daavg[3:],alpha=0.9,label='Acceptor',color=colors[0],width=0.4)#width=0.3
ax.axes.xaxis.set_ticklabels([])
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.set_xticks(x,noolist)
ax.set_ylim((0, 3.5))
#ax.set_ylabel("No. of H bonds")
ax.legend(loc='upper right')
at = AnchoredText("(b)", prop=dict(size=12), frameon=True, loc='upper left')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)

ax = fig.add_subplot(2, 2, 3)
ax.bar(x,datct[:3],alpha=0.9,label='Donor',color=colors[1],width=0.4)#width=0.3
ax.axes.xaxis.set_ticklabels([])
ax.xaxis.set_major_locator(MultipleLocator(1))
#ax.yaxis.set_major_locator(MultipleLocator(1))
ax.set_xticks(x,noolist)
#ax.set_ylim((-0.5, 3.5))
ax.set_ylabel("Time constant of H bonds (ps)")
ax.legend(loc='upper right')
at = AnchoredText("(c)", prop=dict(size=12), frameon=True, loc='upper left')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)

ax = fig.add_subplot(2, 2, 4)
ax.bar(x,datct[3:],alpha=0.9,label='Acceptor',color=colors[1],width=0.4)#width=0.3
ax.axes.xaxis.set_ticklabels([])
ax.xaxis.set_major_locator(MultipleLocator(1))
#ax.yaxis.set_major_locator(MultipleLocator(1))
ax.set_xticks(x,noolist)
#ax.set_ylim((-0.5, 3.5))
#ax.set_ylabel("Time constant (ps)")
ax.legend(loc='upper right')
at = AnchoredText("(d)", prop=dict(size=12), frameon=True, loc='upper left')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)

fig.savefig(os.path.join(save_path, "HB_Time_step{0}_{1}.png".format(mda_step,trj_name.split('.')[0])), dpi=600, bbox_inches='tight')
