import sys
import numpy as np
from scipy.stats import ks_2samp
from scipy.optimize import curve_fit
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt

###################################################################
'''font = {'family' : 'serif',
        'serif'   : 'palatino',
#      'sans-serif'    : 'Computer Modern Sans serif', 
        'style'   : 'normal',
        'variant'   : 'normal',
        'stretch'   : 'normal',
        'weight'   : 'normal',
        'size'   : 20}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rcParams['figure.figsize'] = (12.0, 8.0)'''

colors=['red','green','purple','blue','cyan']
colors=['C1','C2','C3','C4','C5']
markers=['v','^','<','s','o','*','d','h']


##################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, sharey=True)
##########################
# define theoretical CDF
def func(x,tau):
    return 1-np.exp(-x/tau)
###########################
fdir ='./'
#!cat fpt-B25_1.25.dat fpt-B25_1.3.dat >fpt-B25.dat
fesfiles =['times.dat'] #,'fpt-B25_1.25.dat']

Temps =['300']
#Labels = ['Ala2 BARRIER=30']
out=open('results-ks.dat','w')
for i in range(len(fesfiles)) :
    fesfile  = fdir+fesfiles[i] # realtime in S
    Temperature = Temps[i] # in K

    data = np.loadtxt(fesfile,usecols=(0))
    #unrescaled_data = np.loadtxt(fesfile,usecols=(0))*1E-12 #in ps
    unrescaled_data = np.loadtxt(fesfile,usecols=(0))
    acc_factor = np.sum(data)/np.sum(unrescaled_data)

    ntraj = len(data)

    mint = min(data/100)
    maxt = max(data*100)

    ###########################
    # for numerical stability we divide by the mean
    mu = np.mean(data)
    sigma = np.std(data)
    t_m = np.median(data)

    x=data/mu
    # now compute empirical CDF
    ecdf = ECDF(x)
    #x1 = np.linspace(min(x), max(x),10000)
    x1 = np.logspace(np.log10(mint/mu), np.log10(maxt/mu),10000)
    y1 = ecdf(x1)
    # fit to theoretical CDF to obtain tau
    popt,pcov = curve_fit(func,x1,y1)
    tau=popt[0]
    yfit=func(x1,tau)

    ###########################
    # KS-test
    # generate random data with the same exponential distribution
    x2 = np.random.gamma(1,tau,10000000)
    st,pvalue = ks_2samp(x2,x,mode='exact')
    pvalue = '{:.10f}'.format(pvalue)
    str_tau=str(tau*mu)[:4]
    str_p=str(pvalue)[:4]



    #if i == 0:
    #print('#Temperature ntraj mu tau mu_sem sigma t_m K_value mu_sigma_ratio log2mu_median_ratio tau_mu_ratio pvalue_KS_statistic'+'\n')
    ###########################
    #print('%s %s %s %s %s %s %s %s %s %s %s %s'%(Temperature,ntraj,mu,tau*mu,sigma/np.sqrt(len(data)),
                                                  #sigma,t_m,1/(tau*mu),mu/sigma,np.log10(2)*mu/t_m,
                                                  #tau*mu/mu,pvalue)+'\n')
    print('ntraj = %s \n mu = %s \n tau = %s \n mu_sem = %s \n sigma = %s \n t_m = %s \n K_value = %s \n mu_sigma_ratio = %s \n log2mu_median_ratio = %s \n tau_mu_ratio = %s \n pvalue_KS_statistic = %s \n acceleration_factor = %s \n'%(ntraj,mu,tau*mu,sigma/np.sqrt(len(data)),
                                                  sigma,t_m,1/(tau*mu),mu/sigma,np.log10(2)*mu/t_m,
                                                  tau*mu/mu,pvalue,acc_factor)+'\n')
    #print(tau,tau*mu,pvalue)
    #Plot
    ###########################
    #fig, ax = plt.subplots(1, 1)mi
    ax.step(x1*mu, y1,'k-',lw=1. )
    #ax.plot(x1*mu,yfit,'b-',lw=3.,color=colors[i],label=Labels[i])
    ax.plot(x1*mu,yfit,'b-',lw=3.,color=colors[i])
    ax.set_xscale('log')
    ax.set_xlabel('Time')
    ax.set_ylabel('Cumulative Probability')
    ax.legend(loc='upper left', borderaxespad=0.2,ncol=1,fontsize=18)
  #  ax[0,0].title("CDF: "+str(ntraj)+" Trajectories@"+str(Temperature)+"K P-value="+str_p)

#ax[0].set_xlim([0.0,1.0E+11])
#ax[0].set_ylim([0,1.005])
#ax[1].set_xlim([0.0,1.0E+5])
#ax[1].set_ylim([0,2.005])
plt.xscale('log')
plt.ylim(0,1.005)
plt.title('tau = %0.3e; p = %0.3f'%(float(tau*mu),float(pvalue)),fontsize=16)
plt.tight_layout()
plt.savefig('KS-test.png', bbox_inches='tight', dpi=300,transparent=True)
plt.show()
