import numpy as np

#tail -n +2 COLVAR | awk '{a=a+exp(($5+20)/2.481);b=exp(($5+20)/2.481);c=c+0.002*b;{print $1,$2,c/10^12,($1/10^12)*(a/(NR+1))}}' > time


#================= PARAMETERS ==============================

barrier = 45 #from plumed.dat input file for OPES

T = 300 #Temperature in K

biasColumn = 5 #Column in which the bias is printed in the COLVAR file

dt = 0.001 #The frequency at which colvars is deposited (in ps)

time_unit = 1e12 #in ps

#===========================================================
R = 8.314462618 # J/mol/K
kT = R/1000*T #in kJ/mol

#load COLVAR file and do the analysis
cv = np.loadtxt('COLVAR')

f1 = open('time','w')

a = 0
c = 0
for i in range(len(cv)):
    a += np.exp((cv[i,biasColumn-1]+barrier)/kT)
    b = np.exp((cv[i,biasColumn-1]+barrier)/kT)
    c += dt*b

    print(cv[i,0],cv[i,1],c/time_unit,(cv[i,0]/time_unit)*(a/(i+1)),file=f1)


f1.close()


