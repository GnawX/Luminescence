# energy of excited state E = 1/2 omega^2 Q^2 + lambda Q^3                                                                                                               
import sys                                                                                                                                                               
import numpy as np                                                                                                                                                       
from math import sqrt, exp, pi                                                                                                                                           
import matplotlib.pyplot as plt                                                                                                                                          
from functions import fcf, gaussian, fcf_decimal                                                                                                                         

# input parameters
ni    =    30                   # initial
nf    =    200                   # final 
wi    =    0.009937   # in eV            
wf    =    0.010445                      
k     =    9.04728    # in sqrt(amu)*ang 
lam   =    6.88757e-4 # lambda cubic term coeff unit eV/(A sqrt(amu))^3
sigma =    0.225      # broadening                                     
t     =    100                   # K                                   
zpl   =    2.70845    # eV                                             
es    =    1.0                                                         
ed    =    3.0                                                         
ne    =    400                   # energy grid                         
kb    =    8.6173303e-5  # eV/K                                        
# constants                                                            
kb    =    8.6173303e-5  # eV/K                                        
hbar  =    6.582119514e-16 # eV.s/rad                                  
h     =    4.135667662e-15 # eV.s                                      
ev    =    1.6021766208e-19 # J                                        
amu   =    1.660539040e-27 # kg                                        
ang   =    1.0e-10 # m                                                 
fr    =    ev/amu/ang**2                                               
# initialize                                                           
t = t*kb                                                               
f = fcf_decimal(ni,nf,wi,wf,k)                                         
z = sum(exp(-i*wi/t) for i in range(100)) # partition          
beta = sqrt(wi*fr/(wi/hbar)**2/2)  # sqrt(hbar/2omega) unit Q          
e = np.linspace(es,ed,ne)                                              
pl = np.zeros(ne)
ff = np.zeros((ni,nf))
#
for ie in range(ne):
    for i in range(ni-3):
        ei = i*wi - 30*lam**2*beta**6/wi*(i**2+i+11./30)
        bz = exp(-ei/t)/z
        for j in range(nf):
            ef = (0.5 + j)*wf
            if i < 1:
               ff[i,j] = f[i,j]-lam*beta**3/3/wi*(9*sqrt((i+1)**3)*f[i+1,j] \
                    + sqrt((i+1)*(i+2)*(i+3))*f[i+3,j])
            elif i >= 1 and i < 3:
                 ff[i,j] = f[i,j]-lam*beta**3/3/wi*(9*sqrt(i**3)*f[i-1,j] + \
                      9*sqrt((i+1)**3)*f[i+1,j] + \
                      sqrt((i+1)*(i+2)*(i+3))*f[i+3,j])
            else:
                 ff[i,j] = f[i,j]-lam*beta**3/3/wi*(sqrt(i*(i-1)*(i-2))*f[i-3,j] + \
                      9*sqrt(i**3)*f[i-1,j] + \
                      9*sqrt((i+1)**3)*f[i+1,j] + \
                      sqrt((i+1)*(i+2)*(i+3))*f[i+3,j])
            pl[ie] += bz*e[ie]**3*ff[i,j]**2*gaussian(e[ie],zpl+ei-ef,sigma)

tmp = np.vstack((e,pl))
np.savetxt('pl2.dat',tmp.T)
plt.plot(e,pl)
plt.show()
