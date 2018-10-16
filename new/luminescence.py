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
sigma =    0.225      # broadening
t     =    100                   # K
zpl   =    2.70845    # eV
es    =    1.0
ed    =    3.0
ne    =    400                   # energy grid
kb    =    8.6173303e-5  # eV/K
# initialize
t = t*kb
f = fcf_decimal(ni,nf,wi,wf,k)
f = f**2
z = sum(exp(-(0.5 + i)*wi/t) for i in range(100)) # partition
e = np.linspace(es,ed,ne)
pl = np.zeros(ne)
#
for ie in range(ne):
    for i in range(ni):
        ei = (0.5 + i)*wi
        bz = exp(-ei/t)/z
        for j in range(nf):
            ef = (0.5 + j)*wf
            pl[ie] += bz*e[ie]**3*f[i,j]*gaussian(e[ie],zpl+ei-ef,sigma)

tmp = np.vstack((e,pl))
np.savetxt('pl.dat',tmp.T)
plt.plot(e,pl)
plt.show()
