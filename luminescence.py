import sys
import numpy as np
from math import sqrt, exp, pi
from fc import fci

#input parameters
ni,nf = 100, 200
wi,wf = 0.0178, 0.0178  # in eV
k = 4.524  # in sqrt(amu)*ang
sigma = 0.1 # eV
T = 0.025 # eV
zpl = 2.66 # eV

def gaussian(x,u,sigma):
    return 1./sigma/sqrt(2.*pi)*exp(-0.5*(x-u)**2/sigma**2)

def thermal(E,T):
    return exp(-E/T)

f = fci(ni,nf,wi,wf,k)

# Frank-Condon
for en in np.linspace(1.25,3,101,endpoint=True):
    pl = 0.0
    for i in range(ni):
        for j in range(nf):
            es = (i+0.5)*wi
            eg = (j+0.5)*wf
            pl += thermal(es,T)*f[i,j]**2*gaussian(en,zpl+es-eg,sigma)
    print en, pl

# Herzberg-Teller  HT_i
#for en in np.linspace(1.25,3,101,endpoint=True):
#    pl = 0.0
#    for i in range(ni-1):
#        for j in range(nf):
#            es = (i+0.5)*wi
#            eg = (j+0.5)*wf
#            if i - 1 < 0:
#               ht = (sqrt(i)*f[0,j] + sqrt(i+1.)*f[i+1,j])**2
#            else:
#               ht = (sqrt(i)*f[i-1,j] + sqrt(i+1.)*f[i+1,j])**2
#            pl += thermal(es,T)*ht*gaussian(en,zpl+es-eg,sigma)
#    print en, pl
# Herzberg-Teller  HT_f
#for en in np.linspace(1.25,3,101,endpoint=True):
#    pl = 0.0
#    for i in range(ni):
#        for j in range(nf-1):
#            es = (i+0.5)*wi
#            eg = (j+0.5)*wf
#            if j - 1 < 0:
#               ht = (sqrt(j)*f[i,0] + sqrt(j+1.)*f[i,j+1])**2
#            else:
#               ht = (sqrt(j)*f[i,j-1] + sqrt(j+1.)*f[i,j+1])**2
#            pl += thermal(es,T)*ht*gaussian(en,zpl+es-eg,sigma)
#    print en, pl
