import numpy as np                                   
from math import sqrt, exp, pi                       
from functions import fcf_decimal, gaussian          

""" calculate the nonradiative recombination rate """
# kn = 2pi/hbar * W^2 *Sum_n,m p_m*|<n|Q|m>|^2*delta(zpl+e_m-e_n)

# input parameters
ni    =    50                   # initial
nf    =    350                   # final 
wi    =    0.009937   # in eV            
wf    =    0.010445                      
q     =    9.04728    # in sqrt(amu)*ang 
sigma =    0.225      # broadening       
t     =    300                   # K     
zpl   =    2.70845    # eV               
wij   =    6e-4   # <0|dPsi/dQ> unit 1/Q 
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
f = fcf_decimal(ni,nf,wi,wf,q)
z = sum(exp(-(0.5 + i)*wi/t) for i in range(200)) # partition
beta = sqrt(wi*fr/(wi/hbar)**2/2)  # hbar/omega unit Q^2
#

rate = 0.0
for i in range(ni-1):
    ei = (0.5 + i)*wi
    bz = exp(-ei/t)/z
    for j in range(nf):
        ef = (0.5 + j)*wf
        if i < 1:
           rate += 2*pi*beta**2*bz*(sqrt(i+1)*f[i+1,j])**2*zpl**2*wij**2* \
                   gaussian(zpl+ei,ef,sigma)
        else:
           rate += 2*pi*beta**2*bz*(sqrt(i+1)*f[i+1,j]+sqrt(i)*f[i-1,j])**2* \
                   zpl**2*wij**2*gaussian(zpl+ei,ef,sigma)

print rate/h  # Hz
