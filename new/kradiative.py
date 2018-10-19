import sys
import numpy as np
from math import sqrt, exp, pi
import matplotlib.pyplot as plt
from functions import fcf, gaussian, fcf_decimal

# input parameters
ni    =    50                   # initial
nf    =    300                   # final
wi    =    0.008   # in eV
wf    =    0.0075
k     =    8.71    # in sqrt(amu)*ang
sigma =    0.2665      # broadening
t     =    175                   # K
zpl   =    2.74142    # eV
es    =    1.0
ed    =    4.0
ne    =    1000                   # energy grid
eps   =    3.28
kb    =    8.6173303e-5  # eV/K
eh    =    27.211385 # eV
tau   =    2.418884e-17 # hbar/Eh
# initialize
t = t*kb
f = fcf_decimal(ni,nf,wi,wf,k)
f = f**2
z = sum(exp(-(0.5 + i)*wi/t) for i in range(200)) # partition
e = np.linspace(es,ed,ne)
a = sqrt(eps)/3/pi/137**3/tau
pl = np.zeros(ne)
#
for ie in range(ne):
    for i in range(ni):
        ei = (0.5 + i)*wi
        bz = exp(-ei/t)/z
        for j in range(nf):
            ef = (0.5 + j)*wf
            pl[ie] += a*bz*(e[ie]/eh)**3*f[i,j]*gaussian(e[ie],zpl+ei-ef,sigma)

print sum(pl)*(e[1]-e[0])
