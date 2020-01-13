import numpy as np
from math import sqrt, exp, pi
import matplotlib.pyplot as plt

# input
t = 300 # K
sigma = 0.01 # eV
mu = 2.9216 # eV ZPL
q = 11.912 # STE displacement
emin = 1 # eV
emax = 3
xmax = 25
A1,A2,A3 = 0.00651, 1.0835e-4, 0  # ES
B1,B2,B3 = 0.00577, 0, 0        # GS

ne = 100
nx = 1000

# constant
kb = 8.617333e-5 # eV/K

egrid = np.linspace(emin,emax,ne)
xgrid = np.linspace(0,xmax,nx)

def boltzmann(e,t):
    return exp(-e/kb/t)
def gaussian(x,sigma):
    return 1/sigma/sqrt(pi*2)*exp(-(x)**2/2/sigma**2)
def lorentzian(x,sigma):
    return 1/pi*sigma/2/(x**2 + sigma**2/4)
def es(x):
    return A1*(x-q)**2+A2*(x-q)**3+A3*(x-q)**4+mu
def gs(x):
    return B1*x**2+B2*x**3+B3*x**4

pl = np.zeros(ne)
for i in range(ne):
    tmp = 0.0
    for j in range(nx):
        tmp += boltzmann(es(xgrid[j])-mu,t)*gaussian(es(xgrid[j])-gs(xgrid[j])-egrid[i],sigma)
    pl[i] = tmp*(xgrid[1]-xgrid[0])

tmp2 = np.vstack((egrid,pl))
np.savetxt('pl.dat',tmp2.T)
#plt.plot(egrid,pl)
#plt.show()
