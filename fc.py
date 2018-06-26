# Franck-Condon overlap <nf|ni> with n quanta
# P. Ruhoff Chemical Physicas 1994, 186, 355.
# w = omega/h_bar reduced frequency
# q' = q + k normal coordinate shift
# prime for final state
import sys
import numpy as np
from math import sqrt, exp, pi

#input parameters
ni,nf = 10, 100
wi,wf = 0.018, 0.018  # in eV
k = 4.524  # in sqrt(amu)*ang

# constants
hbarev = 6.582119514e-16 # eV.s/rad
hbarj = 1.0545718e-34 # J.s/rad
ev = 1.6021766208e-19 # J
amu = 1.660539040e-27 # kg
ang = 1.0e-10 # m
fr = ev/amu/ang**2

# convert units
wi = (wi/hbarev)**2/fr/wi
wf = (wf/hbarev)**2/fr/wf

# definie parameters
a = (wi - wf)/(wi + wf)
b = 2.*k*sqrt(wi)*wf/(wi + wf)
c = -a
d = -2.*k*sqrt(wf)*wi/(wi + wf)
e = 4.*sqrt(wi*wf)/(wi + wf)

f = np.zeros((ni,nf))
f[0,0] = sqrt(e/2.)*exp(b*d/2./e)
f[1,0] = 1./sqrt(2.)*b*f[0,0]
for i in range(1,nf):
    if i - 2 < 0:
       f[0,i] = 1./sqrt(2.*i)*d*f[0,i-1]
    else:
       f[0,i] = 1./sqrt(2.*i)*d*f[0,i-1] + \
                sqrt((i-1.)/i)*c*f[0,i-2]
for i in range(1,ni):
      for j in range(1,nf):
          if i - 2 < 0:
             f[i,j] = 1./sqrt(2.*i)*b*f[i-1,j] + \
                      0.5*sqrt((j-0.)/i)*e*f[i-1,j-1]
          else:
             f[i,j] = 1./sqrt(2.*i)*b*f[i-1,j] + \
                      sqrt((i-1.)/i)*a*f[i-2,j] + \
                      0.5*sqrt((j-0.)/i)*e*f[i-1,j-1]


