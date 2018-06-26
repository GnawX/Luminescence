from scipy.special import comb,eval_hermite,factorial,factorial2
from math import sqrt,exp
import sys

# Franck-Condon factor |<v|v'>|^2
# J Mol. Spectroscopy 232 (2005) 102-104
# parameters a = omega/h_bar d = x2 - x1 (x normal coordinate)
# a1, x1 for ground state, a2, x2 for excite state

# initialization
d = 4.524;
a1,a2 = 0.018, 0.018
v1 = int(sys.argv[1])
v2 = int(sys.argv[2])
# constants
hbar = 6.582119514e-16 # eV.s/rad
ev = 1.6021766208e-19 # J
amu = 1.660539040e-27 # kg
ang = 1.0e-10 # m
fr = ev/amu/ang**2
# convert frequency to 1/sqrt(amu)/ang omega = omega^2/(hbar*omega)
a1 = (a1/hbar)**2/fr/a1
a2 = (a2/hbar)**2/fr/a2

S = a1*a2*d**2/(a1 + a2)
b1 = -a2*sqrt(a1)*d/(a1 + a2)
b2 = a1*sqrt(a2)*d/(a1 + a2)
A = 2*sqrt(a1*a2)/(a1 + a2)

def fcint(v1,v2):
    aa = 0.0
    for k1 in range(v1+1):
        for k2 in range(v2+1):
            if (k1 + k2) % 2 == 0:
                ik = factorial2(k1+k2-1)/(a1+a2)**((k1+k2)/2)
            else:
                ik = 0
            aa += comb(v1,k1)*comb(v2,k2)*eval_hermite(v1-k1,b1)\
                  *eval_hermite(v2-k2,b2)*(2*sqrt(a1))**k1*(2*sqrt(a2))**k2*ik
    bb = aa*A*exp(-S)/2**(v1+v2)/factorial(v1)/factorial(v2)*aa
    #bb = A*exp(-S)/2**(v1+v2)
    return bb

#v2 = 0
#v1,v2 = 0,0
#for v1 in range(15):
#    e = (v1 + 0.5)*a1
print fcint(v1,v2) 
