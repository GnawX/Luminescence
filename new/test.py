import numpy as np

# define parameters
# E1 = 1/2 w1^2 Q^2
# E2 = 1/2 w2^2 (Q-d)^2 - dE
w1_   =    0.01045   # eV
w2_   =    0.01061
d     =    9.0056    # sqrt(amu)*A
dE    =    0.4219    # eV
# conversion units
kb    =    8.6173303e-5  # eV/K          
hbar  =    6.582119514e-16 # eV.s/rad
h     =    4.135667662e-15 # eV.s
ev    =    1.6021766208e-19 # J
amu   =    1.660539040e-27 # kg
ang   =    1.0e-10 # m
ha2ev    =    27.211385 # eV
amu2au   =    1822.888486192 # me
a2b   =    1.88973 # ang to bohr
fr    =    np.sqrt(ev/amu/ang**2)*hbar # sqrt(ev/amu/ang**2) to eV
w1    =    w1_/fr
w2    =    w2_/fr
# location of the barrier
b     =    (w2**2*d - np.sqrt(w1**2*w2**2*d**2 + 2*dE*(w2**2 - w1**2)))/(w2**2 - w1**2)
v     =    0.5*w1**2*b**2
print b,v

#for n in range(200):
n = 2
e = (n + 0.5)*w1_
a = np.sqrt(2*e/w1**2)
c = d - np.sqrt(2*(e + dE)/w2**2)
l = 0.5*b*np.sqrt(w1**2*b**2-2*e) - e/w1*np.log(w1*np.sqrt(w1**2*b**2-2*e)+w1**2*b) - \
    0.5*a*np.sqrt(w1**2*a**2-2*e) + e/w1*np.log(w1*np.sqrt(w1**2*a**2-2*e)+w1**2*a) 

r = 0.5*(c-d)*np.sqrt(w2**2*(c-d)**2-2*(dE+e)) - \
    0.5*(b-d)*np.sqrt(w2**2*(b-d)**2-2*(dE+e)) + \
    (dE+e)/w2*np.log((np.sqrt(w2**2*(b-d)**2-2*(dE+e))+w2*(b-d))/(np.sqrt(w2**2*(c-d)**2-2*(dE+e))+w2*(c-d))) 

print "a,c,l,r:", a,c,l,r
gamma = (l + r)/np.sqrt(ha2ev)*np.sqrt(amu2au)*a2b
t = np.exp(-2*gamma)
print t
