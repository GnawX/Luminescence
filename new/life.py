from math import sqrt, pi

#parameters

eps = 3.28   # dielectric constant
omega = 2.25  # transition energy eV
mu2 = 0.806   # dipole moment square  e^2 bohr^2
fc = 0.0499
eh = 27.211385 # eV
tau = 2.418884e-17 # hbar/Eh
#       sqrt(eps)*omega^3*mu2
# r =  -----------------------* FC
#        3*pi*eps0*hbar*c^3

# c = 137 in au  light speed

t = 3*pi*137**3/sqrt(eps)/(omega/eh)**3/mu2/fc*tau

print t
