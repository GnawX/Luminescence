# recurrence method to calculate FC integrals                                             
# Schmidt, Molecular Physics 108,1513-1529, 2010.                                         
import numpy as np                                                                        
from math import sqrt, exp, pi                                                            
import sys

# parameters a = omega/hbar in 1/amu/Ang**2, d in sqrt(amu)*Ang
omg1,omg2 = 0.010,0.018  # phonon frequency in eV, ground, excited            
a,a_ = 2.3923, 4.3061                                       
d = 4.524                                          
zpl = 2.66  # zero phonon line eV                              
T = 0.025   # temperature in eV                                
sigma = 0.2 # gauusian broadenning                             
m,n = 50, 50  # vibration states gs, es                       
v1 = int(sys.argv[1])
v2 = int(sys.argv[2])

# eight coefficients
c1 = -a*sqrt(a_)*d/sqrt(2.0)/(a+a_)
c2 = -(a-a_)/(a+a_)                
c3 = sqrt(a)*a_*d/sqrt(2.0)/(a+a_) 
c4 = 2.0*sqrt(a*a_)/(a+a_)         
c5 = -c3                           
c6 = -c2                           
c7 = -c1                           
c8 = 1.0                           

def gaussian(x,u,sigma):
    return 1./sigma/sqrt(2.*pi)*exp(-0.5*(x-u)**2/sigma**2)

def thermal(E,T):
    return exp(-E/T)

# calculate the FC integrals

f = np.zeros((m+1,n+1))
f[0,0] = sqrt(c4)*exp(-a*a_*d**2/2.0/(a + a_))
if m == 0:
   f[0,1] = 2.0*c1*f[0,0]
   for i in range(n-1):
       f[0,i+2] = sqrt((i+1.0)/(i+2.0))*c2*f[0,i] + sqrt(1.0/(i+2.0))*2.0*c1*f[0,i+1]
elif n == 0:
   f[1,0] = 2.0*c3*f[0,0]
   for i in range(m-1):
       f[i+2,0] = sqrt((i+1.0)/(i+2.0))*c6*f[i,0] + sqrt(1.0/(i+2.0))*2.0*c3*f[i+1,0]
else:
   f[0,1] = 2.0*c1*f[0,0]
   f[1,0] = 2.0*c3*f[0,0]
   for i in range(n-1):
       f[0,i+2] = sqrt((i+1.0)/(i+2.0))*c2*f[0,i] + sqrt(1.0/(i+2.0))*2.0*c1*f[0,i+1]
   for i in range(m-1):
       f[i+2,0] = sqrt((i+1.0)/(i+2.0))*c6*f[i,0] + sqrt(1.0/(i+2.0))*2.0*c3*f[i+1,0]
   f[1,1] = c1*f[1,0] + c3*f[0,1] + c4*f[0,0]
   f[1,2] = c1/sqrt(2.0)*f[1,1] + c2/sqrt(2.0)*f[1,0] + c3*f[0,2] + c4/sqrt(2.0)*f[0,1] + c5/sqrt(2.0)*f[0,0]
   f[2,1] = c1*f[2,0] + c3/sqrt(2.0)*f[1,1] + c4/sqrt(2.0)*f[1,0] + c6/sqrt(2.0)*f[0,1] + c7/sqrt(2.0)*f[0,0]
   for i in range(1,m):
       for j in range(1,n):
           f[i+1,j+1] = (c1/sqrt(j+1.0)*f[i+1,j] + c2*sqrt(j/(j+1.0))*f[i+1,j-1] + c3/sqrt(i+1.0)*f[i,j+1] +
                        c4/sqrt((i+1.0)*(j+1.0))*f[i,j] + c5*sqrt(j/(i+1.0)/(j+1.0))*f[i,j-1] +
                        c6*sqrt(i/(i+1.0))*f[i-1,j+1] + c7*sqrt(i/(i+1.0)/(j+1.0))*f[i-1,j] +
                        c8*sqrt(i*j/(i+1.0)/(j+1.0))*f[i-1,j-1])
