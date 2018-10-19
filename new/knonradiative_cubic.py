# energy of excited state E = 1/2 omega^2 Q^2 + lambda Q^3                                                                                                                   
import numpy as np                                                                                                                                                           
from math import sqrt, exp, pi                                                                                                                                               
from functions import fcf_decimal, gaussian                                                                                                                                  

""" calculate the nonradiative recombination rate """

# input parameters
ni    =    50                   # initial
nf    =    350                   # final 
wi    =    0.009937   # in eV            
wf    =    0.010445                      
q     =    9.04728    # in sqrt(amu)*ang 
lam   =    6.88757e-4 # lambda cubic term coeff unit eV/(A sqrt(amu))^3 
sigma =    0.38     # broadening                                        
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
beta = sqrt(wi*fr/(wi/hbar)**2/2)  # sqrt(hbar/2omega) unit Q
#
ff = np.zeros((ni,nf))

rate = 0.0
for i in range(ni-4):
    ei = (0.5 + i)*wi - 30*lam**2*beta**6/wi*(i**2+i+11./30)
    bz = exp(-ei/t)/z
    for j in range(nf):
        ef = (0.5 + j)*wf
        if i < 1:
           ff[i,j] = beta*sqrt(i+1)*f[i+1,j] + lam*beta**4/3/wi* \
                     (-9*(2*i+1)*f[i,j]-2*(5*i+6)*sqrt((i+1)*(i+2))*f[i+2,j] - \
                     sqrt((i+1)*(i+2)*(i+3)*(i+4))*f[i+4,j])
        elif i >= 1 and i < 2:
           ff[i,j] = beta*(sqrt(i+1)*f[i+1,j] + sqrt(i)*f[i-1,j]) + lam*beta**4/3/wi* \
                     (-9*(2*i+1)*f[i,j]-2*(5*i+6)*sqrt((i+1)*(i+2))*f[i+2,j] - \
                     sqrt((i+1)*(i+2)*(i+3)*(i+4))*f[i+4,j])
        elif i >= 2 and i < 4:
           ff[i,j] = beta*(sqrt(i+1)*f[i+1,j] + sqrt(i)*f[i-1,j]) + lam*beta**4/3/wi* \
                     (2*(5*i-1)*sqrt(i*(i-1))*f[i-2,j]-9*(2*i+1)*f[i,j] - \
                     2*(5*i+6)*sqrt((i+1)*(i+2))*f[i+2,j] - \
                     sqrt((i+1)*(i+2)*(i+3)*(i+4))*f[i+4,j])
        else:
           ff[i,j] = beta*(sqrt(i+1)*f[i+1,j] + sqrt(i)*f[i-1,j]) + lam*beta**4/3/wi* \
                     (sqrt(i*(i-1)*(i-2)*(i-3))*f[i-4,j] + \
                     2*(5*i-1)*sqrt(i*(i-1))*f[i-2,j]-9*(2*i+1)*f[i,j] - \
                     2*(5*i+6)*sqrt((i+1)*(i+2))*f[i+2,j] - \
                     sqrt((i+1)*(i+2)*(i+3)*(i+4))*f[i+4,j])

        rate += 2*pi*bz*ff[i,j]**2*zpl**2*wij**2*gaussian(zpl+ei,ef,sigma)

print rate/h  # Hz
