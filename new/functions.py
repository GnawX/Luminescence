from decimal import *                                                                                              
import numpy as np                                                                                                 
from math import exp, sqrt, pi                                                                                     

getcontext().prec = 28

def fcf(ni,nf,wi,wf,k):

   "Franck-Condon overlap <nf|ni> with n quanta \
    P. Ruhoff Chemical Physicas 1994, 186, 355. \
    w = omega/h_bar reduced frequency \          
    q' = q + k normal coordinate shift"          
                                                 
   # constants                                   
   hbar = 6.582119514e-16 # eV.s/rad             
   ev = 1.6021766208e-19 # J                     
   amu = 1.660539040e-27 # kg                    
   ang = 1.0e-10 # m                             
   fr = ev/amu/ang**2                            
   # convert units                               
   wi = (wi/hbar)**2/fr/wi                       
   wf = (wf/hbar)**2/fr/wf                       
                                                 
   # definie parameters                          
   a = (wi - wf)/(wi + wf)                       
   b = 2*k*sqrt(wi)*wf/(wi + wf)                 
   c = -a                                        
   d = -2*k*sqrt(wf)*wi/(wi + wf)                
   e = 4*sqrt(wi*wf)/(wi + wf)                   
                                                 
   f = np.zeros((ni,nf))                         
   f[0,0] = sqrt(e/2)*exp(b*d/2/e)               
   f[1,0] = 1/sqrt(2)*b*f[0,0]                   
   for i in range(1,nf):                         
       if i - 2 < 0:                             
          f[0,i] = 1/sqrt(2*i)*d*f[0,i-1]        
       else:                                     
          f[0,i] = 1/sqrt(2*i)*d*f[0,i-1] + \    
                   sqrt(float(i-1)/i)*c*f[0,i-2] 
   for i in range(1,ni):                         
     for j in range(1,nf):                       
         if i - 2 < 0:                           
            f[i,j] = 1/sqrt(2*i)*b*f[i-1,j] + \  
                     sqrt(float(j)/i)/2*e*f[i-1,j-1]
         else:                                      
            f[i,j] = 1/sqrt(2*i)*b*f[i-1,j] + \     
                     sqrt(float(i-1)/i)*a*f[i-2,j] + \
                     sqrt(float(j)/i)/2*e*f[i-1,j-1]  
   return f                                           

def fcf_decimal(ni,nf,wi,wf,k):

   "Franck-Condon overlap <nf|ni> with n quanta \
    P. Ruhoff Chemical Physicas 1994, 186, 355. \
    w = omega/h_bar reduced frequency \          
    q' = q + k normal coordinate shift"          
   wi = Decimal(wi)                              
   wf = Decimal(wf)                              
   k = Decimal(k)                                
   # constants                                   
   hbar = Decimal('6.582119514e-16') # eV.s/rad  
   ev = Decimal('1.6021766208e-19') # J          
   amu = Decimal('1.660539040e-27') # kg         
   ang = Decimal('1.0e-10') # m                  
   fr = ev/amu/ang**2                            
   # convert units                               
   wi = (wi/hbar)**2/fr/wi                       
   wf = (wf/hbar)**2/fr/wf                       
                                                 
   # definie parameters                          
   a = (wi - wf)/(wi + wf)                       
   b = 2*k*wi.sqrt()*wf/(wi + wf)
   c = -a
   d = -2*k*wf.sqrt()*wi/(wi + wf)
   e = 4*(wi*wf).sqrt()/(wi + wf)

   f = np.zeros((ni,nf),dtype=Decimal)
   f[0,0] = (e/2).sqrt()*(b*d/2/e).exp()
   f[1,0] = 1/Decimal(2).sqrt()*b*f[0,0]
   for i in range(1,nf):
       if i - 2 < 0:
          f[0,i] = 1/Decimal(2*i).sqrt()*d*f[0,i-1]
       else:
          f[0,i] = 1/Decimal(2*i).sqrt()*d*f[0,i-1] + \
                   (Decimal(i-1)/Decimal(i)).sqrt()*c*f[0,i-2]
   for i in range(1,ni):
     for j in range(1,nf):
         if i - 2 < 0:
            f[i,j] = 1/Decimal(2*i).sqrt()*b*f[i-1,j] + \
                     (Decimal(j)/Decimal(i)).sqrt()/2*e*f[i-1,j-1]
         else:
            f[i,j] = 1/Decimal(2*i).sqrt()*b*f[i-1,j] + \
                     (Decimal(i-1)/Decimal(i)).sqrt()*a*f[i-2,j] + \
                     (Decimal(j)/Decimal(i)).sqrt()/2*e*f[i-1,j-1]
   return f.astype(float)

def gaussian(x,u,sigma):
    return 1/sigma/sqrt(pi*2)*exp(-(x-u)**2/2/sigma**2)
