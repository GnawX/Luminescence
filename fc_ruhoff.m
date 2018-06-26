function [ f ] = fc_ruhoff( ni,nf,wi,wf,k )
%Franck-Condon integral <nf|ni> recursion 
%P. Ruhoff Chemical Physicas 1994, 186, 355.
% frequency wi wf in eV
% k normal coordinate in sqrt(amu)*ang
% test
%wi = 0.018;
%wf = 0.018;
%k = 4.524;
%
hbar = 6.582119514e-16;
ev = 1.6021766208e-19;
amu = 1.660539040e-27;
ang = 1.0e-10;
fr = ev/amu/ang^2;
wi = (wi/hbar)^2/fr/wi;
wf = (wf/hbar)^2/fr/wf;
% define parameters
a = (wi - wf)/(wi + wf);
b = 2.*k*sqrt(wi)*wf/(wi + wf);
c = -a;
d = -2.*k*sqrt(wf)*wi/(wi + wf);
e = 4.*sqrt(wi*wf)/(wi + wf);
% calculate overlap
f = zeros([ni,nf]);
f(1,1) = sqrt(e/2.)*exp(b*d/2./e);
f(2,1) = 1./sqrt(2.)*b*f(1,1);
for i=1:nf-1
    if i - 2 < 0
       f(1,i+1) = 1./sqrt(2.*i)*d*f(1,i);
    else
       f(1,i+1) = 1./sqrt(2.*i)*d*f(1,i) + sqrt((i-1.)/i)*c*f(1,i-1);
    end
end
for i=1:ni-1
      for j=1:nf-1
          if i - 2 < 0
             f(i+1,j+1) = 1./sqrt(2.*i)*b*f(i,j+1)...
                 + 0.5*sqrt((j-0.)/i)*e*f(i,j);
          else
             f(i+1,j+1) = 1./sqrt(2.*i)*b*f(i,j+1)...
                 + sqrt((i-1.)/i)*a*f(i-1,j+1)...
                 + 0.5*sqrt((j-0.)/i)*e*f(i,j);
          end
      end
end
                                       
end
 

