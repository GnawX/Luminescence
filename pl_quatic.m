% use fc_ruhoff2
% anharmonic oscillator V = 1/2 omega^2 Q^2 + a Q^3
clc;
n1 = 20;
n2 = 200;
w1 = 0.0174; % eV
w2 = 0.0122;
a1 = 0e-4; % ev/[sqrt(amu)*Ang]^4
a2 = 7.67e-5;
k = 9.157; % sqrt(amu)*Ang
sigma = 0.03;
T = 300; % K
T = T*8.621738e-5; 
zpl = 4.0;
ne = 301;
%unit conversion
hbar = 6.582119514e-16;
ev = 1.6021766208e-19;
amu = 1.660539040e-27;
ang = 1.0e-10;
fr = ev/amu/ang^2;
b1 = sqrt(w1*fr/2/(w1/hbar)^2); % b = sqrt(hbar/2/omega)
b2 = sqrt(w2*fr/2/(w2/hbar)^2);
%
f = fc_ruhoff2(n1,n2,w1,w2,k);

e = linspace(1.3,4.2,ne)';
p = zeros([ne,1]);
for ie=1:ne
    for i=0:n1-4
        for j=0:n2-4
            ei = (i+0.5)*w1+6*a1*b1^4*((i+0.5)^2+0.25);
            ef = w2*(j+0.5)+6*a2*b2^4*((j+0.5)^2+0.25);
            bz = exp(-ei/T);
            if bz < 1.0e-12
                bz = 0;
            end
            fci = f(id(i),id(j))+...
                a2*b2^4/w2*(sqrt(j*(j-1)*(j-2)*(j-3))/4*fc(i,j-4,f)...
                +(2*j-1)*sqrt(j*(j-1))*fc(i,j-2,f)...
                -(2*j+3)*sqrt((j+1)*(j+2))*fc(i,j+2,f)...
                -sqrt((j+1)*(j+2)*(j+3)*(j+4))*fc(i,j+4,f));           
            fci = fci^2;
            p(ie) = p(ie)+e(ie)^3*bz*fci*lorentzian(e(ie),zpl+ei-ef,sigma);
        end
    end
end

function [ r ] = fc( n,m,f )
    if n < 0 || m < 0
        r = 0.0;
    else
        r = f(id(n),id(m));
    end
end