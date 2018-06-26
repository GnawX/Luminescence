function [ f ] = fc_chang( n1,n2 )
%Franck-Condon factor J Mol. Spectroscopy 232 (2005) 102-104
d = 4.524;
a1 = 0.010;
a2 = 0.018;
%constants
hbar = 6.582119514e-16;
ev = 1.6021766208e-19;
amu = 1.660539040e-27;
ang = 1.0e-10;
fr = ev/amu/ang^2;
a1 = (a1/hbar)^2/fr/a1;
a2 = (a2/hbar)^2/fr/a2;
% parameters
S = a1*a2*d^2/(a1+a2);
b1 = -a2*sqrt(a1)*d/(a1+a2);
b2 = a1*sqrt(a2)*d/(a1+a2);
A = 2*sqrt(a1*a2)/(a1+a2);
% calculate the fc factors
temp = 0.0;
for k1 = 0:n1
    for k2 = 0:n2
        if mod(k1+k2,2) == 0
            sq = k1+k2-1:-2:0;
            ik = prod(sq)/(a1+a2)^((k1+k2)/2);
        else
            ik = 0;
        end
        temp = temp + nchoosek(n1,k1)*nchoosek(n2,k2)...
            *hermite(n1-k1,b1)*hermite(n2-k2,b2)*(2*sqrt(a1))^k1...
            *(2*sqrt(a2))^k2*ik;
    end
end
f = temp^2*A*exp(-S)/2^(n1+n2)/factorial(n1)/factorial(n2);
           
end

