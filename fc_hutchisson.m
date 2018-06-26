function [ f ] = fc_hutchisson( n1,n2)
%Phys. Rev. 36, 410 (Hutchisson)
% n1, n2 initial and final quantum number
% w1, w2 in eV (hbar*omega)
% q normal coordinate sqrt(amu)*ang
% test
w1 = 0.010; % omega/hbar
w2 = 0.018;
k = 4.524;
% unit conversion
hbar = 6.582119514e-16;
ev = 1.6021766208e-19;
amu = 1.660539040e-27;
ang = 1.0e-10;
fr = ev/amu/ang^2;
w1 = (w1/hbar)^2/fr/w1; % omega/hbar
w2 = (w2/hbar)^2/fr/w2;
% dimensionless parameters
a = sqrt(w2/w1);
d = k*sqrt(w2);

temp = sqrt(2*a/(a^2+1)*factorial(n1)*factorial(n2)/2^(n1+n2)) ...
    *exp(-d^2/2/(a^2+1));
f = 0;
for l=0:min(n1,n2)
    for i=0:floor((n1-l)/2)
        for j=0:floor((n2-l)/2)
            f = f + (4*a/(1+a^2))^l/factorial(l)*((1-a^2)/(1+a^2))^i/ ...
                factorial(i)*((a^2-1)/(1+a^2))^j/factorial(j)* ...
                (-2*a*d/(1+a^2))^(n1-2*i-l)/factorial(n1-2*i-l)* ...
                (2*d/(1+a^2))^(n2-2*j-l)/factorial(n2-2*j-l);            
        end
    end
end
f = f*temp;
end

