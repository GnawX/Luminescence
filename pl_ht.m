clear;
clc;
ni = 10;
nf = 120;
wi = 0.0181;
wf = 0.0139;
k = 9.157;
sigma = 0.03;
T = 300; % K
T = T*8.621738e-5; 
zpl = 4.0;
ne = 401;
%
f = fc_ruhoff(ni,nf,wi,wf,k);
e = linspace(1.3,4,ne)';
p = zeros([ne,1]);
%for ie=1:ne
%    for i=0:ni-2
%        for j=0:nf-1
%            es = i*wi;
%            eg = j*wf;
%            if i - 1 < 0
%                ht = (sqrt(i+1)*f(i+2,j+1))^2;
%            else
%                ht = (sqrt(i)*f(i,j+1) + sqrt(i+1)*f(i+2,j+1))^2;
%            end
%            p(ie) = p(ie)+e(ie)^3*exp(-es/T)*ht ...
%                *lorentzian(e(ie),zpl+es-eg,sigma);
%        end
%    end
%end
for ie=1:ne
    for i=0:ni-2
        for j=0:nf-2
            es = i*wi;
            eg = j*wf;
            if i == 0
                if j == 0
                    ht = (sqrt(i+1)*f(i+2,j+1) + sqrt(j+1)*f(i+1,j+2))^2;
                else
                    ht = (sqrt(i+1)*f(i+2,j+1) + sqrt(j+1)*f(i+1,j+2) + ...
                        sqrt(j)*f(i+1,j))^2;
                end
            else
                if j == 0
                    ht = (sqrt(i+1)*f(i+2,j+1) + sqrt(j+1)*f(i+1,j+2) + ...
                        sqrt(i)*f(i,j+1))^2;
                else
                    ht = (sqrt(i+1)*f(i+2,j+1) + sqrt(j+1)*f(i+1,j+2) + ...
                        sqrt(i)*f(i,j+1) + sqrt(j)*f(i+1,j))^2;
                end
            end
            bz = exp(-es/T);
            if bz < 1.0e-10
                bz = 0;
            end
            p(ie) = p(ie)+e(ie)^3*bz*ht ...
                *lorentzian(e(ie),zpl+es-eg,sigma);
        end
    end
end