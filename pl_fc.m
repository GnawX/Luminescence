% use fc_ruhoff
clear;
clc;
ni = 20;
nf = 200;
wi = 0.0174;
wf = 0.0122;
k = 9.157;
sigma = 0.03;
T = 300; % K
T = T*8.621738e-5; 
zpl = 4.0;
ne = 401;
%
f = fc_ruhoff(ni,nf,wi,wf,k);
fci = f.^2;
%fci = fci./max(max(fci)); % normalize
e = linspace(1.3,4,ne)';
p = zeros([ne,1]);
for ie=1:ne
    for i=0:ni-1
        for j=0:nf-1
            ei = i*wi;
            ef = j*wf;
            bz = exp(-ei/T);
            if bz < 1.0e-12
                bz = 0;
            end
            p(ie) = p(ie)+e(ie)^3*bz*fci(i+1,j+1) ...
                *lorentzian(e(ie),zpl+ei-ef,sigma);
        end
    end
end