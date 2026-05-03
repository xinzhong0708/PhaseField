function [s] = Calc_S_AllenCahn(phi,p,LL,F,omg)
%Initialize
Np = size(p,3);
ny = size(p,1); 
nx = size(p,2); 
s  = repmat({zeros(ny,nx)},1,Np); 
%Allen Cahn source
for alp = 1:Np
    s{alp}    = -LL.*F.dgdphi(phi(:,:,alp));
    for beta = 1:Np
        s{alp}=  s{alp} - LL.*F.dpdphi(alp,beta,phi).*omg(:,:,beta);
    end
end
end