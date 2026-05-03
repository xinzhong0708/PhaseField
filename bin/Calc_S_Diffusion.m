function [src] = Calc_S_Diffusion(p,po,e,dt_phy)
ny      =  size(p,1);
nx      =  size(p,2);
Np      =  size(p,3);
Ne      =  length(e{1});
src     =  repmat({zeros(ny,nx)},1,Ne);
for ip = 1:Np
    dp_ip        = 0*(p(:,:,ip) - po(:,:,ip)) / dt_phy;
    for ie = 1:Ne
        src{ie}  =  src{ie} + dp_ip .* e{ip}{ie};
    end
end
end