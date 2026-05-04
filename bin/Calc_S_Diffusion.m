function [src] = Calc_S_Diffusion(p,po,e,dt_phy)
%CALC_S_DIFFUSION Source term caused by phase-fraction change.
%
% KKS-style linearized composition equation:
%
%   dE/dt = diffusion - sum_i e_i * dp_i/dt
%
% The diffusion solver uses:
%
%   (E^n - E^o)/dt = diffusion - src
%
% therefore:
%
%   src = sum_i e_i * (p_i^n - p_i^o)/dt

ny = size(p,1);
nx = size(p,2);
Np = size(p,3);
Ne = length(e{1});

src = repmat({zeros(ny,nx)},1,Ne);

for ip = 1:Np
    dp_ip = (p(:,:,ip) - po(:,:,ip)) / dt_phy;

    for ie = 1:Ne
        src{ie} = src{ie} + dp_ip .* e{ip}{ie};
    end
end

end