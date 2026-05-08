function STATE = Calc_S_Diffusion(STATE,STATE_OLD,NUM)
%CALC_S_DIFFUSION Source term caused by phase-fraction change.
%   dE/dt = diffusion - sum_i e_i * dp_i/dt
%   (E^n - E^o)/dt = diffusion - src
%   src = sum_i e_i * (p_i^n - p_i^o)/dt
p      = STATE.p;
po     = STATE_OLD.p;
e      = STATE.e;
dt_phy = NUM.dt_phy;
ny     = size(p,1);
nx     = size(p,2);
Np     = size(p,3);
Ne     = length(e{1});
src    = repmat({zeros(ny,nx)},1,Ne);
for ip = 1:Np
    dp_ip = (p(:,:,ip) - po(:,:,ip)) / dt_phy;
    for ie = 1:Ne
        src{ie} = src{ie} + dp_ip .* e{ip}{ie};
    end
end
STATE.S_DIFF = src;
end