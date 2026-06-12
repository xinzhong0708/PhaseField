function PARAM = Compute_L(STATE,PARAM,PHYS)

%COMPUTE_L_CONSTDCEQ Simple stable Allen-Cahn mobility.
%
% Uses:
%   L = 4*m/(3*kap*(dceq^2/M))
%
% with constant dceq. No local de is used.

dceq = 0.1;

if isfield(PHYS,'dceq')
    dceq = PHYS.dceq;
end
if isfield(PARAM,'dceq')
    dceq = PARAM.dceq;
end

if isfield(PHYS,'kap') && ~isempty(PHYS.kap)
    kap = PHYS.kap;
else
    kap = PHYS.kappa;
end

Msum = zeros(size(PARAM.M{1}));

for ie = 1:numel(PARAM.M)
    Msum = Msum + PARAM.M{ie};
end

Meff = Msum / numel(PARAM.M);

PARAM.L  = PARAM.L_fac*4*PHYS.m*Meff ./ (3*kap*dceq^2);
PARAM.Lm = PARAM.L * PHYS.m;
PARAM.LK = PARAM.L * kap;

end