function PARAM = Calc_GB_Eta(STATE,PARAM,PHYS)

%CALC_GB_ETA Make eta low at grain/interface network and high in interiors.
%
% Low eta means large local storage capacity for off-formula composition:
%   R = E - sum_alpha p_alpha e_alpha
%   mu ~= eta * R
%
% Therefore eta_gb << eta_solid mimics a flexible GB/fluid reservoir.

p = STATE.p;

[ny,nx,~] = size(p);

if isfield(PARAM,'gb_eta_fac')
    gb_eta_fac = PARAM.gb_eta_fac;
else
    gb_eta_fac = 1e-3;
end

if isfield(PARAM,'gb_smooth_iter')
    gb_smooth_iter = PARAM.gb_smooth_iter;
else
    gb_smooth_iter = 2;
end

if isfield(PARAM,'gb_cut')
    gb_cut = PARAM.gb_cut;
else
    gb_cut = 1e-6;
end

% Grain/interface mixedness.
% B = 0 in pure interiors, positive where more than one grain/phase coexists.
B = 1 - sum(p.^2,3);
B = max(B,0);

% Widen/smooth the GB network slightly.
for it = 1:gb_smooth_iter
    B = conv2(B,ones(3)/9,'same');
end

if max(B(:)) > 0
    B = B/max(B(:));
end

B(B < gb_cut) = 0;

% Solid eta.
if isfield(PHYS,'eta')
    eta_solid = PHYS.eta;
else
    eta_solid = PARAM.eta;
end

if isscalar(eta_solid)
    eta_solid = eta_solid*ones(ny,nx);
end

eta_gb = gb_eta_fac .* eta_solid;

% Smooth interpolation: high eta in grain interior, low eta at GB.
PARAM.gb_mask = B;
PARAM.eta     = eta_solid.*(1-B) + eta_gb.*B;

end