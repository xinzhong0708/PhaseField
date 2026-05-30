function PARAM = Calc_Kappa_SpinodalBand(STATE,PARAM,MODEL,GRID,PHYS)
%CALC_KAPPA_SPINODALBAND Spatial kappa switch from det(H_c), expanded by a few cells.
%
% Rule:
%   det(H_c) < 0  -> true spinodal region
%   expanded true spinodal region -> PARAM.kappa_eff = PHYS.kappa
%   elsewhere -> PARAM.kappa_eff = 0
%
% PhaseThermo receives unpacked c{ip}{ic} as 1-by-N.

ny = size(STATE.p,1);
nx = size(STATE.p,2);
N  = ny*nx;

p = UnpackP(STATE.p);     % 1-by-N-by-Nphase
c = UnpackC(STATE.c);     % c{ip}{ic} is 1-by-N

kappa_eff          = zeros(ny,nx);
spinodal_mask_true = false(ny,nx);

if isfield(PARAM,'kappa_spinodal_pmin') && ~isempty(PARAM.kappa_spinodal_pmin)
    pmin = PARAM.kappa_spinodal_pmin;
else
    pmin = 0;
end

if isfield(PARAM,'kappa_spinodal_expand') && ~isempty(PARAM.kappa_spinodal_expand)
    ngrow = PARAM.kappa_spinodal_expand;
else
    ngrow = 0;
end

for ip = 1:numel(c)

    % IMPORTANT: use unpacked c, not STATE.c
    R = PhaseThermo(MODEL.pars{ip},c{ip});

    if isempty(R.H_c)
        continue
    end

    H = 0.5*(R.H_c + permute(R.H_c,[2 1 3]));

    bad = false(1,N);

    for i = 1:N
        bad(i) = det(H(:,:,i)) < 0;
    end

    p_ip = reshape(p(:,:,ip),1,N);
    bad  = bad & (p_ip > pmin);

    bad2 = reshape(bad,ny,nx);

    spinodal_mask_true = spinodal_mask_true | bad2;
end

% Expand the spinodal region by ngrow grid cells
if ngrow > 0
    kernel = ones(2*ngrow + 1, 2*ngrow + 1);
    spinodal_mask_expanded = conv2(double(spinodal_mask_true), kernel, 'same') > 0;
else
    spinodal_mask_expanded = spinodal_mask_true;
end

kappa_eff(spinodal_mask_expanded) = PHYS.kappa;

PARAM.kappa_eff = kappa_eff;

% Diagnostics
PARAM.kappa_spinodal_mask_true     = spinodal_mask_true;
PARAM.kappa_spinodal_mask_expanded = spinodal_mask_expanded;
PARAM.kappa_spinodal_mask          = spinodal_mask_expanded;

end


function p = UnpackP(p)
[ny,nx,np] = size(p);
p = reshape(p,1,ny*nx,np);
end


function c = UnpackC(c)
for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        c{ip}{ic} = reshape(c{ip}{ic},1,[]);
    end
end
end