function PARAM = Calc_Kappa_SpinodalBand(STATE,PARAM,MODEL,GRID,PHYS)
%CALC_KAPPA_SPINODALBAND  Spatially varying kappa for spinodal regions.
%
% Purpose
% -------
% Build PARAM.kappa_eff as a spatial map:
%
%   kappa_eff(x) = kappa0  inside / near nodes where min eig(H_c) < 0
%                = 0       away from spinodal curvature
%
% This is useful when you want the fourth-order Cahn-Hilliard stabilization
% only in thermodynamically unstable regions, instead of applying kappa
% everywhere and creating unwanted uphill diffusion in convex regions.
%
% Typical usage in the main time loop
% -----------------------------------
%   STATE_OLD = LE_Run_Mode_SmoothGP(STATE_OLD,PARAM,MODEL);
%
%   PARAM = Calc_Kappa_SpinodalBand(STATE_OLD,PARAM,MODEL,GRID,PHYS);
%
%   STATE_RAW = PF_Coupled_ACCH_LETangent(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);
%
% Important
% ---------
% This function only detects curvature from the current STATE.c and STATE.p.
% It does not change STATE.  It only writes PARAM.kappa_eff and diagnostics.

% ---------------- user controls ----------------
% Base kappa value
kappa0 = GetOpt(PARAM,'kappa_spinodal_value',PHYS.kappa);

% A phase only contributes to the spinodal mask where its phase weight is
% larger than pmin.  This avoids tiny phase tails activating kappa.
pmin = GetOpt(PARAM,'kappa_spinodal_pmin',1e-3);

% Spinodal criterion:
%   min eig(H_c) < Htol
% Usually Htol = 0.  Use a small positive value to activate kappa slightly
% before H actually becomes negative.
Htol = GetOpt(PARAM,'kappa_spinodal_Htol',0);

% Number of grid cells around the spinodal mask where kappa remains active.
% This gives "within and close to H<0 grids".
halo = GetOpt(PARAM,'kappa_spinodal_halo',2);

% Optional smoothing passes after halo construction.
% 0 = sharp binary band.
% 1 or 2 = smoother transition.
smooth_pass = GetOpt(PARAM,'kappa_spinodal_smooth',1);

% Minimum weight below which kappa is set exactly to zero.
cutoff = GetOpt(PARAM,'kappa_spinodal_cutoff',1e-8);

% ---------------- detect spinodal nodes ----------------
p = STATE.p;
[ny,nx,Ngrain] = size(p);
N = ny*nx;

spinodal = false(ny,nx);
minH_map = inf(ny,nx);

for ig = 1:Ngrain

    % Skip phases that have no composition Hessian.
    R = PhaseThermo(MODEL.pars{ig},STATE.c{ig});
    if isempty(R.H_c)
        continue
    end

    H = 0.5*(R.H_c + permute(R.H_c,[2 1 3]));
    Nc = size(H,1);

    % Compute minimum Hessian eigenvalue at every grid node.
    hmin = inf(1,N);
    for n = 1:N
        hmin(n) = min(eig(H(:,:,n)));
    end
    hmin = reshape(hmin,ny,nx);

    active = p(:,:,ig) > pmin;
    bad    = active & (hmin < Htol);

    spinodal = spinodal | bad;
    minH_map = min(minH_map,hmin);

end

% ---------------- expand to a small band around spinodal nodes ----------------
band = spinodal;

if halo > 0
    kernel = ones(2*halo+1,2*halo+1);
    band = conv2(double(spinodal),kernel,'same') > 0;
end

% ---------------- smooth the band to avoid a sharp kappa jump ----------------
w = double(band);

for it = 1:smooth_pass
    % Simple 3x3 averaging kernel, no toolbox required.
    w = conv2(w,ones(3,3)/9,'same');

    % Do not let smoothing create a global low-amplitude background.
    w(~band) = 0;
end

w(w < cutoff) = 0;
w = min(max(w,0),1);

% ---------------- write PARAM.kappa_eff ----------------
kappa_map = kappa0 .* w;

% Try to preserve the type/shape convention already used by the code.
% Most current code uses PARAM.kappa_eff as a ny-by-nx numeric array.
if isfield(PARAM,'kappa_eff') && iscell(PARAM.kappa_eff)
    for i = 1:numel(PARAM.kappa_eff)
        PARAM.kappa_eff{i} = kappa_map;
    end
else
    PARAM.kappa_eff = kappa_map;
end

% Diagnostics for plotting/debugging
PARAM.kappa_spinodal_mask   = spinodal;
PARAM.kappa_spinodal_band   = band;
PARAM.kappa_spinodal_weight = w;
PARAM.kappa_spinodal_minH   = minH_map;

end


function value = GetOpt(S,name,default)
if isfield(S,name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default;
end
end
