function [STATE_BG,DIAG] = PF_CH_Coarse_Background_CS(STATE_REF,PARAM,MODEL,GRID,PHYS,NUM,STATE_COEF)
%PF_CH_COARSE_BACKGROUND_CS Coarse-grid global CH background correction.
%
% Purpose:
%   Keep the real STATE on the fine grid, but obtain a smooth global
%   chemical-potential update from a temporary coarse-grid fixed-p CH solve.
%
% Output:
%   STATE_BG has the same fine grid as STATE_REF.
%   Its phi and p are unchanged from STATE_REF.
%   Its mu_e and E include the prolonged coarse-grid dmu background.
%
% Controls:
%   NUM.CH_coarse_ratio      default 4
%   NUM.CH_coarse_kappac     default 0
%   NUM.CH_bg_theta          default 1
%   NUM.CH_bg_mean_correct   default 0
%
% This function calls your existing PF_CH_LECorrector_FixedP_Band_CS on the
% coarse grid with CHLE_force_full = 1.

if nargin < 7 || isempty(STATE_COEF)
    STATE_COEF = STATE_REF;
end

r = 1;
if isfield(NUM,'CH_coarse_ratio') && ~isempty(NUM.CH_coarse_ratio)
    r = NUM.CH_coarse_ratio;
end
r = max(1,round(r));

theta = 1;
if isfield(NUM,'CH_bg_theta') && ~isempty(NUM.CH_bg_theta)
    theta = NUM.CH_bg_theta;
end

if r == 1
    NUM_C = NUM;
    NUM_C.CHLE_force_full = 1;
    if isfield(NUM,'CH_coarse_kappac') && NUM.CH_coarse_kappac == 0
        PARAM_C = PARAM;
        PARAM_C.use_kappa_c = 0;
    else
        PARAM_C = PARAM;
    end

    [STATE_C,DIAG_C] = PF_CH_LECorrector_FixedP_Band_CS(STATE_REF,STATE_COEF,PARAM_C,MODEL,GRID,PHYS,NUM_C);

    STATE_BG = STATE_C;
    DIAG.coarse = DIAG_C;
    DIAG.ratio  = 1;
    return
end

% ------------------------------------------------------------
% Restrict fine states and coefficients to coarse grid
% ------------------------------------------------------------
[STATE_TIME_C,GRID_C] = Restrict_State_CH_Local(STATE_REF,GRID,r);
[STATE_IT_C,~]        = Restrict_State_CH_Local(STATE_COEF,GRID,r);

PARAM_C = Restrict_PARAM_CH_Local(PARAM,GRID,GRID_C,r);

if isfield(NUM,'CH_coarse_kappac') && NUM.CH_coarse_kappac == 0
    PARAM_C.use_kappa_c = 0;
end

% ------------------------------------------------------------
% Coarse fixed-p CH solve, full coarse domain
% ------------------------------------------------------------
NUM_C = NUM;
NUM_C.CHLE_force_full = 1;
NUM_C.CHLE_res_rel    = 0;

[STATE_C_CORR,DIAG_C] = PF_CH_LECorrector_FixedP_Band_CS(STATE_TIME_C,STATE_IT_C,PARAM_C,MODEL,GRID_C,PHYS,NUM_C);

% ------------------------------------------------------------
% Prolong only the coarse dmu correction to fine grid
% ------------------------------------------------------------
Ne  = numel(STATE_REF.E);
dmu = cell(1,Ne);

for ie = 1:Ne
    dmu_C   = STATE_C_CORR.mu_e{ie} - STATE_IT_C.mu_e{ie};
    dmu{ie} = theta * Prolong_Field_Local(dmu_C,size(STATE_REF.E{ie}));
end

% ------------------------------------------------------------
% Apply background dmu to fine state
% ------------------------------------------------------------
STATE_BG = STATE_REF;
STATE_BG.mu_e = STATE_REF.mu_e;

for ie = 1:Ne
    STATE_BG.mu_e{ie} = STATE_REF.mu_e{ie} + dmu{ie};
end

% Use the same implicit chi choice as the original solver.
if isfield(PARAM,'use_CS_chi') && PARAM.use_CS_chi == 1
    chi_use = ConvexifyChi_ForImplicit_Local(STATE_COEF.chi,PARAM);
else
    chi_use = STATE_COEF.chi;
end

STATE_BG.E = STATE_REF.E;

for ie = 1:Ne

    En = STATE_REF.E{ie};

    for je = 1:Ne
        En = En + chi_use{ie,je} .* dmu{je};
    end

    STATE_BG.E{ie} = En;

end

if isfield(NUM,'CH_bg_mean_correct') && NUM.CH_bg_mean_correct == 1
    STATE_BG.E = EnforceMeanE_Local(STATE_BG.E,STATE_REF.E);
end

% p/phi are fixed by construction.
STATE_BG.phi = STATE_REF.phi;
STATE_BG.p   = STATE_REF.p;

% Predict omega only as a background diagnostic. Final LE overwrites it.
if isfield(STATE_REF,'omg') && ~isempty(STATE_REF.omg)
    STATE_BG.omg = STATE_REF.omg;

    for ig = 1:size(STATE_REF.omg,3)

        domega = zeros(size(STATE_REF.omg(:,:,ig)));

        for ie = 1:Ne
            domega = domega - STATE_COEF.e{ig}{ie} .* dmu{ie};
        end

        STATE_BG.omg(:,:,ig) = STATE_REF.omg(:,:,ig) + domega;

    end

    omg_mean = mean(STATE_BG.omg,3);

    for ig = 1:size(STATE_BG.omg,3)
        STATE_BG.omg(:,:,ig) = STATE_BG.omg(:,:,ig) - omg_mean;
    end
end

STATE_BG.chi = STATE_COEF.chi;
STATE_BG.e   = STATE_COEF.e;
if isfield(STATE_COEF,'c')
    STATE_BG.c = STATE_COEF.c;
end
if isfield(STATE_COEF,'c_ext')
    STATE_BG.c_ext = STATE_COEF.c_ext;
end

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
max_dmu = 0;
max_dE  = 0;

for ie = 1:Ne
    max_dmu = max(max_dmu,max(abs(dmu{ie}(:))));
    max_dE  = max(max_dE,max(abs(STATE_BG.E{ie}(:)-STATE_REF.E{ie}(:))));
end

DIAG.coarse     = DIAG_C;
DIAG.ratio      = r;
DIAG.theta      = theta;
DIAG.max_dmu_bg = max_dmu;
DIAG.max_dE_bg  = max_dE;
DIAG.ny_coarse  = GRID_C.ny;
DIAG.nx_coarse  = GRID_C.nx;

end


% =========================================================================
% Local helpers
% =========================================================================

function [STATE_C,GRID_C] = Restrict_State_CH_Local(STATE,GRID,r)

[ny,nx,~] = size(STATE.p);

nyc = ceil(ny/r);
nxc = ceil(nx/r);

GRID_C    = GRID;
GRID_C.ny = nyc;
GRID_C.nx = nxc;
GRID_C.dy = GRID.dy*r;
GRID_C.dx = GRID.dx*r;

STATE_C = STATE;

% 2D cell fields
if isfield(STATE,'E')
    for ie = 1:numel(STATE.E)
        STATE_C.E{ie} = BlockMean2D_Local(STATE.E{ie},r);
    end
end

if isfield(STATE,'mu_e')
    for ie = 1:numel(STATE.mu_e)
        STATE_C.mu_e{ie} = BlockMean2D_Local(STATE.mu_e{ie},r);
    end
end

% 3D arrays
if isfield(STATE,'phi')
    STATE_C.phi = BlockMean3D_Local(STATE.phi,r);
end
if isfield(STATE,'p')
    STATE_C.p = BlockMean3D_Local(STATE.p,r);
end
if isfield(STATE,'omg') && ~isempty(STATE.omg)
    STATE_C.omg = BlockMean3D_Local(STATE.omg,r);
end

% Nested phase/element fields
if isfield(STATE,'e')
    for ig = 1:numel(STATE.e)
        for ie = 1:numel(STATE.e{ig})
            STATE_C.e{ig}{ie} = BlockMean2D_Local(STATE.e{ig}{ie},r);
        end
    end
end

if isfield(STATE,'c')
    for ig = 1:numel(STATE.c)
        for ic = 1:numel(STATE.c{ig})
            STATE_C.c{ig}{ic} = BlockMean2D_Local(STATE.c{ig}{ic},r);
        end
    end
end

if isfield(STATE,'c_ext') && ~isempty(STATE.c_ext)
    for ig = 1:numel(STATE.c_ext)
        for ic = 1:numel(STATE.c_ext{ig})
            STATE_C.c_ext{ig}{ic} = BlockMean2D_Local(STATE.c_ext{ig}{ic},r);
        end
    end
end

if isfield(STATE,'chi')
    for ie = 1:size(STATE.chi,1)
        for je = 1:size(STATE.chi,2)
            STATE_C.chi{ie,je} = BlockMean2D_Local(STATE.chi{ie,je},r);
        end
    end
end

end


function PARAM_C = Restrict_PARAM_CH_Local(PARAM,GRID,GRID_C,r)

PARAM_C = PARAM;

if isfield(PARAM,'M')
    for ie = 1:numel(PARAM.M)
        if isscalar(PARAM.M{ie})
            PARAM_C.M{ie} = PARAM.M{ie} * ones(GRID_C.ny,GRID_C.nx);
        else
            PARAM_C.M{ie} = BlockMean2D_Local(PARAM.M{ie},r);
        end
    end
end

names = {'L','LK','Lm','LL','LM','A_ac','eta','kappa_eff','WScale'};

for in = 1:numel(names)
    nm = names{in};
    if isfield(PARAM,nm) && ~isempty(PARAM.(nm)) && ~isscalar(PARAM.(nm))
        if ndims(PARAM.(nm)) == 2
            PARAM_C.(nm) = BlockMean2D_Local(PARAM.(nm),r);
        elseif ndims(PARAM.(nm)) == 3
            PARAM_C.(nm) = BlockMean3D_Local(PARAM.(nm),r);
        end
    end
end

if isfield(PARAM,'kappa_eff') && isscalar(PARAM.kappa_eff)
    PARAM_C.kappa_eff = PARAM.kappa_eff * ones(GRID_C.ny,GRID_C.nx);
end

end


function B = BlockMean2D_Local(A,r)

[ny,nx] = size(A);

nyc = ceil(ny/r);
nxc = ceil(nx/r);

B = zeros(nyc,nxc);

for iy = 1:nyc
    ys = (iy-1)*r + 1;
    ye = min(iy*r,ny);

    for ix = 1:nxc
        xs = (ix-1)*r + 1;
        xe = min(ix*r,nx);

        blk = A(ys:ye,xs:xe);
        v = blk(:);
        v = v(isfinite(v));
        if isempty(v)
            B(iy,ix) = 0;
        else
            B(iy,ix) = mean(v);
        end
    end
end

end


function B = BlockMean3D_Local(A,r)

[ny,nx,nz] = size(A);

nyc = ceil(ny/r);
nxc = ceil(nx/r);

B = zeros(nyc,nxc,nz);

for iz = 1:nz
    B(:,:,iz) = BlockMean2D_Local(A(:,:,iz),r);
end

end


function F = Prolong_Field_Local(C,sz)

ny = sz(1);
nx = sz(2);

[nyc,nxc] = size(C);

if nyc == ny && nxc == nx
    F = C;
    return
end

if nyc == 1 && nxc == 1
    F = C(1)*ones(ny,nx);
    return
end

if nxc == 1
    yc = linspace(1,ny,nyc);
    yf = 1:ny;
    F  = interp1(yc,C(:,1),yf,'linear','extrap').';
    F  = repmat(F,1,nx);
    return
end

if nyc == 1
    xc = linspace(1,nx,nxc);
    xf = 1:nx;
    row = interp1(xc,C(1,:),xf,'linear','extrap');
    F   = repmat(row,ny,1);
    return
end

[xc,yc] = meshgrid(linspace(1,nx,nxc),linspace(1,ny,nyc));
[xf,yf] = meshgrid(1:nx,1:ny);
F       = interp2(xc,yc,C,xf,yf,'linear');

if any(isnan(F(:)))
    Fn = interp2(xc,yc,C,xf,yf,'nearest');
    F(isnan(F)) = Fn(isnan(F));
end

end


function chi_imp = ConvexifyChi_ForImplicit_Local(chi_raw,PARAM)

Ne = size(chi_raw,1);
[ny,nx] = size(chi_raw{1,1});
N = ny*nx;

chi_floor = 1e-8;
chi_cap   = inf;

if isfield(PARAM,'CS_chi_floor')
    chi_floor = PARAM.CS_chi_floor;
end
if isfield(PARAM,'CS_chi_cap')
    chi_cap = PARAM.CS_chi_cap;
end

Cpage = zeros(Ne,Ne,N);

for i = 1:Ne
    for j = 1:Ne
        Cpage(i,j,:) = reshape(chi_raw{i,j},1,1,N);
    end
end

Cpage = 0.5*(Cpage + permute(Cpage,[2 1 3]));

[V,D] = pageeig(Cpage);

lam = zeros(Ne,N);

for i = 1:Ne
    lam(i,:) = reshape(D(i,i,:),1,N);
end

scale = max(1,max(abs(lam),[],1));
lam_imp = max(abs(lam),chi_floor.*scale);

if isfinite(chi_cap)
    lam_imp = min(lam_imp,chi_cap);
end

Dimp = zeros(Ne,Ne,N);

for i = 1:Ne
    Dimp(i,i,:) = reshape(lam_imp(i,:),1,1,N);
end

Cimp_page = pagemtimes(pagemtimes(V,Dimp),permute(V,[2 1 3]));
Cimp_page = 0.5*(Cimp_page + permute(Cimp_page,[2 1 3]));

chi_imp = cell(Ne,Ne);

for i = 1:Ne
    for j = 1:Ne
        chi_imp{i,j} = reshape(Cimp_page(i,j,:),ny,nx);
    end
end

end


function E = EnforceMeanE_Local(E,E_old)

Ne = numel(E);

for ie = 1:Ne
    target_mean = mean(E_old{ie}(:));
    new_mean    = mean(E{ie}(:));
    E{ie}       = E{ie} + target_mean - new_mean;
end

end
