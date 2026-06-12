function [STATE_CORR,DIAG] = PF_CH_LECorrector_FixedP_CoarseFine_Band_CS(STATE_TIME,STATE_IT,PARAM,MODEL,GRID,PHYS,NUM)
%PF_CH_LECORRECTOR_FIXEDP_COARSEFINE_BAND_CS Two-level fixed-p CH corrector.
%
% Idea:
%   1) Restrict the fixed-p CH correction problem to a coarse grid.
%   2) Solve global CH on the coarse grid.
%   3) Prolong the coarse dmu correction back to the fine grid.
%   4) Apply this as a smooth global background correction on the fine grid.
%   5) Optionally call the existing fine-grid band corrector for local error.
%
% This function does not move phi.  It is only a CH chemical correction.
%
% Required extra NUM fields:
%   NUM.CH_coarse_ratio       coarsening factor, e.g. 2, 4
%
% Optional:
%   NUM.CH_do_fine_corrector  1/0, default 1
%   NUM.CH_coarse_kappac      1/0, default 0  (recommended 0)
%   NUM.CH_mean_correct       1/0, default 1

if ~isfield(NUM,'CH_coarse_ratio') || NUM.CH_coarse_ratio <= 1
    [STATE_CORR,DIAG] = PF_CH_LECorrector_FixedP_Band_CS(STATE_TIME,STATE_IT,PARAM,MODEL,GRID,PHYS,NUM);
    return
end

r = NUM.CH_coarse_ratio;

if ~isfield(NUM,'CH_do_fine_corrector')
    do_fine = 1;
else
    do_fine = NUM.CH_do_fine_corrector;
end

if ~isfield(NUM,'CH_coarse_kappac')
    coarse_kappac = 0;
else
    coarse_kappac = NUM.CH_coarse_kappac;
end

if ~isfield(NUM,'CH_mean_correct')
    mean_correct = 1;
else
    mean_correct = NUM.CH_mean_correct;
end

% ------------------------------------------------------------
% Build coarse fixed-time and iterate states
% ------------------------------------------------------------
[STATE_TIME_C,GRID_C] = Restrict_State_CH(STATE_TIME,GRID,r);
[STATE_IT_C,~]        = Restrict_State_CH(STATE_IT,GRID,r);

PARAM_C = Restrict_PARAM_CH(PARAM,GRID,GRID_C,r);

% The coarse solve should give only the smooth global chemical background.
% Let the fine local solve handle composition-kappa and spinodal/interface details.
if coarse_kappac == 0 && isfield(PARAM_C,'use_kappa_c')
    PARAM_C.use_kappa_c = 0;
end

NUM_C = NUM;
NUM_C.CHLE_force_full = 1;      % global on the coarse grid
if isfield(NUM,'CH_band_width') && ~isempty(NUM.CH_band_width)
    NUM_C.CHLE_band_thick = NUM.CH_band_width;
end
NUM_C.norm_E          = 0;      % mean handled on fine grid if requested

% ------------------------------------------------------------
% Coarse global CH solve
% ------------------------------------------------------------
[STATE_C_CORR,DIAG_C] = PF_CH_LECorrector_FixedP_BandLocal_CS( ...
    STATE_TIME_C,STATE_IT_C,PARAM_C,MODEL,GRID_C,PHYS,NUM_C);

% Coarse correction only, not absolute replacement.
dmu_C = cell(1,numel(STATE_IT.mu_e));
for ie = 1:numel(dmu_C)
    dmu_C{ie} = STATE_C_CORR.mu_e{ie} - STATE_IT_C.mu_e{ie};
end

% ------------------------------------------------------------
% Prolong coarse dmu correction to fine grid
% ------------------------------------------------------------
dmu_F = cell(1,numel(dmu_C));
for ie = 1:numel(dmu_C)
    dmu_F{ie} = Prolong2D(dmu_C{ie},GRID_C,GRID,size(STATE_IT.E{ie}));
end

% ------------------------------------------------------------
% Apply smooth global background correction to fine state at fixed p
% ------------------------------------------------------------
STATE_BG = Apply_Dmu_FixedP(STATE_IT,dmu_F);

if mean_correct == 1
    STATE_BG.E = Match_Mean_E(STATE_BG.E,STATE_TIME.E);
end

% ------------------------------------------------------------
% Fine local correction around interface/spinodal band
% ------------------------------------------------------------
if do_fine == 1
    NUM_F = NUM;
    NUM_F.CHLE_force_full = 0;
    if isfield(NUM,'CH_band_width') && ~isempty(NUM.CH_band_width)
        NUM_F.CHLE_band_thick = NUM.CH_band_width;
    end

    PARAM_F = PARAM;

    % Default: local fine repair uses standard CH operator.
    % This avoids the old behavior where use_kappa_c forced full-domain fine solve.
    if ~isfield(NUM,'CH_fine_use_kappa_c') || NUM.CH_fine_use_kappa_c == 0
        PARAM_F.use_kappa_c = 0;
        NUM_F.CH_fast_disable_kappac = 1;
    end

    % For speed, the fast local corrector can disable composition-space
    % kappa in the fine repair. The coarse global solve already supplied
    % the smooth background. The local repair then fixes the main interface
    % chemical residual using the standard CH operator.
    if isfield(NUM,'CH_fast_disable_kappac') && NUM.CH_fast_disable_kappac == 1
        PARAM_F.use_kappa_c = 0;
    end

    if isfield(NUM,'CH_use_fast_local') && NUM.CH_use_fast_local == 0
        [STATE_CORR,DIAG_F] = PF_CH_LECorrector_FixedP_Band_CS( ...
            STATE_TIME,STATE_BG,PARAM,MODEL,GRID,PHYS,NUM_F);
    else
        [STATE_CORR,DIAG_F] = PF_CH_LECorrector_FixedP_BandLocal_CS( ...
            STATE_TIME,STATE_BG,PARAM_F,MODEL,GRID,PHYS,NUM_F);
    end
else
    STATE_CORR = STATE_BG;
    DIAG_F = struct();
end

DIAG = struct();
DIAG.coarse = DIAG_C;
DIAG.fine   = DIAG_F;
DIAG.ratio  = r;
DIAG.coarse_grid = [GRID_C.ny GRID_C.nx];
DIAG.fine_grid   = [GRID.ny GRID.nx];

end


function [STATE_C,GRID_C] = Restrict_State_CH(STATE,GRID,r)
%RESTRICT_STATE_CH Block-average all CH-relevant fields.

[ny,nx,~] = size(STATE.p);
nyc = ceil(ny/r);
nxc = ceil(nx/r);

GRID_C       = GRID;
GRID_C.nx    = nxc;
GRID_C.ny    = nyc;
GRID_C.dx    = GRID.dx*r;
GRID_C.dy    = GRID.dy*r;
GRID_C.x     = linspace(GRID.x(1),GRID.x(end),nxc);
GRID_C.y     = linspace(GRID.y(1),GRID.y(end),nyc);
GRID_C.Lx    = GRID.Lx;
GRID_C.Ly    = GRID.Ly;

STATE_C = STATE;

% 3D grain fields
STATE_C.phi = BlockMean3D(STATE.phi,r);
STATE_C.p   = BlockMean3D(STATE.p,r);

if isfield(STATE,'omg') && ~isempty(STATE.omg)
    STATE_C.omg = BlockMean3D(STATE.omg,r);
end

if isfield(STATE,'mask') && ~isempty(STATE.mask)
    STATE_C.mask = BlockMean3D(STATE.mask,r) > 0;
end

% cell fields
STATE_C.E    = BlockMeanCell1(STATE.E,r);
STATE_C.mu_e = BlockMeanCell1(STATE.mu_e,r);

if isfield(STATE,'chi') && ~isempty(STATE.chi)
    STATE_C.chi = BlockMeanCell2(STATE.chi,r);
end

if isfield(STATE,'e') && ~isempty(STATE.e)
    STATE_C.e = BlockMeanNestedCell(STATE.e,r);
end

if isfield(STATE,'c') && ~isempty(STATE.c)
    STATE_C.c = BlockMeanNestedCell(STATE.c,r);
end

if isfield(STATE,'c_ext') && ~isempty(STATE.c_ext)
    STATE_C.c_ext = BlockMeanNestedCell(STATE.c_ext,r);
end

end


function PARAM_C = Restrict_PARAM_CH(PARAM,GRID,GRID_C,r)
%RESTRICT_PARAM_CH Restrict spatial CH coefficients in PARAM.

PARAM_C = PARAM;

if isfield(PARAM,'M') && ~isempty(PARAM.M)
    PARAM_C.M = BlockMeanCell1(PARAM.M,r);
end

if isfield(PARAM,'eta') && ~isempty(PARAM.eta)
    if isscalar(PARAM.eta)
        PARAM_C.eta = PARAM.eta;
    else
        PARAM_C.eta = BlockMean2D(PARAM.eta,r);
    end
end

if isfield(PARAM,'kappa_eff') && ~isempty(PARAM.kappa_eff)
    if isscalar(PARAM.kappa_eff)
        PARAM_C.kappa_eff = PARAM.kappa_eff;
    elseif ndims(PARAM.kappa_eff) == 2
        PARAM_C.kappa_eff = BlockMean2D(PARAM.kappa_eff,r);
    else
        PARAM_C.kappa_eff = BlockMean3D(PARAM.kappa_eff,r);
    end
end

% Coarse L fields are not essential for fixed-p CH, but keep them consistent.
if isfield(PARAM,'L') && ~isempty(PARAM.L) && ~isscalar(PARAM.L)
    PARAM_C.L = BlockMean2D(PARAM.L,r);
end
if isfield(PARAM,'Lm') && ~isempty(PARAM.Lm) && ~isscalar(PARAM.Lm)
    PARAM_C.Lm = BlockMean2D(PARAM.Lm,r);
end
if isfield(PARAM,'LK') && ~isempty(PARAM.LK) && ~isscalar(PARAM.LK)
    PARAM_C.LK = BlockMean2D(PARAM.LK,r);
end

% Make sure dimensions are sane if any function checks GRID indirectly.
PARAM_C.GRID_size = [GRID_C.ny GRID_C.nx]; %#ok<STRNU>

end


function STATE = Apply_Dmu_FixedP(STATE,dmu)
%APPLY_DMU_FIXEDP Apply dmu and recover E = E + chi*dmu at fixed p.

Ne = numel(STATE.E);

for ie = 1:Ne
    STATE.mu_e{ie} = STATE.mu_e{ie} + dmu{ie};
end

for ie = 1:Ne
    En = STATE.E{ie};
    for je = 1:Ne
        En = En + STATE.chi{ie,je}.*dmu{je};
    end
    STATE.E{ie} = En;
end

if isfield(STATE,'omg') && isfield(STATE,'e')
    for ig = 1:size(STATE.p,3)
        domega = zeros(size(STATE.E{1}));
        for ie = 1:Ne
            domega = domega - STATE.e{ig}{ie}.*dmu{ie};
        end
        STATE.omg(:,:,ig) = STATE.omg(:,:,ig) + domega;
    end
    omg_mean = mean(STATE.omg,3);
    for ig = 1:size(STATE.p,3)
        STATE.omg(:,:,ig) = STATE.omg(:,:,ig) - omg_mean;
    end
end

end


function E = Match_Mean_E(E,Eref)
%MATCH_MEAN_E Preserve global mean E of the fixed-time state.

for ie = 1:numel(E)
    E{ie} = E{ie} + mean(Eref{ie}(:)) - mean(E{ie}(:));
end

end


function B = BlockMean2D(A,r)
%BLOCKMEAN2D Simple block average with edge blocks allowed.

[ny,nx] = size(A);
nyc = ceil(ny/r);
nxc = ceil(nx/r);
B = zeros(nyc,nxc);

for iy = 1:nyc
    iy1 = (iy-1)*r + 1;
    iy2 = min(iy*r,ny);
    for ix = 1:nxc
        ix1 = (ix-1)*r + 1;
        ix2 = min(ix*r,nx);
        blk = A(iy1:iy2,ix1:ix2);
        B(iy,ix) = mean(blk(:));
    end
end

end


function B = BlockMean3D(A,r)

[~,~,nz] = size(A);
B0 = BlockMean2D(A(:,:,1),r);
B = zeros(size(B0,1),size(B0,2),nz);
B(:,:,1) = B0;
for iz = 2:nz
    B(:,:,iz) = BlockMean2D(A(:,:,iz),r);
end

end


function C = BlockMeanCell1(C0,r)

C = C0;
for i = 1:numel(C0)
    if isempty(C0{i}) || isscalar(C0{i})
        C{i} = C0{i};
    else
        C{i} = BlockMean2D(C0{i},r);
    end
end

end


function C = BlockMeanCell2(C0,r)

C = C0;
for i = 1:size(C0,1)
    for j = 1:size(C0,2)
        if isempty(C0{i,j}) || isscalar(C0{i,j})
            C{i,j} = C0{i,j};
        else
            C{i,j} = BlockMean2D(C0{i,j},r);
        end
    end
end

end


function C = BlockMeanNestedCell(C0,r)

C = C0;
for i = 1:numel(C0)
    for j = 1:numel(C0{i})
        if isempty(C0{i}{j}) || isscalar(C0{i}{j})
            C{i}{j} = C0{i}{j};
        else
            C{i}{j} = BlockMean2D(C0{i}{j},r);
        end
    end
end

end


function Af = Prolong2D(Ac,GRID_C,GRID_F,out_size)
%PROLONG2D Linear interpolation from coarse cell centers to fine grid.

nyf = out_size(1);
nxf = out_size(2);

if numel(GRID_C.x) == 1 && numel(GRID_C.y) == 1
    Af = Ac(1)*ones(nyf,nxf);
    return
end

if numel(GRID_C.y) == 1
    Af = interp1(GRID_C.x,Ac,GRID_F.x,'linear','extrap');
    Af = repmat(Af,nyf,1);
    return
end

if numel(GRID_C.x) == 1
    Af = interp1(GRID_C.y,Ac,GRID_F.y,'linear','extrap');
    Af = repmat(Af(:),1,nxf);
    return
end

[Xc,Yc] = meshgrid(GRID_C.x,GRID_C.y);
[Xf,Yf] = meshgrid(GRID_F.x,GRID_F.y);
Af = interp2(Xc,Yc,Ac,Xf,Yf,'linear');

% Fill possible NaN from numerical edge mismatch by nearest.
if any(isnan(Af(:)))
    Af2 = interp2(Xc,Yc,Ac,Xf,Yf,'nearest');
    Af(isnan(Af)) = Af2(isnan(Af));
end

end
