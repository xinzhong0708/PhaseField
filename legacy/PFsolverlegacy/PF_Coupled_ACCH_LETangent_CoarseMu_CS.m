function [STATE_NEW,DIAG] = PF_Coupled_ACCH_LETangent_CoarseMu_CS(STATE_REF,PARAM,MODEL,GRID,PHYS,NUM,STATE_COEF)
%PF_COUPLED_ACCH_LETANGENT_COARSEMUM_CS
%
% Safer coarse-mu predictor.
%
% This version does NOT compute coarse dmu from a discarded coarse dphi.
% Instead:
%   1. solve fine-grid AC first with dmu = 0;
%   2. recover the fine E change from e*dp;
%   3. restrict this fine predicted state to a coarse grid;
%   4. solve a coarse fixed-p CH correction;
%   5. prolong coarse dmu back to the fine grid;
%   6. recover fine E = E_ref + e*dp + chi*dmu.
%
% This is more stable than using a coarse coupled ACCH solve for dmu,
% because the coarse chemical solve now responds to the actual fine-grid
% phase change.
%
% Controls:
%   NUM.CH_coarse_ratio      = 2 or 4. If <=1, calls original fine solver.
%   NUM.CM_coarse_kappac     = 0/1, default 0.
%   NUM.CM_dmu_theta         = under-relax dmu, default 1.
%   NUM.CM_dphi_cap          = cap max |dphi|, optional.
%   NUM.CM_dE_cap            = cap max |E_new-E_ref| by scaling dmu/dphi, optional.
%   NUM.CM_mean_correct      = 0/1, default 0.

if nargin < 7 || isempty(STATE_COEF)
    STATE_COEF = STATE_REF;
end

if ~isfield(NUM,'CH_coarse_ratio') || NUM.CH_coarse_ratio <= 1
    [STATE_NEW,DIAG] = PF_Coupled_ACCH_LETangent_CS(STATE_REF,PARAM,MODEL,GRID,PHYS,NUM,STATE_COEF);
    return
end

r = NUM.CH_coarse_ratio;

coarse_kappac = 0;
if isfield(NUM,'CM_coarse_kappac')
    coarse_kappac = NUM.CM_coarse_kappac;
end

mean_correct = 0;
if isfield(NUM,'CM_mean_correct')
    mean_correct = NUM.CM_mean_correct;
end

dmu_theta = 1.0;
if isfield(NUM,'CM_dmu_theta') && ~isempty(NUM.CM_dmu_theta)
    dmu_theta = NUM.CM_dmu_theta;
end

Ne = numel(STATE_REF.E);
[ny,nx,Ngrain] = size(STATE_REF.p);

% ------------------------------------------------------------
% 1. Fine AC predictor with dmu = 0.
% ------------------------------------------------------------
dmu_zero = cell(1,Ne);
for ie = 1:Ne
    dmu_zero{ie} = zeros(ny,nx);
end

[dphi,DIAG_AC] = Solve_Fine_AC_With_Dmu(STATE_REF,STATE_COEF,dmu_zero,PARAM,MODEL,GRID,PHYS,NUM);

% Optional AC increment cap. This is only a safety limiter for the coarse
% split predictor; the timestep controller still decides acceptance.
if isfield(NUM,'CM_dphi_cap') && ~isempty(NUM.CM_dphi_cap)
    mdp = max(abs(dphi(:)));
    if mdp > NUM.CM_dphi_cap
        dphi = dphi * (NUM.CM_dphi_cap/mdp);
    end
end

STATE_PHI       = STATE_REF;
STATE_PHI.phi   = STATE_REF.phi + dphi;

if isfield(NUM,'norm_phi') && NUM.norm_phi == 1
    STATE_PHI.phi = Normalize_Phi_L2_Local(STATE_PHI.phi);
end

if isfield(NUM,'cut_phi') && NUM.cut_phi == 1
    STATE_PHI.phi = Cut_Phi_Local(STATE_PHI.phi);
end

STATE_PHI.p     = Calc_p(MODEL,STATE_PHI.phi);
STATE_PHI.mu_e  = STATE_REF.mu_e;

% Recover the E change caused by fine-grid phase motion only.
e0 = STATE_COEF.e;
STATE_PHI.E = STATE_REF.E;

for ie = 1:Ne
    En = STATE_REF.E{ie};
    for ig = 1:Ngrain
        dp_ig = STATE_PHI.p(:,:,ig) - STATE_REF.p(:,:,ig);
        En = En + e0{ig}{ie} .* dp_ig;
    end
    STATE_PHI.E{ie} = En;
end

STATE_PHI.chi = STATE_COEF.chi;
STATE_PHI.e   = STATE_COEF.e;
STATE_PHI.c   = STATE_COEF.c;
STATE_PHI.omg = STATE_COEF.omg;
if isfield(STATE_COEF,'c_ext')
    STATE_PHI.c_ext = STATE_COEF.c_ext;
end

% ------------------------------------------------------------
% 2. Coarse fixed-p CH correction responding to the fine AC source.
% ------------------------------------------------------------
[STATE_REF_C,GRID_C] = Restrict_State_CH(STATE_REF,GRID,r);
[STATE_PHI_C,~]      = Restrict_State_CH(STATE_PHI,GRID,r);
PARAM_C              = Restrict_PARAM_CH(PARAM,r);

if coarse_kappac == 0 && isfield(PARAM_C,'use_kappa_c')
    PARAM_C.use_kappa_c = 0;
end

NUM_C = NUM;
NUM_C.CH_coarse_ratio = 1;
NUM_C.norm_E          = 0;
NUM_C.CHLE_force_full = 1;   % coarse solve should be global smooth background

[STATE_C_CORR,DIAG_C] = PF_CH_LECorrector_FixedP_Band_CS( ...
    STATE_REF_C,STATE_PHI_C,PARAM_C,MODEL,GRID_C,PHYS,NUM_C);

% Coarse dmu only.
dmu_C = cell(1,Ne);
for ie = 1:Ne
    dmu_C{ie} = STATE_C_CORR.mu_e{ie} - STATE_PHI_C.mu_e{ie};
end

% ------------------------------------------------------------
% 3. Prolong coarse dmu correction to fine grid.
% ------------------------------------------------------------
dmu_F = cell(1,Ne);
for ie = 1:Ne
    dmu_F{ie} = dmu_theta * Prolong2D(dmu_C{ie},GRID_C,GRID,[ny nx]);
end

% ------------------------------------------------------------
% 4. Build predicted fine state.
% ------------------------------------------------------------
chi0 = STATE_COEF.chi;

if isfield(PARAM,'use_CS_chi') && PARAM.use_CS_chi == 1
    chi_imp = ConvexifyChi_ForImplicit(chi0,PARAM);
else
    chi_imp = chi0;
end

STATE_NEW      = STATE_PHI;
STATE_NEW.mu_e = STATE_REF.mu_e;
for ie = 1:Ne
    STATE_NEW.mu_e{ie} = STATE_REF.mu_e{ie} + dmu_F{ie};
end

STATE_NEW.E = STATE_PHI.E;
for ie = 1:Ne
    En = STATE_PHI.E{ie};
    for je = 1:Ne
        En = En + chi_imp{ie,je} .* dmu_F{je};
    end
    STATE_NEW.E{ie} = En;
end

% Optional safety cap on the total E increment from this split predictor.
% Scale both dphi and dmu together and rebuild if the predictor jump is too large.
if isfield(NUM,'CM_dE_cap') && ~isempty(NUM.CM_dE_cap)
    max_dE = 0;
    for ie = 1:Ne
        max_dE = max(max_dE,max(abs(STATE_NEW.E{ie}(:)-STATE_REF.E{ie}(:))));
    end

    if max_dE > NUM.CM_dE_cap
        lam = NUM.CM_dE_cap/max_dE;
        dphi = lam*dphi;
        for ie = 1:Ne
            dmu_F{ie} = lam*dmu_F{ie};
        end

        STATE_PHI       = STATE_REF;
        STATE_PHI.phi   = STATE_REF.phi + dphi;
        if isfield(NUM,'norm_phi') && NUM.norm_phi == 1
            STATE_PHI.phi = Normalize_Phi_L2_Local(STATE_PHI.phi);
        end
        if isfield(NUM,'cut_phi') && NUM.cut_phi == 1
            STATE_PHI.phi = Cut_Phi_Local(STATE_PHI.phi);
        end
        STATE_PHI.p = Calc_p(MODEL,STATE_PHI.phi);
        STATE_PHI.E = STATE_REF.E;
        for ie = 1:Ne
            En = STATE_REF.E{ie};
            for ig = 1:Ngrain
                dp_ig = STATE_PHI.p(:,:,ig) - STATE_REF.p(:,:,ig);
                En = En + e0{ig}{ie} .* dp_ig;
            end
            STATE_PHI.E{ie} = En;
        end

        STATE_NEW       = STATE_PHI;
        STATE_NEW.mu_e  = STATE_REF.mu_e;
        for ie = 1:Ne
            STATE_NEW.mu_e{ie} = STATE_REF.mu_e{ie} + dmu_F{ie};
        end
        STATE_NEW.E = STATE_PHI.E;
        for ie = 1:Ne
            En = STATE_PHI.E{ie};
            for je = 1:Ne
                En = En + chi_imp{ie,je} .* dmu_F{je};
            end
            STATE_NEW.E{ie} = En;
        end
        STATE_NEW.chi = chi0;
        STATE_NEW.e   = e0;
        STATE_NEW.c   = STATE_COEF.c;
        STATE_NEW.omg = STATE_COEF.omg;
        if isfield(STATE_COEF,'c_ext')
            STATE_NEW.c_ext = STATE_COEF.c_ext;
        end
    end
end

if mean_correct == 1
    STATE_NEW.E = Match_Mean_E(STATE_NEW.E,STATE_REF.E);
end

% Predict omega for diagnostics only. The next LE_Run overwrites it.
STATE_NEW.omg = STATE_REF.omg;
for ig = 1:Ngrain
    domega = zeros(ny,nx);
    for ie = 1:Ne
        domega = domega - e0{ig}{ie} .* dmu_F{ie};
    end
    STATE_NEW.omg(:,:,ig) = STATE_REF.omg(:,:,ig) + domega;
end
omg_mean = mean(STATE_NEW.omg,3);
for ig = 1:Ngrain
    STATE_NEW.omg(:,:,ig) = STATE_NEW.omg(:,:,ig) - omg_mean;
end

STATE_NEW.chi = chi0;
STATE_NEW.e   = e0;

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
DIAG = struct();
DIAG.coarse          = DIAG_C;
DIAG.fine_ac         = DIAG_AC;
DIAG.ratio           = r;
DIAG.coarse_grid     = [GRID_C.ny GRID_C.nx];
DIAG.fine_grid       = [GRID.ny GRID.nx];
DIAG.max_dphi        = max(abs(dphi(:)));
DIAG.dmu_theta       = dmu_theta;

max_dmu = 0;
max_dE  = 0;
for ie = 1:Ne
    max_dmu = max(max_dmu,max(abs(dmu_F{ie}(:))));
    max_dE  = max(max_dE,max(abs(STATE_NEW.E{ie}(:)-STATE_REF.E{ie}(:))));
end
DIAG.max_dmu = max_dmu;
DIAG.max_dE  = max_dE;

end


%% ========================================================================
% Fine AC solve with known dmu
% ========================================================================
function [dphi,DIAG] = Solve_Fine_AC_With_Dmu(STATE_REF,STATE_COEF,dmu,PARAM,MODEL,GRID,PHYS,NUM)

phi_ref = STATE_REF.phi;
[ny,nx,Ngrain] = size(phi_ref);
Ne = numel(STATE_REF.E);
Nnode = ny*nx;

dt = NUM.dt_phy;
dx = GRID.dx;
dy = GRID.dy;
dx2 = dx^2;
dy2 = dy^2;

% Build AC source using old geometry but coefficient omega.
STATE_SRC       = STATE_REF;
STATE_SRC.mu_e  = STATE_COEF.mu_e;
STATE_SRC.chi   = STATE_COEF.chi;
STATE_SRC.e     = STATE_COEF.e;
STATE_SRC.omg   = STATE_COEF.omg;

if isfield(NUM,'use_Aac') && NUM.use_Aac == 1
    if isfield(NUM,'Aac_fac')
        Aac_fac = NUM.Aac_fac;
    else
        Aac_fac = 3;
    end
    PARAM.A_ac = Calc_Aac_FrozenOmega(STATE_SRC,PARAM,MODEL,Aac_fac,1e-6,0,[]);
else
    PARAM.A_ac = zeros(ny,nx);
end

STATE_SRC = Calc_S_AllenCahn(STATE_SRC,PARAM,MODEL);
S_AC      = STATE_SRC.S_AC;

% Interface-active phi mask, copied in spirit from the full coupled solver.
if isfield(NUM,'phi_mask_cut')
    phi_cut = NUM.phi_mask_cut;
else
    phi_cut = 1e-8;
end

if isfield(NUM,'phi_mask_pure_cut')
    pure_cut = NUM.phi_mask_pure_cut;
else
    pure_cut = phi_cut;
end

if isfield(NUM,'phi_mask_thick')
    mask_thick = NUM.phi_mask_thick;
else
    if (isfield(PHYS,'kap') && PHYS.kap ~= 0) || ...
            (isfield(PHYS,'kappa') && PHYS.kappa ~= 0)
        mask_thick = 2;
    else
        mask_thick = 1;
    end
end

maskPhi = Local_Calc_Interface_Mask(phi_ref,phi_cut,pure_cut,mask_thick);

if isfield(NUM,'phi_mask_source_tol') && ~isempty(NUM.phi_mask_source_tol)

    source_tol = NUM.phi_mask_source_tol;

    for alpha = 1:Ngrain
        core_source = abs(S_AC{alpha}) > source_tol;
        mask_source = Local_Dilate_Mask(core_source,mask_thick);
        maskPhi(:,:,alpha) = maskPhi(:,:,alpha) | mask_source;
    end

end

maskPhi = logical(maskPhi);

% Fast tangent coefficient B = sum_i e_i dp_i/dphi_alpha.
e0 = STATE_COEF.e;
eps_phi = 1e-14;
Dphi    = sum(phi_ref.^2,3) + eps_phi;
p_tan   = phi_ref.^2 ./ Dphi;

facPhi = cell(1,Ngrain);
for alpha = 1:Ngrain
    facPhi{alpha} = 2 .* phi_ref(:,:,alpha) ./ Dphi;
end

e_bar = cell(1,Ne);
for ie = 1:Ne
    tmp = zeros(ny,nx);
    for ig = 1:Ngrain
        tmp = tmp + p_tan(:,:,ig) .* e0{ig}{ie};
    end
    e_bar{ie} = tmp;
end

% Reflective neighbour indices.
[Igrid,Jgrid] = ndgrid(1:ny,1:nx);
ii = Igrid(:);
jj = Jgrid(:);

refI = @(i,sh) reflect_index(i + sh, ny);
refJ = @(j,sh) reflect_index(j + sh, nx);

idx_c = sub2ind([ny,nx],ii,jj);
idx_L = sub2ind([ny,nx],ii,refJ(jj,-1));
idx_R = sub2ind([ny,nx],ii,refJ(jj,+1));
idx_U = sub2ind([ny,nx],refI(ii,-1),jj);
idx_D = sub2ind([ny,nx],refI(ii,+1),jj);

% Optional diagonal source tangent.
if isfield(NUM,'use_Jphi') && NUM.use_Jphi == 1
    Jphi_diag = Calc_Jphi_Diag_FD(STATE_SRC,PARAM,MODEL,NUM,maskPhi);
else
    Jphi_diag = [];
end

% Solve one scalar implicit AC problem for each grain.
dphi = zeros(ny,nx,Ngrain);
Nphi = 0;

for alpha = 1:Ngrain

    mask_a = maskPhi(:,:,alpha);
    idx_a  = find(mask_a);

    if isempty(idx_a)
        continue
    end

    idmap = zeros(ny,nx);
    idmap(idx_a) = 1:numel(idx_a);

    row = idmap(idx_a);

    lap_phi_ref = laplacian_reflect(phi_ref(:,:,alpha),dx,dy);

    % Original full coupled AC row has coeff_mu*dmu on the left, where
    % coeff_mu = -L * B.  With known dmu, move it to the RHS:
    %   RHS = LK*lap(phi) + S_AC - coeff_mu*dmu
    rhs_full = PARAM.LK .* lap_phi_ref + S_AC{alpha};

    for ie = 1:Ne
        B_ie_alpha = facPhi{alpha} .* (e0{alpha}{ie} - e_bar{ie});
        coeff_mu   = -PARAM.L .* B_ie_alpha;
        rhs_full   = rhs_full - coeff_mu .* dmu{ie};
    end

    R = rhs_full(idx_a);

    LKc = PARAM.LK(idx_a);
    Aac = PARAM.A_ac(idx_a);

    cC = 1/dt + Aac + 2*LKc/dx2 + 2*LKc/dy2;
    cL = -LKc/dx2;
    cR = -LKc/dx2;
    cU = -LKc/dy2;
    cD = -LKc/dy2;

    if ~isempty(Jphi_diag)
        Jaa = Jphi_diag{alpha};
        cC = cC - Jaa(idx_a);
    end

    rows = zeros(numel(idx_a)*5,1);
    cols = zeros(numel(idx_a)*5,1);
    vals = zeros(numel(idx_a)*5,1);
    k = 1;

    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_a,       cC);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_L(idx_a),cL);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_R(idx_a),cR);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_U(idx_a),cU);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_D(idx_a),cD);

    A = sparse(rows(1:k-1),cols(1:k-1),vals(1:k-1),numel(idx_a),numel(idx_a));

    sol = A\R;

    tmp = zeros(ny,nx);
    tmp(idx_a) = sol;
    dphi(:,:,alpha) = tmp;

    Nphi = Nphi + numel(idx_a);

end

DIAG = struct();
DIAG.Nphi_active = Nphi;
DIAG.Nphi_full   = Ngrain*Nnode;
DIAG.active_phi_ratio = Nphi/max(Ngrain*Nnode,1);
DIAG.max_dphi = max(abs(dphi(:)));
DIAG.maskPhi = maskPhi;

end


%% ========================================================================
% Restriction / prolongation helpers
% ========================================================================
function [STATE_C,GRID_C] = Restrict_State_CH(STATE,GRID,r)

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
STATE_C.phi = BlockMean3D(STATE.phi,r);
STATE_C.p   = BlockMean3D(STATE.p,r);

if isfield(STATE,'omg') && ~isempty(STATE.omg)
    STATE_C.omg = BlockMean3D(STATE.omg,r);
end

if isfield(STATE,'mask') && ~isempty(STATE.mask)
    STATE_C.mask = BlockMean3D(STATE.mask,r) > 0;
end

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


function PARAM_C = Restrict_PARAM_CH(PARAM,r)

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

if isfield(PARAM,'L') && ~isempty(PARAM.L) && ~isscalar(PARAM.L)
    PARAM_C.L = BlockMean2D(PARAM.L,r);
end
if isfield(PARAM,'Lm') && ~isempty(PARAM.Lm) && ~isscalar(PARAM.Lm)
    PARAM_C.Lm = BlockMean2D(PARAM.Lm,r);
end
if isfield(PARAM,'LK') && ~isempty(PARAM.LK) && ~isscalar(PARAM.LK)
    PARAM_C.LK = BlockMean2D(PARAM.LK,r);
end
if isfield(PARAM,'LL') && ~isempty(PARAM.LL) && ~isscalar(PARAM.LL)
    PARAM_C.LL = BlockMean2D(PARAM.LL,r);
end
if isfield(PARAM,'LM') && ~isempty(PARAM.LM) && ~isscalar(PARAM.LM)
    PARAM_C.LM = BlockMean2D(PARAM.LM,r);
end

end


function B = BlockMean2D(A,r)

if isscalar(A)
    B = A;
    return
end

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

if any(isnan(Af(:)))
    Af2 = interp2(Xc,Yc,Ac,Xf,Yf,'nearest');
    Af(isnan(Af)) = Af2(isnan(Af));
end

end


function Af = Prolong3D(Ac,GRID_C,GRID_F,out_size)

Af = zeros(out_size);
for iz = 1:out_size(3)
    Af(:,:,iz) = Prolong2D(Ac(:,:,iz),GRID_C,GRID_F,out_size(1:2));
end

end


function E = Match_Mean_E(E,Eref)

for ie = 1:numel(E)
    E{ie} = E{ie} + mean(Eref{ie}(:)) - mean(E{ie}(:));
end

end


%% ========================================================================
% Small local numerical helpers
% ========================================================================
function idx = reflect_index(idx,n)

if n == 1
    idx = ones(size(idx));
    return
end

period = 2*n - 2;
r = mod(idx - 1,period);
idx = 1 + min(r,period - r);

end


function L = laplacian_reflect(A,dx,dy)

[ny,nx] = size(A);

if nx == 1
    AL = A;
    AR = A;
else
    AL = A(:,[2,1:nx-1]);
    AR = A(:,[2:nx,nx-1]);
end

if ny == 1
    AU = A;
    AD = A;
else
    AU = A([2,1:ny-1],:);
    AD = A([2:ny,ny-1],:);
end

L = (AL - 2*A + AR)/dx^2 + (AU - 2*A + AD)/dy^2;

end


function [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_col,coeff)

col = idmap(idx_col);

row   = row(:);
col   = col(:);
coeff = coeff(:);

keep = col > 0 & isfinite(coeff) & coeff ~= 0;

if ~any(keep)
    return
end

r = row(keep);
c = col(keep);
v = coeff(keep);

n = numel(r);

if k+n-1 > numel(rows)
    grow = max(numel(rows),n + 1000);
    rows = [rows; zeros(grow,1)];
    cols = [cols; zeros(grow,1)];
    vals = [vals; zeros(grow,1)];
end

rows(k:k+n-1) = r(:);
cols(k:k+n-1) = c(:);
vals(k:k+n-1) = v(:);

k = k + n;

end


function phi = Cut_Phi_Local(phi)

phi(phi < 0) = 0;
phi(phi > 1) = 1;

end


function phi = Normalize_Phi_L2_Local(phi)

phi = max(phi,0);
phi = min(phi,1);

s = sqrt(sum(phi.^2,3));
mask = s > eps;

for ip = 1:size(phi,3)
    tmp = phi(:,:,ip);
    tmp(mask) = tmp(mask) ./ s(mask);
    phi(:,:,ip) = tmp;
end

end


function mask = Local_Calc_Interface_Mask(phi,low_cut,pure_cut,thickness)

[ny,nx,Ngrain] = size(phi);
mask = false(ny,nx,Ngrain);

den = sum(phi.^2,3) + eps;

for alpha = 1:Ngrain

    q = phi(:,:,alpha).^2 ./ den;

    core = q > low_cut & q < 1 - pure_cut;

    if nx == 1
        qL = q;
        qR = q;
    else
        qL = q(:,[2,1:nx-1]);
        qR = q(:,[2:nx,nx-1]);
    end

    if ny == 1
        qU = q;
        qD = q;
    else
        qU = q([2,1:ny-1],:);
        qD = q([2:ny,ny-1],:);
    end

    jump = abs(q - qL) > low_cut | ...
           abs(q - qR) > low_cut | ...
           abs(q - qU) > low_cut | ...
           abs(q - qD) > low_cut;

    core = core | jump;

    mask(:,:,alpha) = Local_Dilate_Mask(core,thickness);

end

end


function mask = Local_Dilate_Mask(core,thickness)

if thickness <= 0
    mask = core;
    return
end

ker  = ones(2*thickness+1,2*thickness+1);
mask = conv2(double(core),ker,'same') > 0;

end


function Jdiag = Calc_Jphi_Diag_FD(STATE,PARAM,MODEL,NUM,maskPhi)

phi0 = STATE.phi;
[ny,nx,Ngrain] = size(phi0);

if isfield(NUM,'Jphi_eps') && ~isempty(NUM.Jphi_eps)
    eps_fd = NUM.Jphi_eps;
else
    eps_fd = 1e-6;
end

STATE0 = Calc_S_AllenCahn(STATE,PARAM,MODEL);
S0     = STATE0.S_AC;

Jdiag = cell(1,Ngrain);

for alpha = 1:Ngrain

    Jdiag{alpha} = zeros(ny,nx);
    mask_a = maskPhi(:,:,alpha);

    if ~any(mask_a(:))
        continue
    end

    STATEP = STATE;
    STATEP.phi = phi0;

    tmp = STATEP.phi(:,:,alpha);
    tmp(mask_a) = tmp(mask_a) + eps_fd;
    STATEP.phi(:,:,alpha) = tmp;

    STATEP = Calc_S_AllenCahn(STATEP,PARAM,MODEL);
    Sp = STATEP.S_AC;

    Jtmp = (Sp{alpha} - S0{alpha}) / eps_fd;
    Jdiag{alpha}(mask_a) = Jtmp(mask_a);

end

end


function chi_imp = ConvexifyChi_ForImplicit(chi_raw,PARAM)

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
