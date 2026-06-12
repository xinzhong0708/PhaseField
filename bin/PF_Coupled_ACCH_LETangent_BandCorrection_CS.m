function [STATE_NEW,DIAG] = PF_Coupled_ACCH_LETangent_BandCorrection_CS(STATE_REF,STATE_BG,PARAM,MODEL,GRID,PHYS,NUM,STATE_COEF)
%PF_COUPLED_ACCH_LETANGENT_BANDCORRECTION_CS
%
% Fine-grid local ACCH correction on top of a coarse-grid CH background.
%
% The total chemical update is
%
%   dmu_total = dmu_bg + dmu_loc
%
% where
%
%   dmu_bg  = STATE_BG.mu_e - STATE_REF.mu_e
%   dmu_loc = unknown only inside a dilated interface band
%
% Unknowns:
%   x = [active dphi_grain ; local-band dmu_loc]
%
% Outside the local band, dmu_loc = 0. Therefore the coarse-grid CH result
% remains the boundary/background chemical potential.
%
% Controls:
%   NUM.ACCH_band_width        default 12
%   NUM.ACCH_p_cut             default NUM.CHLE_p_cut or 1e-8
%   NUM.ACCH_include_kappa_mask default 1
%   NUM.ACCH_fine_use_kappa_c  default 0
%
% Note:
%   The local ACCH correction uses the standard E-space CH operator by
%   default. This is intentional because composition-space kappa on a local
%   correction band is much more delicate. The coarse CH background and the
%   final LE still keep the thermodynamic state consistent.

if nargin < 8 || isempty(STATE_COEF)
    STATE_COEF = STATE_BG;
end

% Optional cache for repeated sparse column ordering.
persistent ORDER_CACHE
if isfield(NUM,'clear_order_cache') && NUM.clear_order_cache == 1
    ORDER_CACHE = [];
end

% ------------------------------------------------------------
% Fixed time-reference state
% ------------------------------------------------------------
phi_ref = STATE_REF.phi;
p_ref   = STATE_REF.p;
E_ref   = STATE_REF.E;
mu_ref  = STATE_REF.mu_e;

% ------------------------------------------------------------
% Coefficient/background state
% ------------------------------------------------------------
chi0 = STATE_COEF.chi;
e0   = STATE_COEF.e;

if isfield(PARAM,'use_CS_chi') && PARAM.use_CS_chi == 1
    chi_imp = ConvexifyChi_ForImplicit(chi0,PARAM);
else
    chi_imp = chi0;
end

% ------------------------------------------------------------
% Sizes
% ------------------------------------------------------------
[ny,nx,Ngrain] = size(phi_ref);
Ne             = numel(E_ref);
Nnode          = nx*ny;

dt = NUM.dt_phy;
dx = GRID.dx;
dy = GRID.dy;

dx2 = dx^2;
dy2 = dy^2;
dx4 = dx2^2;
dy4 = dy2^2;

% ------------------------------------------------------------
% Spatially varying kappa for CH gradient stabilization
% ------------------------------------------------------------
if isfield(PARAM,'kappa_eff') && ~isempty(PARAM.kappa_eff)
    kappa_eff = PARAM.kappa_eff;
    if isscalar(kappa_eff)
        kappa_eff = kappa_eff * ones(ny,nx);
    end
else
    if isfield(PHYS,'kap') && ~isempty(PHYS.kap)
        kappa0 = PHYS.kap;
    elseif isfield(PHYS,'kappa') && ~isempty(PHYS.kappa)
        kappa0 = PHYS.kappa;
    else
        kappa0 = 0;
    end
    kappa_eff = kappa0 * ones(ny,nx);
end

% Local band correction uses standard CH operator unless explicitly enabled.
use_kappa_c = 0;
if isfield(NUM,'ACCH_fine_use_kappa_c') && NUM.ACCH_fine_use_kappa_c == 1
    use_kappa_c = isfield(PARAM,'use_kappa_c') && PARAM.use_kappa_c == 1;
end

if use_kappa_c
    error('PF_Coupled_ACCH_LETangent_BandCorrection_CS: local kappa_c is not implemented in this first band ACCH version. Set NUM.ACCH_fine_use_kappa_c = 0.')
end

% ------------------------------------------------------------
% Known coarse-grid CH background dmu
% ------------------------------------------------------------
dmu_bg = cell(1,Ne);

for ie = 1:Ne
    dmu_bg{ie} = STATE_BG.mu_e{ie} - STATE_REF.mu_e{ie};
end

% ------------------------------------------------------------
% Build AC source using old geometry but background/coefficient omega
% ------------------------------------------------------------
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

% ------------------------------------------------------------
% Grid indexing with reflective boundary
% ------------------------------------------------------------
[Igrid,Jgrid] = ndgrid(1:ny,1:nx);
ii = Igrid(:);
jj = Jgrid(:);

refI = @(i,sh) reflect_index(i + sh, ny);
refJ = @(j,sh) reflect_index(j + sh, nx);

jjL  = refJ(jj,-1);
jjR  = refJ(jj,+1);
iiU  = refI(ii,-1);
iiD  = refI(ii,+1);

jjL2 = refJ(jj,-2);
jjR2 = refJ(jj,+2);
iiU2 = refI(ii,-2);
iiD2 = refI(ii,+2);

iiUR = refI(ii,-1); jjUR = refJ(jj,+1);
iiDR = refI(ii,+1); jjDR = refJ(jj,+1);
iiUL = refI(ii,-1); jjUL = refJ(jj,-1);
iiDL = refI(ii,+1); jjDL = refJ(jj,-1);

idx_c  = sub2ind([ny,nx], ii,   jj);
idx_L  = sub2ind([ny,nx], ii,   jjL);
idx_R  = sub2ind([ny,nx], ii,   jjR);
idx_U  = sub2ind([ny,nx], iiU,  jj);
idx_D  = sub2ind([ny,nx], iiD,  jj);

idx_L2 = sub2ind([ny,nx], ii,   jjL2);
idx_R2 = sub2ind([ny,nx], ii,   jjR2);
idx_U2 = sub2ind([ny,nx], iiU2, jj);
idx_D2 = sub2ind([ny,nx], iiD2, jj);

idx_UR = sub2ind([ny,nx], iiUR, jjUR);
idx_DR = sub2ind([ny,nx], iiDR, jjDR);
idx_UL = sub2ind([ny,nx], iiUL, jjUL);
idx_DL = sub2ind([ny,nx], iiDL, jjDL);

% ------------------------------------------------------------
% Active phi mask: solve only around interfaces
% ------------------------------------------------------------
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
active_per_cell = sum(maskPhi,3);

% ------------------------------------------------------------
% Mu correction band
% ------------------------------------------------------------
if isfield(NUM,'ACCH_p_cut') && ~isempty(NUM.ACCH_p_cut)
    p_cut = NUM.ACCH_p_cut;
elseif isfield(NUM,'CHLE_p_cut') && ~isempty(NUM.CHLE_p_cut)
    p_cut = NUM.CHLE_p_cut;
else
    p_cut = 1e-8;
end

if isfield(NUM,'ACCH_band_width') && ~isempty(NUM.ACCH_band_width)
    band_width = NUM.ACCH_band_width;
elseif isfield(NUM,'CH_band_width') && ~isempty(NUM.CH_band_width)
    band_width = NUM.CH_band_width;
elseif isfield(NUM,'CHLE_band_thick') && ~isempty(NUM.CHLE_band_thick)
    band_width = NUM.CHLE_band_thick;
else
    band_width = 12;
end

if band_width < 2
    error('NUM.ACCH_band_width must be at least 2 because the CH operator contains a +/-2 stencil.')
end

mask_core = Local_Interface_Core_From_p(STATE_REF.p,p_cut);
mask_core = mask_core | any(maskPhi,3);

include_kappa = 1;
if isfield(NUM,'ACCH_include_kappa_mask')
    include_kappa = NUM.ACCH_include_kappa_mask;
elseif isfield(NUM,'CHLE_include_kappa_mask')
    include_kappa = NUM.CHLE_include_kappa_mask;
end

if include_kappa == 1
    mask_core = mask_core | (kappa_eff > 0);
end

maskMu = Local_Dilate_Mask(mask_core,band_width);

% ------------------------------------------------------------
% Optional real diagonal AC source Jacobian
% ------------------------------------------------------------
if isfield(NUM,'use_Jphi') && NUM.use_Jphi == 1
    Jphi_diag = Calc_Jphi_Diag_FD(STATE_SRC,PARAM,MODEL,NUM,maskPhi);
else
    Jphi_diag = [];
end

% ------------------------------------------------------------
% Unknown ids
% ------------------------------------------------------------
idPhiMap = cell(1,Ngrain);

counter = 0;

for alpha = 1:Ngrain

    idmap = zeros(ny,nx);
    ids   = find(maskPhi(:,:,alpha));
    nids  = numel(ids);

    idmap(ids) = counter + (1:nids);

    idPhiMap{alpha} = idmap;
    counter = counter + nids;

end

Nphi = counter;

idsMuNode = find(maskMu);
NmuNode   = numel(idsMuNode);
Nmu       = Ne*NmuNode;
Ntot      = Nphi + Nmu;

idMuMap = cell(1,Ne);

for ie = 1:Ne
    idmap = zeros(ny,nx);
    idmap(idsMuNode) = Nphi + (ie-1)*NmuNode + (1:NmuNode);
    idMuMap{ie} = idmap;
end

% ------------------------------------------------------------
% Fast grain-aware tangent ingredients
% ------------------------------------------------------------
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

% ------------------------------------------------------------
% Sparse matrix allocation
% ------------------------------------------------------------
max_nnz = ...
    Nphi * (6 + Ne) + ...
    NmuNode * Ne * (5 + 13*Ne + 13*Ngrain) + ...
    1000;

rows = zeros(max_nnz,1);
cols = zeros(max_nnz,1);
vals = zeros(max_nnz,1);
R    = zeros(Ntot,1);

k = 1;

% ============================================================
% 1. Allen-Cahn block
% ============================================================
for alpha = 1:Ngrain

    mask_a = maskPhi(:,:,alpha);
    idx_a  = find(mask_a);

    if isempty(idx_a)
        continue
    end

    row = idPhiMap{alpha}(idx_a);

    phi_a = phi_ref(:,:,alpha);

    lap_phi_ref = laplacian_reflect(phi_a,dx,dy);

    rhs_full = PARAM.LK .* lap_phi_ref + S_AC{alpha};

    LKc = PARAM.LK(idx_a);
    Aac = PARAM.A_ac(idx_a);

    cC = 1/dt + Aac + 2*LKc/dx2 + 2*LKc/dy2;
    cL = -LKc/dx2;
    cR = -LKc/dx2;
    cU = -LKc/dy2;
    cD = -LKc/dy2;

    idxa_L = idx_L(idx_a);
    idxa_R = idx_R(idx_a);
    idxa_U = idx_U(idx_a);
    idxa_D = idx_D(idx_a);

    R(row) = rhs_full(idx_a);

    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idPhiMap{alpha},idx_a,  cC);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idPhiMap{alpha},idxa_L,cL);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idPhiMap{alpha},idxa_R,cR);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idPhiMap{alpha},idxa_U,cU);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idPhiMap{alpha},idxa_D,cD);

    if ~isempty(Jphi_diag)
        Jaa = Jphi_diag{alpha};
        [rows,cols,vals,k] = add_active_block( ...
            rows,cols,vals,k,row,idPhiMap{alpha},idx_a,-Jaa(idx_a));
    end

    % Coupling to local dmu correction. The known background dmu is moved
    % to the RHS.
    for ie = 1:Ne

        B_ie_alpha = facPhi{alpha} .* (e0{alpha}{ie} - e_bar{ie});
        coeff_mu   = -PARAM.L .* B_ie_alpha;

        % S_AC is already evaluated with the background state.
        % Therefore only the unknown local correction dmu_loc enters here.
        [rows,cols,vals,k] = add_block( ...
            rows,cols,vals,k,row,idMuMap{ie}(idx_a),coeff_mu(idx_a));

    end

end

% ============================================================
% 2. Local Cahn-Hilliard rows for dmu_loc
% ============================================================
idx_m = idsMuNode(:);

for l = 1:Ne

    row = idMuMap{l}(idx_m);

    Ml   = PARAM.M{l};
    M_c  = Ml(idx_m);
    M_L  = Ml(idx_L(idx_m));
    M_R  = Ml(idx_R(idx_m));
    M_U  = Ml(idx_U(idx_m));
    M_D  = Ml(idx_D(idx_m));

    d_L = -(M_L + M_c)/2/dx2;
    d_R = -(M_R + M_c)/2/dx2;
    d_U = -(M_U + M_c)/2/dy2;
    d_D = -(M_D + M_c)/2/dy2;
    d_C = -(d_L + d_R + d_U + d_D);

    [q_C,q_L,q_R,q_U,q_D,q_L2,q_R2,q_U2,q_D2,q_UR,q_DR,q_UL,q_DL] = ...
        Kappa4_Coeffs(Ml.*kappa_eff,idx_m,idx_L(idx_m),idx_R(idx_m),idx_U(idx_m),idx_D(idx_m), ...
                      idx_L2(idx_m),idx_R2(idx_m),idx_U2(idx_m),idx_D2(idx_m), ...
                      idx_UR(idx_m),idx_DR(idx_m),idx_UL(idx_m),idx_DL(idx_m), ...
                      dx2,dy2,dx4,dy4);

    a_C  = 1/dt + q_C;
    a_L  = q_L;
    a_R  = q_R;
    a_U  = q_U;
    a_D  = q_D;

    a_L2 = q_L2;
    a_R2 = q_R2;
    a_U2 = q_U2;
    a_D2 = q_D2;

    a_UR = q_UR;
    a_DR = q_DR;
    a_UL = q_UL;
    a_DL = q_DL;

    E_l  = E_ref{l};
    mu_l = mu_ref{l};

    E_c  = E_l(idx_m);
    E_L  = E_l(idx_L(idx_m));
    E_R  = E_l(idx_R(idx_m));
    E_U  = E_l(idx_U(idx_m));
    E_D  = E_l(idx_D(idx_m));
    E_L2 = E_l(idx_L2(idx_m));
    E_R2 = E_l(idx_R2(idx_m));
    E_U2 = E_l(idx_U2(idx_m));
    E_D2 = E_l(idx_D2(idx_m));
    E_UR = E_l(idx_UR(idx_m));
    E_DR = E_l(idx_DR(idx_m));
    E_UL = E_l(idx_UL(idx_m));
    E_DL = E_l(idx_DL(idx_m));

    mu_c = mu_l(idx_m);
    mu_L = mu_l(idx_L(idx_m));
    mu_R = mu_l(idx_R(idx_m));
    mu_U = mu_l(idx_U(idx_m));
    mu_D = mu_l(idx_D(idx_m));

    AE_Eref = ...
        a_C  .* E_c  + ...
        a_L  .* E_L  + a_R  .* E_R  + a_U  .* E_U  + a_D  .* E_D  + ...
        a_L2 .* E_L2 + a_R2 .* E_R2 + a_U2 .* E_U2 + a_D2 .* E_D2 + ...
        a_UR .* E_UR + a_DR .* E_DR + a_UL .* E_UL + a_DL .* E_DL;

    D_muref = ...
        d_C .* mu_c + d_L .* mu_L + d_R .* mu_R + d_U .* mu_U + d_D .* mu_D;

    R(row) = E_c/dt - AE_Eref - D_muref;

    % Move known background dmu through D*dmu_bg_l to the RHS.
    db_l = dmu_bg{l};
    D_bg = ...
        d_C .* db_l(idx_m) + ...
        d_L .* db_l(idx_L(idx_m)) + d_R .* db_l(idx_R(idx_m)) + ...
        d_U .* db_l(idx_U(idx_m)) + d_D .* db_l(idx_D(idx_m));

    R(row) = R(row) - D_bg;

    % Diffusion block on local dmu_l
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{l}(idx_m),d_C);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{l}(idx_L(idx_m)),d_L);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{l}(idx_R(idx_m)),d_R);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{l}(idx_U(idx_m)),d_U);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{l}(idx_D(idx_m)),d_D);

    % A_E * chi * dmu_loc, and subtract A_E * chi * dmu_bg
    for m = 1:Ne

        Chi_lm = chi_imp{l,m};

        Chi_c  = Chi_lm(idx_m);
        Chi_L  = Chi_lm(idx_L(idx_m));
        Chi_R  = Chi_lm(idx_R(idx_m));
        Chi_U  = Chi_lm(idx_U(idx_m));
        Chi_D  = Chi_lm(idx_D(idx_m));

        Chi_L2 = Chi_lm(idx_L2(idx_m));
        Chi_R2 = Chi_lm(idx_R2(idx_m));
        Chi_U2 = Chi_lm(idx_U2(idx_m));
        Chi_D2 = Chi_lm(idx_D2(idx_m));

        Chi_UR = Chi_lm(idx_UR(idx_m));
        Chi_DR = Chi_lm(idx_DR(idx_m));
        Chi_UL = Chi_lm(idx_UL(idx_m));
        Chi_DL = Chi_lm(idx_DL(idx_m));

        bg_m = dmu_bg{m};

        BG_AEChi = ...
            a_C  .* Chi_c  .* bg_m(idx_m) + ...
            a_L  .* Chi_L  .* bg_m(idx_L(idx_m)) + ...
            a_R  .* Chi_R  .* bg_m(idx_R(idx_m)) + ...
            a_U  .* Chi_U  .* bg_m(idx_U(idx_m)) + ...
            a_D  .* Chi_D  .* bg_m(idx_D(idx_m)) + ...
            a_L2 .* Chi_L2 .* bg_m(idx_L2(idx_m)) + ...
            a_R2 .* Chi_R2 .* bg_m(idx_R2(idx_m)) + ...
            a_U2 .* Chi_U2 .* bg_m(idx_U2(idx_m)) + ...
            a_D2 .* Chi_D2 .* bg_m(idx_D2(idx_m)) + ...
            a_UR .* Chi_UR .* bg_m(idx_UR(idx_m)) + ...
            a_DR .* Chi_DR .* bg_m(idx_DR(idx_m)) + ...
            a_UL .* Chi_UL .* bg_m(idx_UL(idx_m)) + ...
            a_DL .* Chi_DL .* bg_m(idx_DL(idx_m));

        R(row) = R(row) - BG_AEChi;

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_m),         a_C  .* Chi_c);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_L(idx_m)),  a_L  .* Chi_L);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_R(idx_m)),  a_R  .* Chi_R);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_U(idx_m)),  a_U  .* Chi_U);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_D(idx_m)),  a_D  .* Chi_D);

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_L2(idx_m)), a_L2 .* Chi_L2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_R2(idx_m)), a_R2 .* Chi_R2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_U2(idx_m)), a_U2 .* Chi_U2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_D2(idx_m)), a_D2 .* Chi_D2);

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_UR(idx_m)), a_UR .* Chi_UR);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_DR(idx_m)), a_DR .* Chi_DR);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_UL(idx_m)), a_UL .* Chi_UL);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMuMap{m}(idx_DL(idx_m)), a_DL .* Chi_DL);

    end

    % A_E * B * dphi
    for alpha = 1:Ngrain

        idmap = idPhiMap{alpha};

        if Nphi > 0 && ~any(idmap(:))
            continue
        end

        B_la = facPhi{alpha} .* (e0{alpha}{l} - e_bar{l});

        B_c  = B_la(idx_m);
        B_L  = B_la(idx_L(idx_m));
        B_R  = B_la(idx_R(idx_m));
        B_U  = B_la(idx_U(idx_m));
        B_D  = B_la(idx_D(idx_m));

        B_L2 = B_la(idx_L2(idx_m));
        B_R2 = B_la(idx_R2(idx_m));
        B_U2 = B_la(idx_U2(idx_m));
        B_D2 = B_la(idx_D2(idx_m));

        B_UR = B_la(idx_UR(idx_m));
        B_DR = B_la(idx_DR(idx_m));
        B_UL = B_la(idx_UL(idx_m));
        B_DL = B_la(idx_DL(idx_m));

        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_m,          a_C  .* B_c);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_L(idx_m),   a_L  .* B_L);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_R(idx_m),   a_R  .* B_R);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_U(idx_m),   a_U  .* B_U);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_D(idx_m),   a_D  .* B_D);

        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_L2(idx_m),  a_L2 .* B_L2);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_R2(idx_m),  a_R2 .* B_R2);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_U2(idx_m),  a_U2 .* B_U2);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_D2(idx_m),  a_D2 .* B_D2);

        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_UR(idx_m),  a_UR .* B_UR);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_DR(idx_m),  a_DR .* B_DR);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_UL(idx_m),  a_UL .* B_UL);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_DL(idx_m),  a_DL .* B_DL);

    end

end

% ------------------------------------------------------------
% Assemble and solve
% ------------------------------------------------------------
rows = rows(1:k-1);
cols = cols(1:k-1);
vals = vals(1:k-1);

A = sparse(rows,cols,vals,Ntot,Ntot);

use_order_cache = isfield(NUM,'use_order_cache') && NUM.use_order_cache == 1;

if use_order_cache
    [q,order_cache_hit,ORDER_CACHE] = CachedColamdOrdering(ORDER_CACHE,A,rows,cols,Ntot);
else
    q = colamd(A);
    order_cache_hit = 0;
end

y = A(:,q) \ R;

sol = zeros(size(R));
sol(q) = y;

flag   = 0;
relres = norm(A*sol - R)/max(norm(R),eps);
iter   = [0 0];

% ------------------------------------------------------------
% Unpack dphi and local dmu
% ------------------------------------------------------------
dphi = zeros(ny,nx,Ngrain);

for alpha = 1:Ngrain

    idmap  = idPhiMap{alpha};
    mask_a = maskPhi(:,:,alpha);

    ids = idmap(mask_a);

    tmp = zeros(ny,nx);

    if ~isempty(ids)
        tmp(mask_a) = sol(ids);
    end

    dphi(:,:,alpha) = tmp;

end

dmu_loc   = cell(1,Ne);
dmu_total = cell(1,Ne);

for ie = 1:Ne
    tmp = zeros(ny,nx);
    tmp(idsMuNode) = sol(idMuMap{ie}(idsMuNode));
    dmu_loc{ie}   = tmp;
    dmu_total{ie} = dmu_bg{ie} + dmu_loc{ie};
end

% ------------------------------------------------------------
% Update phi and p
% ------------------------------------------------------------
STATE_NEW     = STATE_REF;
STATE_NEW.phi = phi_ref + dphi;

if isfield(NUM,'norm_phi') && NUM.norm_phi == 1
    STATE_NEW.phi = Normalize_Phi_L2_Local(STATE_NEW.phi);
end

if isfield(NUM,'cut_phi') && NUM.cut_phi == 1
    STATE_NEW.phi = Cut_Phi_Local(STATE_NEW.phi);
end

STATE_NEW.p = Calc_p(MODEL,STATE_NEW.phi);

% ------------------------------------------------------------
% Update mu and recover E from total dmu and dp
% ------------------------------------------------------------
STATE_NEW.mu_e = mu_ref;

for ie = 1:Ne
    STATE_NEW.mu_e{ie} = mu_ref{ie} + dmu_total{ie};
end

STATE_NEW.E = E_ref;

for ie = 1:Ne

    En = E_ref{ie};

    for je = 1:Ne
        En = En + chi_imp{ie,je} .* dmu_total{je};
    end

    for ig = 1:Ngrain
        dp_ig = STATE_NEW.p(:,:,ig) - p_ref(:,:,ig);
        En = En + e0{ig}{ie} .* dp_ig;
    end

    STATE_NEW.E{ie} = En;

end

if isfield(NUM,'norm_E') && NUM.norm_E == 1
    STATE_NEW.E = EnforceMeanE_Local(STATE_NEW.E,E_ref);
end

% ------------------------------------------------------------
% Predict omega for diagnostics only.
% ------------------------------------------------------------
STATE_NEW.omg = STATE_REF.omg;

for ig = 1:Ngrain

    domega = zeros(ny,nx);

    for ie = 1:Ne
        domega = domega - e0{ig}{ie} .* dmu_total{ie};
    end

    STATE_NEW.omg(:,:,ig) = STATE_REF.omg(:,:,ig) + domega;

end

omg_mean = mean(STATE_NEW.omg,3);

for ig = 1:Ngrain
    STATE_NEW.omg(:,:,ig) = STATE_NEW.omg(:,:,ig) - omg_mean;
end

STATE_NEW.chi = chi0;
STATE_NEW.e   = e0;
if isfield(STATE_COEF,'c')
    STATE_NEW.c = STATE_COEF.c;
end
if isfield(STATE_COEF,'c_ext')
    STATE_NEW.c_ext = STATE_COEF.c_ext;
end

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
max_dmu_loc   = 0;
max_dmu_bg    = 0;
max_dmu_total = 0;
max_dE        = 0;

for ie = 1:Ne
    max_dmu_loc   = max(max_dmu_loc,max(abs(dmu_loc{ie}(:))));
    max_dmu_bg    = max(max_dmu_bg,max(abs(dmu_bg{ie}(:))));
    max_dmu_total = max(max_dmu_total,max(abs(dmu_total{ie}(:))));
    max_dE        = max(max_dE,max(abs(STATE_NEW.E{ie}(:)-STATE_REF.E{ie}(:))));
end

% Edge diagnostic for local correction only.
if band_width > 2
    mask_inner = Local_Dilate_Mask(mask_core,band_width-2);
    mask_edge  = maskMu & ~mask_inner;
else
    mask_edge  = maskMu & ~mask_core;
end

edge_dmu = 0;

if any(mask_edge(:))
    for ie = 1:Ne
        tmp = abs(dmu_loc{ie}(mask_edge));
        if ~isempty(tmp)
            edge_dmu = max(edge_dmu,max(tmp));
        end
    end
end

DIAG.flag              = flag;
DIAG.relres            = relres;
DIAG.iter              = iter;
DIAG.matrix_size       = Ntot;
DIAG.nnz               = nnz(A);
DIAG.use_order_cache   = use_order_cache;
DIAG.order_cache_hit   = order_cache_hit;
DIAG.Nphi_active       = Nphi;
DIAG.Nphi_full         = Ngrain*Nnode;
DIAG.Nmu_local         = Nmu;
DIAG.Nmu_full          = Ne*Nnode;
DIAG.local_mu_ratio    = Nmu/max(Ne*Nnode,1);
DIAG.Ngrain            = Ngrain;
DIAG.Ne                = Ne;
DIAG.active_phi_ratio  = Nphi / max(Ngrain*Nnode,1);
DIAG.active_cell_mean  = mean(active_per_cell(:));
DIAG.active_cell_max   = max(active_per_cell(:));
DIAG.band_width        = band_width;
DIAG.n_core_nodes      = nnz(mask_core);
DIAG.n_mu_nodes        = nnz(maskMu);
DIAG.mask_core         = mask_core;
DIAG.maskMu            = maskMu;
DIAG.max_dphi          = max(abs(dphi(:)));
DIAG.max_dmu_bg        = max_dmu_bg;
DIAG.max_dmu_loc       = max_dmu_loc;
DIAG.max_dmu           = max_dmu_total;
DIAG.max_dE            = max_dE;
DIAG.edge_dmu_ratio    = edge_dmu/max(max_dmu_loc,eps);
DIAG.kappa_on_frac     = mean(kappa_eff(:) > 0);
DIAG.kappa_max         = max(kappa_eff(:));
DIAG.kappa_mean        = mean(kappa_eff(:));
DIAG.use_Jphi          = ~isempty(Jphi_diag);
DIAG.use_CS_chi        = isfield(PARAM,'use_CS_chi') && PARAM.use_CS_chi == 1;
DIAG.use_kappa_c       = use_kappa_c;

end


% =========================================================================
% Local helper functions copied/simplified from the original solver
% =========================================================================

function mask_core = Local_Interface_Core_From_p(p,p_cut)

[ny,nx,Ng] = size(p);

mask_core = false(ny,nx);

for ig = 1:Ng

    p_ig  = p(:,:,ig);
    mixed = p_ig > p_cut & p_ig < 1-p_cut;

    if nx == 1
        p_L = p_ig;
        p_R = p_ig;
    else
        p_L = p_ig(:,[2,1:nx-1]);
        p_R = p_ig(:,[2:nx,nx-1]);
    end

    if ny == 1
        p_U = p_ig;
        p_D = p_ig;
    else
        p_U = p_ig([2,1:ny-1],:);
        p_D = p_ig([2:ny,ny-1],:);
    end

    jump = abs(p_ig-p_L) > p_cut | abs(p_ig-p_R) > p_cut | ...
           abs(p_ig-p_U) > p_cut | abs(p_ig-p_D) > p_cut;

    mask_core = mask_core | mixed | jump;

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


function [q,hit,CACHE] = CachedColamdOrdering(CACHE,A,rows,cols,Ntot)

hit = 0;

if ~isempty(CACHE) && isfield(CACHE,'Ntot') && CACHE.Ntot == Ntot && ...
        isfield(CACHE,'rows') && numel(CACHE.rows) == numel(rows) && ...
        isequal(CACHE.rows,rows) && isequal(CACHE.cols,cols)

    q   = CACHE.q;
    hit = 1;
    return

end

q = colamd(A);

CACHE        = struct();
CACHE.Ntot   = Ntot;
CACHE.rows   = rows;
CACHE.cols   = cols;
CACHE.q      = q;

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


function [q_C,q_L,q_R,q_U,q_D,q_L2,q_R2,q_U2,q_D2,q_UR,q_DR,q_UL,q_DL] = ...
    Kappa4_Coeffs(K_l,idx_c,idx_L,idx_R,idx_U,idx_D, ...
                  idx_L2,idx_R2,idx_U2,idx_D2,idx_UR,idx_DR,idx_UL,idx_DL, ...
                  dx2,dy2,dx4,dy4)

K_c = K_l(idx_c);
K_L = K_l(idx_L);
K_R = K_l(idx_R);
K_U = K_l(idx_U);
K_D = K_l(idx_D);

q_L  = -2*(K_L+K_c)/dx4 - 2*(K_L+K_c)/(dx2*dy2);
q_R  = -2*(K_R+K_c)/dx4 - 2*(K_R+K_c)/(dx2*dy2);
q_U  = -2*(K_U+K_c)/dy4 - 2*(K_U+K_c)/(dx2*dy2);
q_D  = -2*(K_D+K_c)/dy4 - 2*(K_D+K_c)/(dx2*dy2);

q_L2 = K_L/dx4;
q_R2 = K_R/dx4;
q_U2 = K_U/dy4;
q_D2 = K_D/dy4;

q_UR = (K_U+K_R)/(dx2*dy2);
q_DR = (K_D+K_R)/(dx2*dy2);
q_UL = (K_U+K_L)/(dx2*dy2);
q_DL = (K_D+K_L)/(dx2*dy2);

q_C  = (K_L+K_R)/dx4 + (K_U+K_D)/dy4 + ...
       4*K_c/dx4 + 4*K_c/dy4 + 8*K_c/(dx2*dy2);

end


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


function [rows,cols,vals,k] = add_block(rows,cols,vals,k,r,c,v)

row   = r(:);
col   = c(:);
coeff = v(:);

keep = col > 0 & isfinite(coeff) & coeff ~= 0;

if ~any(keep)
    return
end

row = row(keep);
col = col(keep);
val = coeff(keep);

n = numel(row);

if k+n-1 > numel(rows)
    grow = max(numel(rows),n + 1000);
    rows = [rows; zeros(grow,1)];
    cols = [cols; zeros(grow,1)];
    vals = [vals; zeros(grow,1)];
end

rows(k:k+n-1) = row;
cols(k:k+n-1) = col;
vals(k:k+n-1) = val;

k = k + n;

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


function E = EnforceMeanE_Local(E,E_old)

Ne = numel(E);

for ie = 1:Ne
    target_mean = mean(E_old{ie}(:));
    new_mean    = mean(E{ie}(:));
    E{ie}       = E{ie} + target_mean - new_mean;
end

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
