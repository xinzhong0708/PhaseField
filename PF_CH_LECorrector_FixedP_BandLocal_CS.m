function [STATE_CORR,DIAG] = PF_CH_LECorrector_FixedP_BandLocal_CS(STATE_TIME,STATE_IT,PARAM,MODEL,GRID,PHYS,NUM)
%PF_CH_LECORRECTOR_FIXEDP_BANDLOCAL_CS Fixed-p CH corrector with local assembly.
%
% This is the fast version of PF_CH_LECorrector_FixedP_Band_CS.
% It builds the correction band first, then assembles only the active rows
% and active columns. Outside the band, dmu_corr = 0.
%
% Main limitation:
%   Fast local assembly currently supports the standard E-kappa CH operator.
%   If PARAM.use_kappa_c = 1, this function falls back to the original
%   PF_CH_LECorrector_FixedP_Band_CS unless NUM.CH_fast_disable_kappac = 1.
%
% Recommended first use inside PF_CH_CoarseFine_Corrector_CS_fast:
%   NUM.CH_fast_disable_kappac = 1;   % local fine repair uses standard CH
%
% The coarse global solve can still be done by the original corrector.

if nargin < 7
    error('PF_CH_LECorrector_FixedP_BandLocal_CS requires STATE_TIME, STATE_IT, PARAM, MODEL, GRID, PHYS and NUM.')
end

% Backward compatibility with the temporary signature:
%   (STATE_TIME,STATE_IT,MODEL,PARAM,GRID,PHYS,NUM)
if isfield(PARAM,'pars') && isfield(MODEL,'M')
    tmp   = PARAM;
    PARAM = MODEL;
    MODEL = tmp;
end

use_kappa_c = isfield(PARAM,'use_kappa_c') && PARAM.use_kappa_c == 1;

if use_kappa_c
    if isfield(NUM,'CH_fast_disable_kappac') && NUM.CH_fast_disable_kappac == 1
        PARAM.use_kappa_c = 0;
        use_kappa_c = 0;
    else
        % Safe fallback.  This avoids silently changing the equation.
        [STATE_CORR,DIAG] = PF_CH_LECorrector_FixedP_Band_CS(STATE_TIME,STATE_IT,PARAM,MODEL,GRID,PHYS,NUM);
        DIAG.fast_local_fallback = 1;
        return
    end
end

E_time = STATE_TIME.E;
E_it   = STATE_IT.E;
mu_it  = STATE_IT.mu_e;
chi_it = STATE_IT.chi;

if numel(E_time) ~= numel(E_it) || numel(E_it) ~= numel(mu_it)
    error('STATE_TIME and STATE_IT have inconsistent numbers of elemental fields.')
end

[ny,nx,~] = size(STATE_IT.p);
Ne        = numel(E_it);
Nnode     = ny*nx;

dt        = NUM.dt_phy;
dx        = GRID.dx;
dy        = GRID.dy;

if dt <= 0 || dx <= 0 || dy <= 0
    error('NUM.dt_phy, GRID.dx and GRID.dy must be positive.')
end

dx2       = dx^2;
dy2       = dy^2;
dx4       = dx2^2;
dy4       = dy2^2;

% ------------------------------------------------------------
% Spatially varying kappa
% ------------------------------------------------------------
if isfield(PARAM,'kappa_eff') && ~isempty(PARAM.kappa_eff)
    kappa_eff = PARAM.kappa_eff;
    if isscalar(kappa_eff)
        kappa_eff = kappa_eff*ones(ny,nx);
    elseif ndims(kappa_eff) == 3
        kappa_eff = max(kappa_eff,[],3);
    end
else
    if isfield(PHYS,'kap') && ~isempty(PHYS.kap)
        kappa0 = PHYS.kap;
    elseif isfield(PHYS,'kappa') && ~isempty(PHYS.kappa)
        kappa0 = PHYS.kappa;
    else
        kappa0 = 0;
    end
    kappa_eff = kappa0*ones(ny,nx);
end

% ------------------------------------------------------------
% Convex-split implicit capacity
% ------------------------------------------------------------
if isfield(PARAM,'use_CS_chi') && PARAM.use_CS_chi == 1
    chi_imp = ConvexifyChi_ForImplicit(chi_it,PARAM);
else
    chi_imp = chi_it;
end

% ------------------------------------------------------------
% Build correction band before assembly
% ------------------------------------------------------------
if isfield(NUM,'CHLE_p_cut') && ~isempty(NUM.CHLE_p_cut)
    p_cut = NUM.CHLE_p_cut;
else
    p_cut = 1e-8;
end

if isfield(NUM,'CH_band_width') && ~isempty(NUM.CH_band_width)
    band_thick = NUM.CH_band_width;
elseif isfield(NUM,'CHLE_band_thick') && ~isempty(NUM.CHLE_band_thick)
    band_thick = NUM.CHLE_band_thick;
else
    band_thick = 12;
end

if band_thick < 2
    error('NUM.CHLE_band_thick must be at least 2 because the CH operator contains a +/-2 stencil.')
end

if isfield(NUM,'CHLE_force_full') && NUM.CHLE_force_full == 1
    mask_core = true(ny,nx);
    mask_corr = true(ny,nx);
else
    mask_core = Build_Interface_Core(STATE_IT.p,p_cut);

    % Include kappa/spinodal region in the local fine repair by default.
    include_kappa = 1;
    if isfield(NUM,'CHLE_include_kappa_mask')
        include_kappa = NUM.CHLE_include_kappa_mask;
    end
    if include_kappa == 1 && any(kappa_eff(:) > 0)
        mask_core = mask_core | (kappa_eff > 0);
    end

    mask_corr = Dilate_Mask_Local(mask_core,band_thick);
end

if ~any(mask_corr(:))
    STATE_CORR = STATE_IT;
    DIAG = Empty_Diag(Ne,Nnode,mask_core,mask_corr,band_thick);
    return
end

active_node = mask_corr(:);
active_ids  = find(active_node);
Na          = numel(active_ids);

node_id = zeros(Nnode,1);
node_id(active_ids) = (1:Na).';

idMu = cell(1,Ne);
for ie = 1:Ne
    tmp = zeros(Nnode,1);
    tmp(active_ids) = (ie-1)*Na + (1:Na).';
    idMu{ie} = tmp;
end

Nmu = Ne*Na;

% ------------------------------------------------------------
% Reflective-neighbour indices for active row nodes only
% ------------------------------------------------------------
[ii,jj] = ind2sub([ny,nx],active_ids);
ii = ii(:);
jj = jj(:);

refI = @(i,sh) reflect_index(i+sh,ny);
refJ = @(j,sh) reflect_index(j+sh,nx);

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

idx_c  = active_ids;
idx_L  = sub2ind([ny,nx],ii,  jjL);
idx_R  = sub2ind([ny,nx],ii,  jjR);
idx_U  = sub2ind([ny,nx],iiU, jj);
idx_D  = sub2ind([ny,nx],iiD, jj);

idx_L2 = sub2ind([ny,nx],ii,   jjL2);
idx_R2 = sub2ind([ny,nx],ii,   jjR2);
idx_U2 = sub2ind([ny,nx],iiU2, jj);
idx_D2 = sub2ind([ny,nx],iiD2, jj);

idx_UR = sub2ind([ny,nx],iiUR,jjUR);
idx_DR = sub2ind([ny,nx],iiDR,jjDR);
idx_UL = sub2ind([ny,nx],iiUL,jjUL);
idx_DL = sub2ind([ny,nx],iiDL,jjDL);

% ------------------------------------------------------------
% Sparse allocation and residual, active rows only
% ------------------------------------------------------------
max_nnz = Na*Ne*(5 + 13*Ne) + 1000;
rows    = zeros(max_nnz,1);
cols    = zeros(max_nnz,1);
vals    = zeros(max_nnz,1);
R       = zeros(Nmu,1);
k       = 1;

% ============================================================
% CH-only Newton correction block at fixed updated p
% ============================================================
for l = 1:Ne

    row = ((l-1)*Na + (1:Na)).';

    Ml = PARAM.M{l};
    if isscalar(Ml)
        Ml = Ml*ones(ny,nx);
    end

    M_c = Ml(idx_c);
    M_L = Ml(idx_L);
    M_R = Ml(idx_R);
    M_U = Ml(idx_U);
    M_D = Ml(idx_D);

    d_L = -(M_L+M_c)/2/dx2;
    d_R = -(M_R+M_c)/2/dx2;
    d_U = -(M_U+M_c)/2/dy2;
    d_D = -(M_D+M_c)/2/dy2;
    d_C = -(d_L+d_R+d_U+d_D);

    [q_C,q_L,q_R,q_U,q_D,q_L2,q_R2,q_U2,q_D2,q_UR,q_DR,q_UL,q_DL] = ...
        Kappa4_Coeffs(Ml.*kappa_eff,idx_c,idx_L,idx_R,idx_U,idx_D, ...
                      idx_L2,idx_R2,idx_U2,idx_D2,idx_UR,idx_DR,idx_UL,idx_DL, ...
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

    E_time_l = E_time{l};
    E_l      = E_it{l};
    mu_l     = mu_it{l};

    Etime_c = E_time_l(idx_c);

    E_c  = E_l(idx_c);
    E_L  = E_l(idx_L);
    E_R  = E_l(idx_R);
    E_U  = E_l(idx_U);
    E_D  = E_l(idx_D);
    E_L2 = E_l(idx_L2);
    E_R2 = E_l(idx_R2);
    E_U2 = E_l(idx_U2);
    E_D2 = E_l(idx_D2);
    E_UR = E_l(idx_UR);
    E_DR = E_l(idx_DR);
    E_UL = E_l(idx_UL);
    E_DL = E_l(idx_DL);

    mu_c = mu_l(idx_c);
    mu_L = mu_l(idx_L);
    mu_R = mu_l(idx_R);
    mu_U = mu_l(idx_U);
    mu_D = mu_l(idx_D);

    AE_Eit = ...
        a_C .*E_c  + ...
        a_L .*E_L  + a_R .*E_R  + a_U .*E_U  + a_D .*E_D  + ...
        a_L2.*E_L2 + a_R2.*E_R2 + a_U2.*E_U2 + a_D2.*E_D2 + ...
        a_UR.*E_UR + a_DR.*E_DR + a_UL.*E_UL + a_DL.*E_DL;

    D_muit = ...
        d_C.*mu_c + d_L.*mu_L + d_R.*mu_R + d_U.*mu_U + d_D.*mu_D;

    R(row) = Etime_c/dt - AE_Eit - D_muit;

    % D*dmu_l
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_c),d_C);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_L),d_L);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_R),d_R);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_U),d_U);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_D),d_D);

    % A_E*chi*dmu
    for m = 1:Ne

        Chi_lm = chi_imp{l,m};
        if isscalar(Chi_lm)
            Chi_lm = Chi_lm*ones(ny,nx);
        end

        Chi_c  = Chi_lm(idx_c);
        Chi_L  = Chi_lm(idx_L);
        Chi_R  = Chi_lm(idx_R);
        Chi_U  = Chi_lm(idx_U);
        Chi_D  = Chi_lm(idx_D);
        Chi_L2 = Chi_lm(idx_L2);
        Chi_R2 = Chi_lm(idx_R2);
        Chi_U2 = Chi_lm(idx_U2);
        Chi_D2 = Chi_lm(idx_D2);
        Chi_UR = Chi_lm(idx_UR);
        Chi_DR = Chi_lm(idx_DR);
        Chi_UL = Chi_lm(idx_UL);
        Chi_DL = Chi_lm(idx_DL);

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_c), a_C .*Chi_c);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_L), a_L .*Chi_L);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_R), a_R .*Chi_R);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_U), a_U .*Chi_U);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_D), a_D .*Chi_D);

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_L2),a_L2.*Chi_L2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_R2),a_R2.*Chi_R2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_U2),a_U2.*Chi_U2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_D2),a_D2.*Chi_D2);

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_UR),a_UR.*Chi_UR);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_DR),a_DR.*Chi_DR);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_UL),a_UL.*Chi_UL);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_DL),a_DL.*Chi_DL);

    end
end

% ------------------------------------------------------------
% Assemble and solve reduced system
% ------------------------------------------------------------
rows = rows(1:k-1);
cols = cols(1:k-1);
vals = vals(1:k-1);

A = sparse(rows,cols,vals,Nmu,Nmu);

if isempty(R) || Nmu == 0
    sol = zeros(Nmu,1);
    relres = 0;
else
    q = colamd(A);
    y = A(:,q)\R;
    sol = zeros(size(R));
    sol(q) = y;
    relres = norm(A*sol - R)/max(norm(R),eps);
end

% ------------------------------------------------------------
% Unpack chemical-potential correction to full fine grid
% ------------------------------------------------------------
dmu = cell(1,Ne);

for ie = 1:Ne
    tmp = zeros(Nnode,1);
    tmp(active_ids) = sol((ie-1)*Na + (1:Na));
    dmu{ie} = reshape(tmp,ny,nx);
end

% ------------------------------------------------------------
% Fixed-p corrected state
% ------------------------------------------------------------
STATE_CORR      = STATE_IT;
STATE_CORR.phi  = STATE_IT.phi;
STATE_CORR.p    = STATE_IT.p;

for ie = 1:Ne
    STATE_CORR.mu_e{ie} = STATE_IT.mu_e{ie} + dmu{ie};
end

for ie = 1:Ne
    En = STATE_IT.E{ie};
    for je = 1:Ne
        En = En + chi_imp{ie,je}.*dmu{je};
    end
    STATE_CORR.E{ie} = En;
end

if isfield(NUM,'norm_E') && NUM.norm_E == 1
    STATE_CORR.E = EnforceMeanE_Local(STATE_CORR.E,STATE_TIME.E);
end

STATE_CORR.omg = STATE_IT.omg;

for ig = 1:size(STATE_IT.omg,3)
    domega = zeros(ny,nx);
    for ie = 1:Ne
        domega = domega - STATE_IT.e{ig}{ie}.*dmu{ie};
    end
    STATE_CORR.omg(:,:,ig) = STATE_IT.omg(:,:,ig) + domega;
end

omg_mean = mean(STATE_CORR.omg,3);
for ig = 1:size(STATE_CORR.omg,3)
    STATE_CORR.omg(:,:,ig) = STATE_CORR.omg(:,:,ig) - omg_mean;
end

STATE_CORR.chi = STATE_IT.chi;
STATE_CORR.e   = STATE_IT.e;
STATE_CORR.c   = STATE_IT.c;
if isfield(STATE_IT,'c_ext')
    STATE_CORR.c_ext = STATE_IT.c_ext;
end

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
max_dmu = 0;
max_dE  = 0;
for ie = 1:Ne
    max_dmu = max(max_dmu,max(abs(dmu{ie}(:))));
    max_dE  = max(max_dE,max(abs(STATE_CORR.E{ie}(:)-STATE_IT.E{ie}(:))));
end

if band_thick > 2
    mask_inner = Dilate_Mask_Local(mask_core,band_thick-2);
    mask_edge  = mask_corr & ~mask_inner;
else
    mask_edge  = mask_corr & ~mask_core;
end

edge_dmu = 0;
if any(mask_edge(:))
    for ie = 1:Ne
        tmp = abs(dmu{ie}(mask_edge));
        if ~isempty(tmp)
            edge_dmu = max(edge_dmu,max(tmp));
        end
    end
end

DIAG.relres            = relres;
DIAG.matrix_size_full  = Ne*Nnode;
DIAG.matrix_size       = Nmu;
DIAG.nnz               = nnz(A);
DIAG.solve_ratio       = Nmu/max(Ne*Nnode,1);
DIAG.mask_core         = mask_core;
DIAG.mask_corr         = mask_corr;
DIAG.Ne                = Ne;
DIAG.band_thick        = band_thick;
DIAG.full_coverage     = all(mask_corr(:));
DIAG.n_core_nodes      = nnz(mask_core);
DIAG.n_corr_nodes      = nnz(mask_corr);
DIAG.max_dmu           = max_dmu;
DIAG.max_dE            = max_dE;
DIAG.edge_dmu_ratio    = edge_dmu/max(max_dmu,eps);
DIAG.fixed_phi_error   = max(abs(STATE_CORR.phi(:)-STATE_IT.phi(:)));
DIAG.fixed_p_error     = max(abs(STATE_CORR.p(:)-STATE_IT.p(:)));
DIAG.use_CS_chi        = isfield(PARAM,'use_CS_chi') && PARAM.use_CS_chi == 1;
DIAG.use_kappa_c       = 0;
DIAG.fast_local        = 1;
DIAG.fast_local_fallback = 0;

end


function mask_core = Build_Interface_Core(p,p_cut)

[ny,nx,Ng] = size(p);
mask_core = false(ny,nx);

for ig = 1:Ng

    p_ig = p(:,:,ig);
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


function DIAG = Empty_Diag(Ne,Nnode,mask_core,mask_corr,band_thick)

DIAG.relres            = 0;
DIAG.matrix_size_full  = Ne*Nnode;
DIAG.matrix_size       = 0;
DIAG.nnz               = 0;
DIAG.solve_ratio       = 0;
DIAG.mask_core         = mask_core;
DIAG.mask_corr         = mask_corr;
DIAG.Ne                = Ne;
DIAG.band_thick        = band_thick;
DIAG.full_coverage     = false;
DIAG.n_core_nodes      = nnz(mask_core);
DIAG.n_corr_nodes      = nnz(mask_corr);
DIAG.max_dmu           = 0;
DIAG.max_dE            = 0;
DIAG.edge_dmu_ratio    = 0;
DIAG.fixed_phi_error   = 0;
DIAG.fixed_p_error     = 0;
DIAG.fast_local        = 1;
DIAG.fast_local_fallback = 0;

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


function mask = Dilate_Mask_Local(core,thickness)

if thickness <= 0
    mask = logical(core);
    return
end

ker  = ones(2*thickness+1,2*thickness+1);
mask = conv2(double(core),ker,'same') > 0;

end


function [rows,cols,vals,k] = add_block(rows,cols,vals,k,r,c,v)

row   = r(:);
col   = c(:);
coeff = v(:);

keep  = col > 0 & isfinite(coeff) & coeff ~= 0;

if ~any(keep)
    return
end

row = row(keep);
col = col(keep);
val = coeff(keep);
n   = numel(row);

if k+n-1 > numel(rows)
    grow = max(numel(rows),n+1000);
    rows = [rows; zeros(grow,1)];
    cols = [cols; zeros(grow,1)];
    vals = [vals; zeros(grow,1)];
end

rows(k:k+n-1) = row;
cols(k:k+n-1) = col;
vals(k:k+n-1) = val;
k = k+n;

end


function E = EnforceMeanE_Local(E,E_old)

Ne = numel(E);

for ie = 1:Ne
    target_mean = mean(E_old{ie}(:));
    new_mean    = mean(E{ie}(:));
    E{ie}       = E{ie}+target_mean-new_mean;
end

end
