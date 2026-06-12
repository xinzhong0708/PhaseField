function [STATE_CORR,DIAG] = PF_CH_LECorrector_FixedP_Band_Mixed(STATE_TIME,STATE_IT,PARAM,GRID,PHYS,NUM)
%PF_CH_LECORRECTOR_FIXEDP_BAND_MIXED
%
% Mixed fixed-p CH corrector.
%
% Unknowns:
%   x = [dmu ; dE]
%
% Equations:
%
%   1. CH residual:
%        A_E*dE + D*dmu = E_TIME/dt - A_E*E_IT - D*mu_IT
%
%   2. Fixed-p tangent closure:
%        dE - chi_IT*dmu = 0
%
% This avoids inv(chi) and avoids recovering E after the solve.

E_time = STATE_TIME.E;
E_it   = STATE_IT.E;
mu_it  = STATE_IT.mu_e;
chi_it = STATE_IT.chi;

[ny,nx,~] = size(STATE_IT.p);
Ne        = numel(E_it);
Nnode     = ny*nx;

dt        = NUM.dt_phy;
dx        = GRID.dx;
dy        = GRID.dy;

dx2       = dx^2;
dy2       = dy^2;
dx4       = dx2^2;
dy4       = dy2^2;

if isfield(PARAM,'kappa_eff') && ~isempty(PARAM.kappa_eff)
    kappa_eff = PARAM.kappa_eff;
    if isscalar(kappa_eff)
        kappa_eff = kappa_eff*ones(ny,nx);
    end
else
    kappa_eff = PHYS.kappa*ones(ny,nx);
end

% ------------------------------------------------------------
% Reflective neighbour indices
% ------------------------------------------------------------
[Igrid,Jgrid] = ndgrid(1:ny,1:nx);

ii = Igrid(:);
jj = Jgrid(:);

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

idx_c  = sub2ind([ny,nx],ii,jj);
idx_L  = sub2ind([ny,nx],ii,jjL);
idx_R  = sub2ind([ny,nx],ii,jjR);
idx_U  = sub2ind([ny,nx],iiU,jj);
idx_D  = sub2ind([ny,nx],iiD,jj);

idx_L2 = sub2ind([ny,nx],ii,jjL2);
idx_R2 = sub2ind([ny,nx],ii,jjR2);
idx_U2 = sub2ind([ny,nx],iiU2,jj);
idx_D2 = sub2ind([ny,nx],iiD2,jj);

idx_UR = sub2ind([ny,nx],iiUR,jjUR);
idx_DR = sub2ind([ny,nx],iiDR,jjDR);
idx_UL = sub2ind([ny,nx],iiUL,jjUL);
idx_DL = sub2ind([ny,nx],iiDL,jjDL);

% ------------------------------------------------------------
% Unknown and row ids
% ------------------------------------------------------------
Nchem = Ne*Nnode;
Ntot  = 2*Nchem;

idMu = cell(1,Ne);
idE  = cell(1,Ne);

for ie = 1:Ne
    idMu{ie} = ((ie-1)*Nnode+(1:Nnode)).';
    idE{ie}  = Nchem + ((ie-1)*Nnode+(1:Nnode)).';
end

% Row ids:
%   CH rows use idMu
%   closure rows use idE

max_nnz = Nnode*Ne*(18+Ne) + 1000;

rows = zeros(max_nnz,1);
cols = zeros(max_nnz,1);
vals = zeros(max_nnz,1);
R    = zeros(Ntot,1);

k = 1;

% ============================================================
% 1. CH rows
% ============================================================
for l = 1:Ne

    row = idMu{l}(idx_c);

    Ml  = PARAM.M{l};
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

    kappa_c = kappa_eff(idx_c);
    MK_c    = M_c.*kappa_c;

    q_L  = MK_c.*(-4/dx4-4/(dx2*dy2));
    q_R  = MK_c.*(-4/dx4-4/(dx2*dy2));
    q_U  = MK_c.*(-4/dy4-4/(dx2*dy2));
    q_D  = MK_c.*(-4/dy4-4/(dx2*dy2));

    q_L2 = MK_c.*(1/dx4);
    q_R2 = MK_c.*(1/dx4);
    q_U2 = MK_c.*(1/dy4);
    q_D2 = MK_c.*(1/dy4);

    q_UR = MK_c.*(2/(dx2*dy2));
    q_DR = MK_c.*(2/(dx2*dy2));
    q_UL = MK_c.*(2/(dx2*dy2));
    q_DL = MK_c.*(2/(dx2*dy2));

    q_C  = MK_c.*(6/dx4+6/dy4+8/(dx2*dy2));

    a_C  = 1/dt+q_C;
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
        a_C.*E_c + ...
        a_L.*E_L + a_R.*E_R + a_U.*E_U + a_D.*E_D + ...
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

    % A_E*dE_l
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_c), a_C);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_L), a_L);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_R), a_R);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_U), a_U);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_D), a_D);

    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_L2),a_L2);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_R2),a_R2);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_U2),a_U2);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_D2),a_D2);

    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_UR),a_UR);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_DR),a_DR);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_UL),a_UL);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_DL),a_DL);
end

% ============================================================
% 2. Closure rows: dE - chi*dmu = 0
% ============================================================
for l = 1:Ne

    row = idE{l}(idx_c);

    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_c),ones(Nnode,1));

    for m = 1:Ne
        [rows,cols,vals,k] = add_block( ...
            rows,cols,vals,k,row,idMu{m}(idx_c),-chi_it{l,m}(idx_c));
    end
end

% ------------------------------------------------------------
% Assemble
% ------------------------------------------------------------
rows = rows(1:k-1);
cols = cols(1:k-1);
vals = vals(1:k-1);

A = sparse(rows,cols,vals,Ntot,Ntot);

% ------------------------------------------------------------
% Build correction band
% ------------------------------------------------------------
if isfield(NUM,'CHLE_p_cut') && ~isempty(NUM.CHLE_p_cut)
    p_cut = NUM.CHLE_p_cut;
else
    p_cut = 1e-8;
end

if isfield(NUM,'CHLE_band_thick') && ~isempty(NUM.CHLE_band_thick)
    band_thick = NUM.CHLE_band_thick;
else
    band_thick = 12;
end

if isfield(NUM,'CHLE_force_full') && NUM.CHLE_force_full == 1

    mask_core = true(ny,nx);
    mask_corr = true(ny,nx);

else

    mask_core = false(ny,nx);

    for ig = 1:size(STATE_IT.p,3)

        p_ig  = STATE_IT.p(:,:,ig);
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

    res_node = zeros(ny,nx);

    for ie = 1:Ne
        res_node = max(res_node,abs(reshape(R(idMu{ie}),ny,nx)));
    end

    if isfield(NUM,'CHLE_res_rel') && ~isempty(NUM.CHLE_res_rel) && NUM.CHLE_res_rel > 0
        res_max = max(res_node(:));

        if res_max > 0
            mask_core = mask_core | res_node > NUM.CHLE_res_rel*res_max;
        end
    end

    mask_corr = Dilate_Mask_Local(mask_core,band_thick);
end

% ------------------------------------------------------------
% Solve
% ------------------------------------------------------------
if ~any(mask_corr(:))

    sol = zeros(size(R));
    relres = 0;
    active = zeros(0,1);
    Ared = sparse(0,0);
    Rred = zeros(0,1);
    outside_res_ratio = norm(R)/max(norm(R),eps);

else

    active_node = mask_corr(:);
    active      = zeros(2*Ne*nnz(active_node),1);
    kact        = 1;

    for ie = 1:Ne
        ids = idMu{ie}(active_node);
        active(kact:kact+numel(ids)-1) = ids;
        kact = kact+numel(ids);
    end

    for ie = 1:Ne
        ids = idE{ie}(active_node);
        active(kact:kact+numel(ids)-1) = ids;
        kact = kact+numel(ids);
    end

    Ared = A(active,active);
    Rred = R(active);

    q = colamd(Ared);
    y = Ared(:,q)\Rred;

    sol_red = zeros(size(Rred));
    sol_red(q) = y;

    sol = zeros(size(R));
    sol(active) = sol_red;

    relres = norm(Ared*sol_red-Rred)/max(norm(Rred),eps);

    outside = true(Ntot,1);
    outside(active) = false;
    outside_res_ratio = norm(R(outside))/max(norm(R),eps);
end

% ------------------------------------------------------------
% Unpack dmu and dE
% ------------------------------------------------------------
dmu = cell(1,Ne);
dE  = cell(1,Ne);

for ie = 1:Ne
    dmu{ie} = reshape(sol(idMu{ie}),ny,nx);
    dE{ie}  = reshape(sol(idE{ie}),ny,nx);
end

% ------------------------------------------------------------
% Corrected state
% ------------------------------------------------------------
STATE_CORR     = STATE_IT;
STATE_CORR.phi = STATE_IT.phi;
STATE_CORR.p   = STATE_IT.p;

for ie = 1:Ne
    STATE_CORR.mu_e{ie} = STATE_IT.mu_e{ie}+dmu{ie};
    STATE_CORR.E{ie}    = STATE_IT.E{ie}+dE{ie};
end

if isfield(NUM,'norm_E') && NUM.norm_E == 1
    STATE_CORR.E = EnforceMeanE_Local(STATE_CORR.E,STATE_TIME.E);
end

STATE_CORR.omg = STATE_IT.omg;

for ig = 1:size(STATE_IT.omg,3)

    domega = zeros(ny,nx);

    for ie = 1:Ne
        domega = domega-STATE_IT.e{ig}{ie}.*dmu{ie};
    end

    STATE_CORR.omg(:,:,ig) = STATE_IT.omg(:,:,ig)+domega;
end

omg_mean = mean(STATE_CORR.omg,3);

for ig = 1:size(STATE_CORR.omg,3)
    STATE_CORR.omg(:,:,ig) = STATE_CORR.omg(:,:,ig)-omg_mean;
end

STATE_CORR.chi = STATE_IT.chi;
STATE_CORR.e   = STATE_IT.e;
STATE_CORR.c   = STATE_IT.c;

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
DIAG.relres             = relres;
DIAG.matrix_size_full   = Ntot;
DIAG.nnz_full           = nnz(A);
DIAG.matrix_size        = numel(active);
DIAG.nnz                = nnz(Ared);
DIAG.solve_ratio        = numel(active)/max(Ntot,1);
DIAG.outside_res_ratio  = outside_res_ratio;
DIAG.mask_core          = mask_core;
DIAG.mask_corr          = mask_corr;
DIAG.Ne                 = Ne;
DIAG.band_thick         = band_thick;
DIAG.full_coverage      = all(mask_corr(:));
DIAG.n_core_nodes       = nnz(mask_core);
DIAG.n_corr_nodes       = nnz(mask_corr);

max_dmu = 0;
max_dE  = 0;

for ie = 1:Ne
    max_dmu = max(max_dmu,max(abs(dmu{ie}(:))));
    max_dE  = max(max_dE, max(abs(dE{ie}(:))));
end

DIAG.max_dmu = max_dmu;
DIAG.max_dE  = max_dE;
DIAG.fixed_phi_error = max(abs(STATE_CORR.phi(:)-STATE_IT.phi(:)));
DIAG.fixed_p_error   = max(abs(STATE_CORR.p(:)-STATE_IT.p(:)));

end



%% ========================================================================
% Local helper functions
% ========================================================================


function mask = Dilate_Mask_Local(core,thickness)
%DILATE_MASK_LOCAL Dilate a node mask by a square stencil.

if thickness <= 0
    mask = logical(core);
    return
end

ker  = ones(2*thickness+1,2*thickness+1);
mask = conv2(double(core),ker,'same') > 0;

end

function idx = reflect_index(idx,n)

if n == 1
    idx = ones(size(idx));
    return
end

period = 2*n-2;
r      = mod(idx-1,period);
idx    = 1+min(r,period-r);

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
k             = k+n;

end


function E = EnforceMeanE_Local(E,E_old)

Ne = numel(E);

for ie = 1:Ne

    target_mean = mean(E_old{ie}(:));
    new_mean    = mean(E{ie}(:));

    E{ie}       = E{ie}+target_mean-new_mean;

end

end
