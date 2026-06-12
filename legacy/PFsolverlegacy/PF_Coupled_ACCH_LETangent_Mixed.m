function [STATE_NEW,DIAG] = PF_Coupled_ACCH_LETangent_Mixed(STATE_REF,PARAM,MODEL,GRID,PHYS,NUM,STATE_COEF)
%PF_COUPLED_ACCH_LETANGENT_MIXED
%
% Mixed AC-CH tangent solve.
%
% Unknowns:
%   x = [active dphi_grain ; global dmu ; global dE]
%
% Equations:
%
%   1. Allen-Cahn:
%        AC(dphi) - L*B*dmu = RHS_AC
%
%   2. Cahn-Hilliard:
%        A_E*dE + D*dmu = E_ref/dt - A_E*E_ref - D*mu_ref
%
%   3. Local tangent closure:
%        dE - chi*dmu - B*dphi = 0
%
% This avoids post-solve E reconstruction and avoids inv(chi).
% Final LE_Run_Mode_New should overwrite c, e, chi, omega.

if nargin < 7 || isempty(STATE_COEF)
    STATE_COEF = STATE_REF;
end

% ------------------------------------------------------------
% Fixed reference state
% ------------------------------------------------------------
phi_ref = STATE_REF.phi;
p_ref   = STATE_REF.p;
E_ref   = STATE_REF.E;
mu_ref  = STATE_REF.mu_e;

% ------------------------------------------------------------
% Coefficients
% ------------------------------------------------------------
chi0 = STATE_COEF.chi;
e0   = STATE_COEF.e;

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
% Spatially varying kappa
% ------------------------------------------------------------
if isfield(PARAM,'kappa_eff') && ~isempty(PARAM.kappa_eff)
    kappa_eff = PARAM.kappa_eff;
    if isscalar(kappa_eff)
        kappa_eff = kappa_eff*ones(ny,nx);
    end
else
    kappa_eff = PHYS.kappa*ones(ny,nx);
end

% ------------------------------------------------------------
% Build AC source using coefficient thermodynamics
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
% Active phi mask
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
    if isfield(PHYS,'kappa') && PHYS.kappa ~= 0
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
    counter = counter+nids;
end

Nphi  = counter;
Nchem = Ne*Nnode;
Ntot  = Nphi + 2*Nchem;

idMu = cell(1,Ne);
idE  = cell(1,Ne);

for ie = 1:Ne
    idMu{ie} = Nphi + ((ie-1)*Nnode+(1:Nnode)).';
    idE{ie}  = Nphi + Nchem + ((ie-1)*Nnode+(1:Nnode)).';
end

% Row ids:
%   AC rows use idPhiMap
%   CH rows use idMu
%   closure rows use idE

% ------------------------------------------------------------
% Tangent ingredients B_ie_alpha
% ------------------------------------------------------------
eps_phi = 1e-14;
Dphi    = sum(phi_ref.^2,3) + eps_phi;
p_tan   = phi_ref.^2 ./ Dphi;

facPhi = cell(1,Ngrain);

for alpha = 1:Ngrain
    facPhi{alpha} = 2.*phi_ref(:,:,alpha)./Dphi;
end

e_bar = cell(1,Ne);

for ie = 1:Ne

    tmp = zeros(ny,nx);

    for ig = 1:Ngrain
        tmp = tmp + p_tan(:,:,ig).*e0{ig}{ie};
    end

    e_bar{ie} = tmp;
end

B = cell(Ne,Ngrain);

for ie = 1:Ne
    for alpha = 1:Ngrain
        B{ie,alpha} = facPhi{alpha}.*(e0{alpha}{ie}-e_bar{ie});
    end
end

% ------------------------------------------------------------
% Sparse allocation
% ------------------------------------------------------------
max_nnz = ...
    Nphi*(6+Ne) + ...
    Nnode*Ne*(18+Ne+Ngrain) + ...
    1000;

rows = zeros(max_nnz,1);
cols = zeros(max_nnz,1);
vals = zeros(max_nnz,1);
R    = zeros(Ntot,1);

k = 1;

% ============================================================
% 1. Allen-Cahn rows
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

    rhs_full = PARAM.LK.*lap_phi_ref + S_AC{alpha};
    R(row)   = rhs_full(idx_a);

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

    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idPhiMap{alpha},idx_a, cC);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idPhiMap{alpha},idxa_L,cL);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idPhiMap{alpha},idxa_R,cR);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idPhiMap{alpha},idxa_U,cU);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idPhiMap{alpha},idxa_D,cD);

    if ~isempty(Jphi_diag)
        Jaa = Jphi_diag{alpha};
        [rows,cols,vals,k] = add_active_block( ...
            rows,cols,vals,k,row,idPhiMap{alpha},idx_a,-Jaa(idx_a));
    end

    % Coupling to dmu:
    %   -L * B_ie_alpha * dmu_ie
    for ie = 1:Ne

        coeff_mu = -PARAM.L.*B{ie,alpha};

        [rows,cols,vals,k] = add_block( ...
            rows,cols,vals,k,row,idMu{ie}(idx_a),coeff_mu(idx_a));
    end
end

% ============================================================
% 2. CH rows: A_E*dE + D*dmu = residual
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

    E_l  = E_ref{l};
    mu_l = mu_ref{l};

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

    AE_Eref = ...
        a_C.*E_c + ...
        a_L.*E_L + a_R.*E_R + a_U.*E_U + a_D.*E_D + ...
        a_L2.*E_L2 + a_R2.*E_R2 + a_U2.*E_U2 + a_D2.*E_D2 + ...
        a_UR.*E_UR + a_DR.*E_DR + a_UL.*E_UL + a_DL.*E_DL;

    D_muref = ...
        d_C.*mu_c + d_L.*mu_L + d_R.*mu_R + d_U.*mu_U + d_D.*mu_D;

    R(row) = E_c/dt - AE_Eref - D_muref;

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
% 3. Closure rows: dE - chi*dmu - B*dphi = 0
% ============================================================
for l = 1:Ne

    row = idE{l}(idx_c);

    % dE_l
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_c),ones(Nnode,1));

    % -chi_lm*dmu_m
    for m = 1:Ne
        [rows,cols,vals,k] = add_block( ...
            rows,cols,vals,k,row,idMu{m}(idx_c),-chi0{l,m}(idx_c));
    end

    % -B_lalpha*dphi_alpha
    for alpha = 1:Ngrain

        idmap = idPhiMap{alpha};

        if ~any(idmap(:))
            continue
        end

        [rows,cols,vals,k] = add_active_block( ...
            rows,cols,vals,k,row,idmap,idx_c,-B{l,alpha}(idx_c));
    end
end

% ------------------------------------------------------------
% Assemble and solve
% ------------------------------------------------------------
rows = rows(1:k-1);
cols = cols(1:k-1);
vals = vals(1:k-1);

A = sparse(rows,cols,vals,Ntot,Ntot);

q = colamd(A);
y = A(:,q)\R;

sol = zeros(size(R));
sol(q) = y;

flag   = 0;
relres = norm(A*sol-R)/max(norm(R),eps);
iter   = [0 0];

% ------------------------------------------------------------
% Unpack dphi, dmu, dE
% ------------------------------------------------------------
dphi = zeros(ny,nx,Ngrain);

for alpha = 1:Ngrain

    idmap  = idPhiMap{alpha};
    mask_a = maskPhi(:,:,alpha);
    ids    = idmap(mask_a);

    tmp = zeros(ny,nx);

    if ~isempty(ids)
        tmp(mask_a) = sol(ids);
    end

    dphi(:,:,alpha) = tmp;
end

dmu = cell(1,Ne);
dE  = cell(1,Ne);

for ie = 1:Ne
    dmu{ie} = reshape(sol(idMu{ie}),ny,nx);
    dE{ie}  = reshape(sol(idE{ie}),ny,nx);
end

% ------------------------------------------------------------
% Update state
% ------------------------------------------------------------
STATE_NEW     = STATE_REF;
STATE_NEW.phi = phi_ref+dphi;

if isfield(NUM,'norm_phi') && NUM.norm_phi == 1
    STATE_NEW.phi = Normalize_Phi_L2_Local(STATE_NEW.phi);
end

if isfield(NUM,'cut_phi') && NUM.cut_phi == 1
    STATE_NEW.phi = Cut_Phi_Local(STATE_NEW.phi);
end

STATE_NEW.p = Calc_p(MODEL,STATE_NEW.phi);

STATE_NEW.mu_e = mu_ref;
STATE_NEW.E    = E_ref;

for ie = 1:Ne
    STATE_NEW.mu_e{ie} = mu_ref{ie}+dmu{ie};
    STATE_NEW.E{ie}    = E_ref{ie}+dE{ie};
end

if isfield(NUM,'norm_E') && NUM.norm_E == 1
    STATE_NEW.E = EnforceMeanE_Local(STATE_NEW.E,E_ref);
end

% Diagnostics-only omega prediction
STATE_NEW.omg = STATE_REF.omg;

for ig = 1:Ngrain

    domega = zeros(ny,nx);

    for ie = 1:Ne
        domega = domega-e0{ig}{ie}.*dmu{ie};
    end

    STATE_NEW.omg(:,:,ig) = STATE_REF.omg(:,:,ig)+domega;
end

omg_mean = mean(STATE_NEW.omg,3);

for ig = 1:Ngrain
    STATE_NEW.omg(:,:,ig) = STATE_NEW.omg(:,:,ig)-omg_mean;
end

STATE_NEW.chi = chi0;
STATE_NEW.e   = e0;

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
DIAG.flag     = flag;
DIAG.relres   = relres;
DIAG.iter     = iter;
DIAG.max_dphi = max(abs(dphi(:)));

max_dmu = 0;
max_dE  = 0;

for ie = 1:Ne
    max_dmu = max(max_dmu,max(abs(dmu{ie}(:))));
    max_dE  = max(max_dE, max(abs(dE{ie}(:))));
end

DIAG.max_dmu          = max_dmu;
DIAG.max_dE           = max_dE;
DIAG.matrix_size      = Ntot;
DIAG.nnz              = nnz(A);
DIAG.Nphi_active      = Nphi;
DIAG.Nphi_full        = Ngrain*Nnode;
DIAG.Nmu              = Nchem;
DIAG.NE               = Nchem;
DIAG.Ngrain           = Ngrain;
DIAG.Ne               = Ne;
DIAG.active_phi_ratio = Nphi/max(Ngrain*Nnode,1);
DIAG.active_cell_mean = mean(active_per_cell(:));
DIAG.active_cell_max  = max(active_per_cell(:));
DIAG.kappa_on_frac    = mean(kappa_eff(:) > 0);
DIAG.kappa_max        = max(kappa_eff(:));
DIAG.kappa_mean       = mean(kappa_eff(:));
DIAG.maskPhi          = maskPhi;
DIAG.use_Jphi         = ~isempty(Jphi_diag);

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
