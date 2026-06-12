function [STATE_NEW,DIAG] = PF_Coupled_ACCH_LETangent_EPrimary_Kappa0(STATE_REF,PARAM,MODEL,GRID,PHYS,NUM,STATE_COEF)
%PF_COUPLED_ACCH_LETANGENT_EPRIMARY_KAPPA0
%
% Masked coupled AC-CH tangent solve with E as primary chemical unknown.
%
% This is a kappa=0 test version.
%
% Unknowns:
%   x = [active dphi_grain ; global dE]
%
% Local tangent:
%   dE = chi*dmu + B*dphi
%
% Therefore:
%   dmu = xi*(dE - B*dphi), where xi = inv(chi)
%
% CH row:
%   dE/dt + D*dmu = -D*mu_ref
%
% AC row:
%   AC(dphi) - L*B*dmu = RHS
%
% This avoids reconstructing E after solving dmu. E is solved directly.

if nargin < 7 || isempty(STATE_COEF)
    STATE_COEF = STATE_REF;
end

% ------------------------------------------------------------
% Fixed time-reference state
% ------------------------------------------------------------
phi_ref = STATE_REF.phi;
p_ref   = STATE_REF.p;
E_ref   = STATE_REF.E;
mu_ref  = STATE_REF.mu_e;

% ------------------------------------------------------------
% Coefficient state
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

% ------------------------------------------------------------
% Invert chi locally
% ------------------------------------------------------------
if isfield(PARAM,'xi_floor')
    xi_floor = PARAM.xi_floor;
else
    xi_floor = 1e-10;
end

xi0 = InvertChi_Local(chi0,xi_floor);

% ------------------------------------------------------------
% Build AC source using old geometry but coefficient omega
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

idx_c = sub2ind([ny,nx], ii,  jj);
idx_L = sub2ind([ny,nx], ii,  jjL);
idx_R = sub2ind([ny,nx], ii,  jjR);
idx_U = sub2ind([ny,nx], iiU, jj);
idx_D = sub2ind([ny,nx], iiD, jj);

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
    mask_thick = 1;
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
    counter = counter + nids;
end

Nphi = counter;
NEunk = Ne*Nnode;
Ntot = Nphi + NEunk;

idE = cell(1,Ne);

for ie = 1:Ne
    idE{ie} = Nphi + ((ie-1)*Nnode + (1:Nnode)).';
end

% ------------------------------------------------------------
% Grain-aware tangent ingredients B_ie_alpha
%
% B_ie_alpha = sum_i e_i * dp_i/dphi_alpha
%            = 2*phi_alpha/Dphi * (e_alpha - e_bar)
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
        tmp = tmp + p_tan(:,:,ig).*e0{ig}{ie};
    end

    e_bar{ie} = tmp;
end

B = cell(Ne,Ngrain);

for ie = 1:Ne
    for alpha = 1:Ngrain
        B{ie,alpha} = facPhi{alpha} .* (e0{alpha}{ie} - e_bar{ie});
    end
end

% ------------------------------------------------------------
% Sparse matrix allocation
% ------------------------------------------------------------
max_nnz = ...
    Nphi * (6 + Ne + Ngrain) + ...
    Nnode * Ne * (5*Ne + 5*Ngrain + 1) + ...
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

    % --------------------------------------------------------
    % Substitute dmu = xi*(dE - B*dphi) into AC coupling:
    %
    % old AC coupling:
    %   sum_i (-L*B_i_alpha) * dmu_i
    %
    % new dE coupling:
    %   sum_i,m (-L*B_i_alpha*xi_i_m) * dE_m
    %
    % new dphi coupling:
    %   -sum_i,m,beta (-L*B_i_alpha*xi_i_m*B_m_beta) * dphi_beta
    % --------------------------------------------------------

    % Coupling to dE
    for m = 1:Ne

        coeff_E = zeros(ny,nx);

        for ie = 1:Ne
            coeff_E = coeff_E - PARAM.L .* B{ie,alpha} .* xi0{ie,m};
        end

        [rows,cols,vals,k] = add_block( ...
            rows,cols,vals,k,row,idE{m}(idx_a),coeff_E(idx_a));
    end

    % Extra local dphi coupling from dmu(dphi)
    for beta = 1:Ngrain

        idmap_beta = idPhiMap{beta};

        if ~any(idmap_beta(:))
            continue
        end

        coeff_phi = zeros(ny,nx);

        for ie = 1:Ne
            for m = 1:Ne
                coeff_phi = coeff_phi + PARAM.L .* B{ie,alpha} .* xi0{ie,m} .* B{m,beta};
            end
        end

        [rows,cols,vals,k] = add_active_block( ...
            rows,cols,vals,k,row,idmap_beta,idx_a,coeff_phi(idx_a));
    end
end

% ============================================================
% 2. Cahn-Hilliard block, E-primary, kappa = 0
% ============================================================
for l = 1:Ne

    row = idE{l}(idx_c);

    Ml  = PARAM.M{l};
    M_c = Ml(idx_c);
    M_L = Ml(idx_L);
    M_R = Ml(idx_R);
    M_U = Ml(idx_U);
    M_D = Ml(idx_D);

    d_L = -(M_L + M_c)/2/dx2;
    d_R = -(M_R + M_c)/2/dx2;
    d_U = -(M_U + M_c)/2/dy2;
    d_D = -(M_D + M_c)/2/dy2;
    d_C = -(d_L + d_R + d_U + d_D);

    mu_l = mu_ref{l};

    mu_c = mu_l(idx_c);
    mu_L = mu_l(idx_L);
    mu_R = mu_l(idx_R);
    mu_U = mu_l(idx_U);
    mu_D = mu_l(idx_D);

    D_muref = d_C.*mu_c + d_L.*mu_L + d_R.*mu_R + d_U.*mu_U + d_D.*mu_D;

    % Increment equation:
    %   dE_l/dt + D*dmu_l = -D*mu_ref_l
    R(row) = -D_muref;

    % dE_l / dt
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{l}(idx_c),ones(Nnode,1)/dt);

    % D * xi * dE
    for m = 1:Ne

        Xi_lm = xi0{l,m};

        Xi_c = Xi_lm(idx_c);
        Xi_L = Xi_lm(idx_L);
        Xi_R = Xi_lm(idx_R);
        Xi_U = Xi_lm(idx_U);
        Xi_D = Xi_lm(idx_D);

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{m}(idx_c),d_C.*Xi_c);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{m}(idx_L),d_L.*Xi_L);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{m}(idx_R),d_R.*Xi_R);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{m}(idx_U),d_U.*Xi_U);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idE{m}(idx_D),d_D.*Xi_D);
    end

    % -D * xi * B * dphi
    for beta = 1:Ngrain

        idmap_beta = idPhiMap{beta};

        if ~any(idmap_beta(:))
            continue
        end

        C_lb = zeros(ny,nx);

        for m = 1:Ne
            C_lb = C_lb + xi0{l,m} .* B{m,beta};
        end

        C_c = C_lb(idx_c);
        C_L = C_lb(idx_L);
        C_R = C_lb(idx_R);
        C_U = C_lb(idx_U);
        C_D = C_lb(idx_D);

        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap_beta,idx_c,-d_C.*C_c);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap_beta,idx_L,-d_L.*C_L);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap_beta,idx_R,-d_R.*C_R);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap_beta,idx_U,-d_U.*C_U);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap_beta,idx_D,-d_D.*C_D);
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
y = A(:,q) \ R;

sol = zeros(size(R));
sol(q) = y;

flag   = 0;
relres = norm(A*sol - R)/max(norm(R),eps);
iter   = [0 0];

% ------------------------------------------------------------
% Unpack dphi and dE
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

dE = cell(1,Ne);

for ie = 1:Ne
    dE{ie} = reshape(sol(idE{ie}),ny,nx);
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
% Update E directly
% ------------------------------------------------------------
STATE_NEW.E = E_ref;

for ie = 1:Ne
    STATE_NEW.E{ie} = E_ref{ie} + dE{ie};
end

if isfield(NUM,'norm_E') && NUM.norm_E == 1
    STATE_NEW.E = EnforceMeanE_Local(STATE_NEW.E,E_ref);
end

% ------------------------------------------------------------
% Recover dmu from tangent:
%
%   dmu = xi*(dE - B*dphi)
%
% Use the solved dphi, not nonlinear dp.
% Final LE_Run should overwrite mu_e anyway.
% ------------------------------------------------------------
qE = cell(1,Ne);

for m = 1:Ne

    tmp = dE{m};

    for alpha = 1:Ngrain
        tmp = tmp - B{m,alpha} .* dphi(:,:,alpha);
    end

    qE{m} = tmp;
end

dmu = cell(1,Ne);

for ie = 1:Ne

    tmp = zeros(ny,nx);

    for m = 1:Ne
        tmp = tmp + xi0{ie,m} .* qE{m};
    end

    dmu{ie} = tmp;
end

STATE_NEW.mu_e = mu_ref;

for ie = 1:Ne
    STATE_NEW.mu_e{ie} = mu_ref{ie} + dmu{ie};
end

% ------------------------------------------------------------
% Predict omega for diagnostics only.
% Final LE_Run should overwrite this.
% ------------------------------------------------------------
STATE_NEW.omg = STATE_REF.omg;

for ig = 1:Ngrain

    domega = zeros(ny,nx);

    for ie = 1:Ne
        domega = domega - e0{ig}{ie} .* dmu{ie};
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
DIAG.flag     = flag;
DIAG.relres   = relres;
DIAG.iter     = iter;
DIAG.max_dphi = max(abs(dphi(:)));

max_dE = 0;
max_dmu = 0;

for ie = 1:Ne
    max_dE  = max(max_dE,max(abs(dE{ie}(:))));
    max_dmu = max(max_dmu,max(abs(dmu{ie}(:))));
end

DIAG.max_dE           = max_dE;
DIAG.max_dmu          = max_dmu;
DIAG.matrix_size      = Ntot;
DIAG.nnz              = nnz(A);
DIAG.Nphi_active      = Nphi;
DIAG.Nphi_full        = Ngrain*Nnode;
DIAG.NE               = NEunk;
DIAG.Ngrain           = Ngrain;
DIAG.Ne               = Ne;
DIAG.active_phi_ratio = Nphi / max(Ngrain*Nnode,1);
DIAG.active_cell_mean = mean(active_per_cell(:));
DIAG.active_cell_max  = max(active_per_cell(:));
DIAG.maskPhi          = maskPhi;
DIAG.use_Jphi         = ~isempty(Jphi_diag);

end


function xi = InvertChi_Local(chi,chi_floor)

Ne = size(chi,1);
[ny,nx] = size(chi{1,1});
N = ny*nx;

xi = cell(Ne,Ne);

for i = 1:Ne
    for j = 1:Ne
        xi{i,j} = zeros(ny,nx);
    end
end

for id = 1:N

    C = zeros(Ne,Ne);

    for i = 1:Ne
        for j = 1:Ne
            tmp = chi{i,j};
            C(i,j) = tmp(id);
        end
    end

    C = 0.5*(C+C.');

    [V,D] = eig(C);
    lam = diag(D);

    sc = max(1,max(abs(lam)));
    lam = max(lam,chi_floor*sc);

    X = V*diag(1./lam)*V.';
    X = 0.5*(X+X.');

    for i = 1:Ne
        for j = 1:Ne
            tmp = xi{i,j};
            tmp(id) = X(i,j);
            xi{i,j} = tmp;
        end
    end
end

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
    tmp(mask) = tmp(mask)./s(mask);
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