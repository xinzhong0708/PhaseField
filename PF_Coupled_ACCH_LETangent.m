function [STATE_NEW,DIAG] = PF_Coupled_ACCH_LETangent(STATE_REF,PARAM,MODEL,GRID,PHYS,NUM,STATE_COEF)
%PF_COUPLED_ACCH_LETANGENT
%
% Masked coupled AC-CH tangent solve.
%
% dmu is solved globally on all grid cells.
% dphi is solved only where the corresponding phi field is active.
%
% STATE_REF:
%   fixed old accepted state at time t^n.
%
% STATE_COEF:
%   optional coefficient state for chi/e/omega/mu_e.
%   If omitted, STATE_REF is used.
%
% Unknowns:
%   x = [active dphi ; global dmu]
%
% LE tangent:
%   dE = chi_coef*dmu + sum_i e_i_coef*dp_i
%
% This function does not call LE_Run.

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
[ny,nx,Nop] = size(phi_ref);      % number of order parameters
Nphase      = size(p_ref,3);      % number of thermodynamic phases
Ne          = numel(E_ref);
Nnode       = nx*ny;

dt = NUM.dt_phy;
dx = GRID.dx;
dy = GRID.dy;

dx2 = dx^2;
dy2 = dy^2;
dx4 = dx2^2;
dy4 = dy2^2;

% ------------------------------------------------------------
% Build AC source using old geometry but coefficient omega
% ------------------------------------------------------------
STATE_SRC       = STATE_REF;
STATE_SRC.mu_e  = STATE_COEF.mu_e;
STATE_SRC.chi   = STATE_COEF.chi;
STATE_SRC.e     = STATE_COEF.e;
STATE_SRC.omg   = STATE_COEF.omg;

PARAM.A_ac = 0*Calc_Aac_FrozenOmega(STATE_SRC,PARAM,MODEL,3,1e-6,0,[]);
STATE_SRC  = Calc_S_AllenCahn(STATE_SRC,PARAM,MODEL);
S_AC       = STATE_SRC.S_AC;

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

% % ------------------------------------------------------------
% % Active phi mask
% % ------------------------------------------------------------
% if isfield(NUM,'phi_mask_cut')
%     phi_cut = NUM.phi_mask_cut;
% else
%     phi_cut = 1e-8;
% end
% 
% if isfield(NUM,'phi_mask_thick')
%     mask_thick = NUM.phi_mask_thick;
% else
%     if isfield(PHYS,'kappa') && PHYS.kappa ~= 0
%         mask_thick = 2;
%     else
%         mask_thick = 1;
%     end
% end
% 
% % Remove tiny numerical tails before masking.
% % This prevents all phases/grains becoming active everywhere.
% phi_for_mask = phi_ref;
% phi_for_mask(phi_for_mask < phi_cut) = 0;
% 
% if exist('Calc_Mask','file') == 2
%     maskPhi = Calc_Mask(phi_for_mask,mask_thick) > 0;
% else
%     maskPhi = Local_Calc_Mask(phi_for_mask,phi_cut,mask_thick);
% end
% 
% maskPhi = logical(maskPhi);
% active_per_cell = sum(maskPhi,3);



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
    if isfield(PHYS,'kappa') && PHYS.kappa ~= 0
        mask_thick = 2;
    else
        mask_thick = 1;
    end
end

% Interface-only mask:
%   do not solve where phi is fully absent
%   do not solve where phi is fully pure
%   solve only near 0 < p_phi < 1 and its stencil neighborhood
maskPhi = Local_Calc_Interface_Mask(phi_ref,phi_cut,pure_cut,mask_thick);

% Optional safety: include cells where AC source is unexpectedly nonzero.
% Usually keep this off unless debugging.
if isfield(NUM,'phi_mask_source_tol') && ~isempty(NUM.phi_mask_source_tol)

    source_tol = NUM.phi_mask_source_tol;

    for alpha = 1:Nop

        core_source = abs(S_AC{alpha}) > source_tol;
        mask_source = Local_Dilate_Mask(core_source,mask_thick);

        maskPhi(:,:,alpha) = maskPhi(:,:,alpha) | mask_source;

    end

end

maskPhi = logical(maskPhi);
active_per_cell = sum(maskPhi,3);


% ------------------------------------------------------------
% Unknown ids
%
% dphi: active mask only
% dmu : global everywhere
% ------------------------------------------------------------
idPhiMap = cell(1,Nop);

counter = 0;

for alpha = 1:Nop

    idmap = zeros(ny,nx);

    ids  = find(maskPhi(:,:,alpha));
    nids = numel(ids);

    idmap(ids) = counter + (1:nids);

    idPhiMap{alpha} = idmap;
    counter = counter + nids;
end

Nphi = counter;
Nmu  = Ne*Nnode;
Ntot = Nphi + Nmu;

idMu = cell(1,Ne);

for l = 1:Ne
    idMu{l} = Nphi + ((l-1)*Nnode + (1:Nnode)).';
end

% ------------------------------------------------------------
% Tangents dp_i/dphi_alpha at fixed old geometry
% ------------------------------------------------------------
dpdphi = cell(Nop,Nphase);

for alpha = 1:Nop
    for ip = 1:Nphase
        dpdphi{alpha,ip} = MODEL.dpdphi(alpha,ip,phi_ref);
    end
end

% ------------------------------------------------------------
% B{ie,alpha} = dE_ie/dphi_alpha
%              = sum_ip e_ip_ie * dp_ip/dphi_alpha
% ------------------------------------------------------------
B = cell(Ne,Nop);

for ie = 1:Ne
    for alpha = 1:Nop

        tmp = zeros(ny,nx);

        for ip = 1:Nphase
            tmp = tmp + e0{ip}{ie} .* dpdphi{alpha,ip};
        end

        B{ie,alpha} = tmp;
    end
end

% ------------------------------------------------------------
% Sparse matrix allocation
% ------------------------------------------------------------
max_nnz = ...
    Nphi * (5 + Ne) + ...
    Nnode * Ne * (5 + 13*Ne) + ...
    13 * Ne * max(Nphi,1) + ...
    1000;

rows = zeros(max_nnz,1);
cols = zeros(max_nnz,1);
vals = zeros(max_nnz,1);
R    = zeros(Ntot,1);

k = 1;

% ============================================================
% 1. Allen-Cahn block
% ============================================================
for alpha = 1:Nop

    mask_a = maskPhi(:,:,alpha);
    idx_a  = find(mask_a);

    if isempty(idx_a)
        continue
    end

    row = idPhiMap{alpha}(idx_a);

    phi_a = phi_ref(:,:,alpha);

    % Reference Laplacian for RHS
    lap_phi_ref = laplacian_reflect(phi_a,dx,dy);

    rhs_full = PARAM.LK .* lap_phi_ref + S_AC{alpha};
    R(row) = rhs_full(idx_a);

    LKc = PARAM.LK(idx_a);
    Aac = PARAM.A_ac(idx_a);

    % dphi_alpha block:
    % (1/dt)dphi - LK*Lap(dphi) + Aac*dphi
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

    % Coupling to dmu through domega_i = -e_i*dmu.
    for ie = 1:Ne

        coeff_mu = zeros(ny,nx);

        for ip = 1:Nphase
            coeff_mu = coeff_mu - PARAM.L .* dpdphi{alpha,ip} .* e0{ip}{ie};
        end

        [rows,cols,vals,k] = add_block( ...
            rows,cols,vals,k,row,idMu{ie}(idx_a),coeff_mu(idx_a));
    end
end

% ============================================================
% 2. Cahn-Hilliard block with implicit LE closure
% ============================================================
for l = 1:Ne

    row = idMu{l}(idx_c);

    Ml   = PARAM.M{l};
    M_c  = Ml(idx_c);
    M_L  = Ml(idx_L);
    M_R  = Ml(idx_R);
    M_U  = Ml(idx_U);
    M_D  = Ml(idx_D);

    % Diffusion operator on dmu_l: -div(M grad dmu_l)
    d_L = -(M_L + M_c)/2/dx2;
    d_R = -(M_R + M_c)/2/dx2;
    d_U = -(M_U + M_c)/2/dy2;
    d_D = -(M_D + M_c)/2/dy2;
    d_C = -(d_L + d_R + d_U + d_D);

    % A_E operator on E:
    % A_E = 1/dt + M*kappa*biLaplacian
    q_L  = M_c * PHYS.kappa .* (-4/dx4 - 4/(dx2*dy2));
    q_R  = M_c * PHYS.kappa .* (-4/dx4 - 4/(dx2*dy2));
    q_U  = M_c * PHYS.kappa .* (-4/dy4 - 4/(dx2*dy2));
    q_D  = M_c * PHYS.kappa .* (-4/dy4 - 4/(dx2*dy2));

    q_L2 = M_c * PHYS.kappa .* (1/dx4);
    q_R2 = M_c * PHYS.kappa .* (1/dx4);
    q_U2 = M_c * PHYS.kappa .* (1/dy4);
    q_D2 = M_c * PHYS.kappa .* (1/dy4);

    q_UR = M_c * PHYS.kappa .* (2/(dx2*dy2));
    q_DR = M_c * PHYS.kappa .* (2/(dx2*dy2));
    q_UL = M_c * PHYS.kappa .* (2/(dx2*dy2));
    q_DL = M_c * PHYS.kappa .* (2/(dx2*dy2));

    q_C  = M_c * PHYS.kappa .* (6/dx4 + 6/dy4 + 8/(dx2*dy2));

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

    % --------------------------------------------------------
    % RHS uses fixed reference state
    % --------------------------------------------------------
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
        a_C  .* E_c  + ...
        a_L  .* E_L  + a_R  .* E_R  + a_U  .* E_U  + a_D  .* E_D  + ...
        a_L2 .* E_L2 + a_R2 .* E_R2 + a_U2 .* E_U2 + a_D2 .* E_D2 + ...
        a_UR .* E_UR + a_DR .* E_DR + a_UL .* E_UL + a_DL .* E_DL;

    D_muref = ...
        d_C .* mu_c + d_L .* mu_L + d_R .* mu_R + d_U .* mu_U + d_D .* mu_D;

    R(row) = E_c/dt - AE_Eref - D_muref;

    % --------------------------------------------------------
    % Diffusion block on dmu_l
    % --------------------------------------------------------
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_c),d_C);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_L),d_L);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_R),d_R);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_U),d_U);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_D),d_D);

    % --------------------------------------------------------
    % A_E * chi * dmu
    % --------------------------------------------------------
    for m = 1:Ne

        Chi_lm = chi0{l,m};

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

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_c),  a_C  .* Chi_c);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_L),  a_L  .* Chi_L);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_R),  a_R  .* Chi_R);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_U),  a_U  .* Chi_U);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_D),  a_D  .* Chi_D);

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_L2), a_L2 .* Chi_L2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_R2), a_R2 .* Chi_R2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_U2), a_U2 .* Chi_U2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_D2), a_D2 .* Chi_D2);

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_UR), a_UR .* Chi_UR);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_DR), a_DR .* Chi_DR);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_UL), a_UL .* Chi_UL);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_DL), a_DL .* Chi_DL);
    end

    % --------------------------------------------------------
    % A_E * B * dphi
    %
    % dphi columns exist only inside phi mask.
    % --------------------------------------------------------
    for alpha = 1:Nop

        B_la = B{l,alpha};

        B_c  = B_la(idx_c);
        B_L  = B_la(idx_L);
        B_R  = B_la(idx_R);
        B_U  = B_la(idx_U);
        B_D  = B_la(idx_D);

        B_L2 = B_la(idx_L2);
        B_R2 = B_la(idx_R2);
        B_U2 = B_la(idx_U2);
        B_D2 = B_la(idx_D2);

        B_UR = B_la(idx_UR);
        B_DR = B_la(idx_DR);
        B_UL = B_la(idx_UL);
        B_DL = B_la(idx_DL);

        idmap = idPhiMap{alpha};

        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_c,  a_C  .* B_c);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_L,  a_L  .* B_L);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_R,  a_R  .* B_R);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_U,  a_U  .* B_U);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_D,  a_D  .* B_D);

        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_L2, a_L2 .* B_L2);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_R2, a_R2 .* B_R2);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_U2, a_U2 .* B_U2);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_D2, a_D2 .* B_D2);

        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_UR, a_UR .* B_UR);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_DR, a_DR .* B_DR);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_UL, a_UL .* B_UL);
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_DL, a_DL .* B_DL);
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
% Solve
% ------------------------------------------------------------
direct_mode = 'colamd';

if isfield(NUM,'direct_mode')
    direct_mode = NUM.direct_mode;
end

switch lower(direct_mode)

    case 'none'

        sol = A \ R;

    case 'colamd'

        q = colamd(A);
        y = A(:,q) \ R;

        sol = zeros(size(R));
        sol(q) = y;

    case 'symamd'

        S = spones(A) + spones(A');
        p = symamd(S);

        y = A(p,p) \ R(p);

        sol = zeros(size(R));
        sol(p) = y;

    case 'symrcm'

        S = spones(A) + spones(A');
        p = symrcm(S);

        y = A(p,p) \ R(p);

        sol = zeros(size(R));
        sol(p) = y;

    otherwise

        error('Unknown NUM.direct_mode: %s',direct_mode)

end

flag   = 0;
relres = norm(A*sol - R)/max(norm(R),eps);
iter   = [0 0];




% ------------------------------------------------------------
% Unpack dphi and dmu
% ------------------------------------------------------------
dphi = zeros(ny,nx,Nop);

for alpha = 1:Nop

    idmap  = idPhiMap{alpha};
    mask_a = maskPhi(:,:,alpha);

    ids = idmap(mask_a);

    tmp = zeros(ny,nx);

    if ~isempty(ids)
        tmp(mask_a) = sol(ids);
    end

    dphi(:,:,alpha) = tmp;
end

dmu = cell(1,Ne);

for ie = 1:Ne
    dmu{ie} = reshape(sol(idMu{ie}),ny,nx);
end

% ------------------------------------------------------------
% Update phi and p from fixed reference state
% ------------------------------------------------------------
STATE_NEW          = STATE_REF;
STATE_NEW.phi      = phi_ref + dphi;
if NUM.norm_phi == 1
    STATE_NEW.phi  = Normalize_Phi_L2_Local(STATE_NEW.phi);
end
if NUM.cut_phi == 1
    STATE_NEW.phi  = Cut_Phi_Local(STATE_NEW.phi);
end
STATE_NEW.p        = Calc_p(MODEL,STATE_NEW.phi);

% ------------------------------------------------------------
% Update mu from fixed reference state
% ------------------------------------------------------------
STATE_NEW.mu_e = mu_ref;

for ie = 1:Ne
    STATE_NEW.mu_e{ie} = mu_ref{ie} + dmu{ie};
end

% ------------------------------------------------------------
% Recover E from LE tangent:
%
%   E_new = E_ref + chi_coef*dmu + sum_i e_i_coef*dp_i
% ------------------------------------------------------------
STATE_NEW.E = E_ref;

for ie = 1:Ne

    En = E_ref{ie};

    % chi*dmu
    for je = 1:Ne
        En = En + chi0{ie,je} .* dmu{je};
    end

    % sum_i e_i*dp_i
    for ip = 1:Nphase
        dp_ip = STATE_NEW.p(:,:,ip) - p_ref(:,:,ip);
        En = En + e0{ip}{ie} .* dp_ip;
    end

    STATE_NEW.E{ie} = En;
end

% Keep global mass exactly conserved relative to fixed old state
if NUM.norm_E==1
    STATE_NEW.E = EnforceMeanE_Local(STATE_NEW.E,E_ref);
end

% ------------------------------------------------------------
% Predict omega for diagnostics only.
% Final LE_Run should overwrite this.
% ------------------------------------------------------------
STATE_NEW.omg = STATE_REF.omg;

for ip = 1:Nphase

    domega = zeros(ny,nx);

    for ie = 1:Ne
        domega = domega - e0{ip}{ie} .* dmu{ie};
    end

    STATE_NEW.omg(:,:,ip) = STATE_REF.omg(:,:,ip) + domega;
end

% Remove common omega offset
omg_mean = mean(STATE_NEW.omg,3);

for ip = 1:Nphase
    STATE_NEW.omg(:,:,ip) = STATE_NEW.omg(:,:,ip) - omg_mean;
end

% chi and e are coefficient tangents until final nonlinear LE
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

for ie = 1:Ne
    max_dmu = max(max_dmu,max(abs(dmu{ie}(:))));
end

DIAG.max_dmu          = max_dmu;
DIAG.matrix_size      = Ntot;
DIAG.nnz              = nnz(A);
DIAG.Nphi_active      = Nphi;
DIAG.Nphi_full        = Nop*Nnode;
DIAG.Nmu              = Nmu;
DIAG.active_phi_ratio = Nphi / max(Nop*Nnode,1);
DIAG.active_cell_mean = mean(active_per_cell(:));
DIAG.active_cell_max  = max(active_per_cell(:));

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


function mask = Local_Calc_Mask(phi,phi_cut,thickness)

[ny,nx,Nop] = size(phi);
mask = false(ny,nx,Nop);

ker = ones(2*thickness+1,2*thickness+1);

for ip = 1:Nop
    core = phi(:,:,ip) > phi_cut;
    mask(:,:,ip) = conv2(double(core),ker,'same') > 0;
end

end


function [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_col,coeff)
%ADD_ACTIVE_BLOCK
% Add sparse entries only where the phi unknown exists.
%
% idmap(idx_col) == 0 means dphi is outside the active mask and is
% treated as zero correction.

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

function [phi] = Cut_Phi_Local(phi)
phi(phi<0) = 0;
phi(phi>1) = 1;
end


function phi = Normalize_Phi_L2_Local(phi)
%NORMALIZE_PHI_L2_LOCAL
% Natural normalization for p_i = phi_i^2 / sum_j(phi_j^2).
% If all phi are zero at a grid, keep them zero.
% Remove numerical overshoot
phi  = max(phi,0);
phi  = min(phi,1);
% L2 amplitude
s    = sqrt(sum(phi.^2,3));
% Only normalize where amplitude is nonzero
mask = s > eps;
for ip = 1:size(phi,3)
    tmp         = phi(:,:,ip);
    tmp(mask)   = tmp(mask) ./ s(mask);
    phi(:,:,ip) = tmp;
end
end

function mask = Local_Calc_Interface_Mask(phi,low_cut,pure_cut,thickness)
%LOCAL_CALC_INTERFACE_MASK
%
% Build mask for solving dphi only around interfaces.
%
% A cell is active for phi_alpha if the local normalized phase weight
%
%       q_alpha = phi_alpha^2 / sum_beta phi_beta^2
%
% is neither absent nor pure:
%
%       low_cut < q_alpha < 1 - pure_cut
%
% Sharp 0/1 interfaces are also detected by neighbor jumps in q_alpha.
% The resulting core is then dilated by 'thickness' cells.

[ny,nx,Nop] = size(phi);

mask = false(ny,nx,Nop);

den = sum(phi.^2,3) + eps;

for alpha = 1:Nop

    q = phi(:,:,alpha).^2 ./ den;

    % Diffuse-interface cells
    core = q > low_cut & q < 1 - pure_cut;

    % Also detect sharp 0/1 jumps, important for initial maps
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

    % Dilation gives stencil support
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