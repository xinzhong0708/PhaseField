function [STATE_NEW,DIAG] = PF_Coupled_ACCH_LETangent_CS_offdiagM(STATE_REF,PARAM,MODEL,GRID,PHYS,NUM,STATE_COEF)
%PF_COUPLED_ACCH_LETANGENT
%
% Masked coupled AC-CH tangent solve for many grains.
%
% Unknowns:
%   x = [active dphi_grain ; global dmu]
%
% STATE_REF.phi, STATE_REF.p, STATE_REF.e, STATE_REF.omg are grain-sized.
% MODEL.phase_index is not used directly here; LE_Run should already expand
% e and omega back to grain size.
%
% This function does not call LE_Run.

if nargin < 7 || isempty(STATE_COEF)
    STATE_COEF = STATE_REF;
end

% Optional cache for repeated sparse column ordering.
% Matrix values change every step, but the sparsity pattern often stays the same.
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
% Coefficient state
% ------------------------------------------------------------
chi0 = STATE_COEF.chi;
e0   = STATE_COEF.e;

% Convex-split implicit capacity.
% Raw mu_e is still used as the physical driving force.
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
use_full_M     = Is_Full_M_Local(PARAM.M,Ne);

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


% ------------------------------------------------------------
% Optional composition-space kappa.
%
% If enabled, the fourth-order regularization is applied to phase
% composition c through dc = inv(H_c)*J'*dmu, instead of to total E.
% This keeps the existing [dphi; dmu] unknown structure.
% ------------------------------------------------------------
use_kappa_c = isfield(PARAM,'use_kappa_c') && PARAM.use_kappa_c == 1;

if use_kappa_c
    KAPC = Build_KappaC_Data(STATE_COEF,PARAM,MODEL,kappa_eff);
else
    KAPC = [];
end

% ------------------------------------------------------------
% Build AC source using old geometry but coefficient omega
% ------------------------------------------------------------
STATE_SRC       = STATE_REF;
STATE_SRC.mu_e  = STATE_COEF.mu_e;
STATE_SRC.chi   = STATE_COEF.chi;
STATE_SRC.e     = STATE_COEF.e;
STATE_SRC.omg   = STATE_COEF.omg;

% ------------------------------------------------------------
% AC source stabilizer
% Default: off, because it changes effective AC kinetics.
% ------------------------------------------------------------
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
if isfield(PARAM,'maskPhi_ACCH') && ~isempty(PARAM.maskPhi_ACCH)

    maskPhi = logical(PARAM.maskPhi_ACCH);

else

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

end

active_per_cell = sum(maskPhi,3);

% ------------------------------------------------------------
% Optional real diagonal AC source Jacobian
%
% Jphi_diag{alpha} = dS_alpha/dphi_alpha, with mu/omega frozen.
% This is different from A_ac: it is the true local source tangent.
% In the AC row it enters as -Jphi_diag*dphi.
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
Nmu  = Ne*Nnode;
Ntot = Nphi + Nmu;

idMu = cell(1,Ne);

for l = 1:Ne
    idMu{l} = Nphi + ((l-1)*Nnode + (1:Nnode)).';
end

% ------------------------------------------------------------
% Fast grain-aware tangent ingredients
%
% For p_i = phi_i^2 / sum_j phi_j^2:
%
%   sum_i e_i * dp_i/dphi_alpha
%     = 2*phi_alpha/D * (e_alpha - sum_i p_i e_i)
%
% This avoids building dpdphi{alpha,ip}.
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
    Nnode * Ne * (5 + 13*Ne) + ...
    13 * Ne * max(Nphi,1) + ...
    1000;

% Extra allocation for aggregated composition-space kappa block.
% The old allocation did not include these entries, causing repeated growth
% of rows/cols/vals when use_kappa_c = 1.
if use_kappa_c && ~isempty(KAPC)
    max_nnz = max_nnz + Nnode * Ne * 13 * Ne;
end

if use_full_M
    max_nnz = max_nnz + Nnode * Ne * 5 * (Ne-1);
end

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

    if isfield(PARAM,'LK_AC') && size(PARAM.LK_AC,3) >= alpha
        LK_a = PARAM.LK_AC(:,:,alpha);
    else
        LK_a = PARAM.LK;
    end

    rhs_full = LK_a .* lap_phi_ref + S_AC{alpha};
    R(row) = rhs_full(idx_a);

    LKc = LK_a(idx_a);

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

    % Real diagonal phi-Jacobian of AC source:
    %   S_alpha(phi+dphi) ~= S_alpha + Jphi*dphi
    % Move to left:
    %   -Jphi*dphi
    if ~isempty(Jphi_diag)
        Jaa = Jphi_diag{alpha};
        [rows,cols,vals,k] = add_active_block( ...
            rows,cols,vals,k,row,idPhiMap{alpha},idx_a,-Jaa(idx_a));
    end

    % Coupling to dmu:
    % coeff_mu = -L * sum_i e_i * dp_i/dphi_alpha
    if isfield(PARAM,'L_AC') && ~isempty(PARAM.L_AC) && size(PARAM.L_AC,3) >= alpha
        L_a = PARAM.L_AC(:,:,alpha);
    else
        L_a = PARAM.L;
    end

    for ie = 1:Ne

        B_ie_alpha = facPhi{alpha} .* (e0{alpha}{ie} - e_bar{ie});
        coeff_mu   = -L_a .* B_ie_alpha;

        [rows,cols,vals,k] = add_block( ...
            rows,cols,vals,k,row,idMu{ie}(idx_a),coeff_mu(idx_a));

    end

end

% ============================================================
% 2. Cahn-Hilliard block with implicit LE closure
% ============================================================
if use_kappa_c

    % --------------------------------------------------------
    % Fast composition-kappa path.
    %
    % In composition-kappa mode A_E = 1/dt only. The fourth-order
    % contribution is assembled separately by Add_KappaC_Block.
    % Therefore all 13-neighbour A_E*chi and A_E*B entries are zero
    % except the center entries. Avoid building and inserting these
    % zero neighbour blocks.
    % --------------------------------------------------------
    invdt = 1/dt;

    for l = 1:Ne

        row = idMu{l}(idx_c);

        % Diagonal mobility is still used by the fourth-order kappa term.
        Ml   = Get_Mobility_Local(PARAM,l,l,ny,nx);

        % Diffusion residual with optional off-diagonal mobility:
        %   D_lm mu_m = -div(M_lm grad(mu_m))
        D_muref = zeros(size(idx_c));
        Dcoef   = cell(Ne,5);

        for m = 1:Ne

            Mlm = Get_Mobility_Local(PARAM,l,m,ny,nx);

            if ~any(Mlm(:))
                continue
            end

            M_c  = Mlm(idx_c);
            M_L  = Mlm(idx_L);
            M_R  = Mlm(idx_R);
            M_U  = Mlm(idx_U);
            M_D  = Mlm(idx_D);

            d_L = -(M_L + M_c)/2/dx2;
            d_R = -(M_R + M_c)/2/dx2;
            d_U = -(M_U + M_c)/2/dy2;
            d_D = -(M_D + M_c)/2/dy2;
            d_C = -(d_L + d_R + d_U + d_D);

            Dcoef{m,1} = d_C;
            Dcoef{m,2} = d_L;
            Dcoef{m,3} = d_R;
            Dcoef{m,4} = d_U;
            Dcoef{m,5} = d_D;

            mu_m = mu_ref{m};

            mu_c = mu_m(idx_c);
            mu_L = mu_m(idx_L);
            mu_R = mu_m(idx_R);
            mu_U = mu_m(idx_U);
            mu_D = mu_m(idx_D);

            D_muref = D_muref + ...
                d_C .* mu_c + d_L .* mu_L + d_R .* mu_R + d_U .* mu_U + d_D .* mu_D;

        end

        % Since A_E = 1/dt, E_ref/dt - A_E*E_ref cancels exactly.
        R(row) = -D_muref;

        % Composition-space fourth-order term.
        % This adds the implicit operator from
        %   dc_phase = inv(H_c)*J'*dmu
        % and the explicit old c contribution, while p is kept fixed.
        if ~isempty(KAPC)
            [R,rows,cols,vals,k] = Add_KappaC_Block( ...
                R,rows,cols,vals,k,row,idMu,l,Ml,kappa_eff,KAPC, ...
                idx_c,idx_L,idx_R,idx_U,idx_D,idx_L2,idx_R2,idx_U2,idx_D2, ...
                idx_UR,idx_DR,idx_UL,idx_DL,dx2,dy2,dx4,dy4);
        end

        % Diffusion block on dmu_m.  Off-diagonal M only enters here.
        for m = 1:Ne

            if isempty(Dcoef{m,1})
                continue
            end

            [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_c),Dcoef{m,1});
            [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_L),Dcoef{m,2});
            [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_R),Dcoef{m,3});
            [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_U),Dcoef{m,4});
            [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_D),Dcoef{m,5});

        end

        % A_E * chi * dmu.  In c-kappa mode only the center term remains.
        for m = 1:Ne
            Chi_lm = chi_imp{l,m};
            [rows,cols,vals,k] = add_block( ...
                rows,cols,vals,k,row,idMu{m}(idx_c),invdt.*Chi_lm(idx_c));
        end

        % A_E * B * dphi.  In c-kappa mode only the center term remains.
        for alpha = 1:Ngrain

            idmap = idPhiMap{alpha};

            % If this grain has no active dphi, skip immediately.
            if Nphi > 0 && ~any(idmap(:))
                continue
            end

            B_la = facPhi{alpha} .* (e0{alpha}{l} - e_bar{l});

            [rows,cols,vals,k] = add_active_block( ...
                rows,cols,vals,k,row,idmap,idx_c,invdt.*B_la(idx_c));

        end
    end

else

    for l = 1:Ne

        row = idMu{l}(idx_c);

        % Diagonal mobility is still used by the fourth-order kappa term.
        Ml   = Get_Mobility_Local(PARAM,l,l,ny,nx);

        % Diffusion residual with optional off-diagonal mobility:
        %   D_lm mu_m = -div(M_lm grad(mu_m))
        D_muref = zeros(size(idx_c));
        Dcoef   = cell(Ne,5);

        for m = 1:Ne

            Mlm = Get_Mobility_Local(PARAM,l,m,ny,nx);

            if ~any(Mlm(:))
                continue
            end

            M_c  = Mlm(idx_c);
            M_L  = Mlm(idx_L);
            M_R  = Mlm(idx_R);
            M_U  = Mlm(idx_U);
            M_D  = Mlm(idx_D);

            d_L = -(M_L + M_c)/2/dx2;
            d_R = -(M_R + M_c)/2/dx2;
            d_U = -(M_U + M_c)/2/dy2;
            d_D = -(M_D + M_c)/2/dy2;
            d_C = -(d_L + d_R + d_U + d_D);

            Dcoef{m,1} = d_C;
            Dcoef{m,2} = d_L;
            Dcoef{m,3} = d_R;
            Dcoef{m,4} = d_U;
            Dcoef{m,5} = d_D;

            mu_m = mu_ref{m};

            mu_c = mu_m(idx_c);
            mu_L = mu_m(idx_L);
            mu_R = mu_m(idx_R);
            mu_U = mu_m(idx_U);
            mu_D = mu_m(idx_D);

            D_muref = D_muref + ...
                d_C .* mu_c + d_L .* mu_L + d_R .* mu_R + d_U .* mu_U + d_D .* mu_D;

        end

        % --------------------------------------------------------
        % A_E operator on E.
        %
        % Standard mode:
        %   A_E = 1/dt + lap(M*kappa_eff*lap(E))
        % --------------------------------------------------------
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

        E_l  = E_ref{l};

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

        AE_Eref = ...
            a_C  .* E_c  + ...
            a_L  .* E_L  + a_R  .* E_R  + a_U  .* E_U  + a_D  .* E_D  + ...
            a_L2 .* E_L2 + a_R2 .* E_R2 + a_U2 .* E_U2 + a_D2 .* E_D2 + ...
            a_UR .* E_UR + a_DR .* E_DR + a_UL .* E_UL + a_DL .* E_DL;

        R(row) = E_c/dt - AE_Eref - D_muref;

        % Diffusion block on dmu_m.  Off-diagonal M only enters here.
        for m = 1:Ne

            if isempty(Dcoef{m,1})
                continue
            end

            [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_c),Dcoef{m,1});
            [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_L),Dcoef{m,2});
            [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_R),Dcoef{m,3});
            [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_U),Dcoef{m,4});
            [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_D),Dcoef{m,5});

        end

        % A_E * chi * dmu
        for m = 1:Ne

            Chi_lm = chi_imp{l,m};

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

        % A_E * B * dphi
        for alpha = 1:Ngrain

            idmap = idPhiMap{alpha};

            % If this grain has no active dphi, skip immediately.
            if Nphi > 0 && ~any(idmap(:))
                continue
            end

            B_la = facPhi{alpha} .* (e0{alpha}{l} - e_bar{l});

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

end

% ------------------------------------------------------------
% Assemble
% ------------------------------------------------------------
rows = rows(1:k-1);
cols = cols(1:k-1);
vals = vals(1:k-1);

A = sparse(rows,cols,vals,Ntot,Ntot);



% ------------------------------------------------------------
% Solve linear system
% ------------------------------------------------------------
NUM_L = NUM;
if ~isfield(NUM_L,'linear_cache_id') || isempty(NUM_L.linear_cache_id)
    NUM_L.linear_cache_id = 'ACCH';
end
if isfield(NUM,'ilu_reuse_steps_ACCH') && ~isempty(NUM.ilu_reuse_steps_ACCH)
    NUM_L.ilu_reuse_steps = NUM.ilu_reuse_steps_ACCH;
end
[sol,LDIAG] = PF_LinearSolve(A,R,NUM_L);
flag        = LDIAG.flag;
relres      = LDIAG.relres;
iter        = LDIAG.iter;

% ------------------------------------------------------------
% Unpack dphi and dmu
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

dmu = cell(1,Ne);

for ie = 1:Ne
    dmu{ie} = reshape(sol(idMu{ie}),ny,nx);
end

% ------------------------------------------------------------
% Update phi and p from fixed reference state
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
% Update mu from fixed reference state
% ------------------------------------------------------------
STATE_NEW.mu_e = mu_ref;

for ie = 1:Ne
    STATE_NEW.mu_e{ie} = mu_ref{ie} + dmu{ie};
end

% ------------------------------------------------------------
% Recover E from LE tangent
% ------------------------------------------------------------
STATE_NEW.E = E_ref;

for ie = 1:Ne

    En = E_ref{ie};

    for je = 1:Ne
        En = En + chi_imp{ie,je} .* dmu{je};
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
DIAG.flag             = flag;
DIAG.relres           = relres;
DIAG.iter             = iter;
DIAG.linear_solver    = LDIAG.linear_solver;
DIAG.linear_cache_id  = LDIAG.linear_cache_id;
DIAG.solve_flag       = LDIAG.flag;
DIAG.solve_relres     = LDIAG.relres;
DIAG.solve_iter       = LDIAG.iter;
DIAG.solve_time       = LDIAG.solve_time;
DIAG.prec_time        = LDIAG.prec_time;
DIAG.fallback_used    = LDIAG.fallback_used;
DIAG.perm_method      = LDIAG.perm_method;
DIAG.ilu_cache_hit    = LDIAG.ilu_cache_hit;
DIAG.ilu_cache_age    = LDIAG.ilu_cache_age;
DIAG.ilu_rebuilt      = LDIAG.ilu_rebuilt;
DIAG.max_dphi         = max(abs(dphi(:)));
DIAG.KAPC             = KAPC;
max_dmu               = 0;

for ie = 1:Ne
    max_dmu = max(max_dmu,max(abs(dmu{ie}(:))));
end

DIAG.max_dmu          = max_dmu;
DIAG.matrix_size      = Ntot;
DIAG.nnz              = nnz(A);
% DIAG.order_cache_hit  = order_cache_hit;
DIAG.Nphi_active      = Nphi;
DIAG.Nphi_full        = Ngrain*Nnode;
DIAG.Nmu              = Nmu;
DIAG.Ngrain           = Ngrain;
DIAG.Ne               = Ne;
DIAG.active_phi_ratio = Nphi / max(Ngrain*Nnode,1);
DIAG.active_cell_mean = mean(active_per_cell(:));
DIAG.active_cell_max  = max(active_per_cell(:));
DIAG.kappa_on_frac    = mean(kappa_eff(:) > 0);
DIAG.kappa_max        = max(kappa_eff(:));
DIAG.kappa_mean       = mean(kappa_eff(:));

% Fields required by the equation checker
DIAG.maskPhi          = maskPhi;
DIAG.use_Jphi         = ~isempty(Jphi_diag);
DIAG.use_CS_chi       = isfield(PARAM,'use_CS_chi') && PARAM.use_CS_chi == 1;
DIAG.use_kappa_c      = use_kappa_c;
DIAG.use_full_M_diffusion = use_full_M;
DIAG.full_M_kappa_diag    = use_full_M;
DIAG.N_kappa_c       = numel(KAPC);
if isempty(KAPC)
    DIAG.kappa_c_grains = [];
else
    DIAG.kappa_c_grains = [KAPC.ig];
end

if ~isempty(Jphi_diag)

    DIAG.Jphi_diag = Jphi_diag;

    max_Jphi = 0;

    for alpha = 1:Ngrain
        max_Jphi = max(max_Jphi,max(abs(Jphi_diag{alpha}(:))));
    end

    DIAG.max_Jphi = max_Jphi;

else

    DIAG.Jphi_diag = cell(1,Ngrain);

    for alpha = 1:Ngrain
        DIAG.Jphi_diag{alpha} = zeros(ny,nx);
    end

    DIAG.max_Jphi = 0;

end

end



function KAPC = Build_KappaC_Data(STATE,PARAM,MODEL,kappa_eff)
%BUILD_KAPPAC_DATA Build local dc/dmu maps for composition-space kappa.
%
% For each grain/phase:
%   dc = inv(H_c_use)*J'*dmu
%
% p is kept fixed in the fourth-order term through the coefficient
% M*kappa_eff*p.

[ny,nx,Ng] = size(STATE.p);
N          = ny*nx;
Ne         = numel(STATE.E);

% Use ghost-extended phase composition for the kappa(c) stencil if present.
% Physical thermodynamics outside this block still uses STATE.c.
if isfield(STATE,'c_ext') && ~isempty(STATE.c_ext)
    c_kap = STATE.c_ext;
else
    c_kap = STATE.c;
end

% Use composition-space kappa for every grain/phase.
% Kappa itself can still be heterogeneous through kappa_eff.
grains = 1:Ng;

KAPC = struct('ig',{},'p',{},'c',{},'J',{},'X',{},'Nc',{});

for ii = 1:numel(grains)

    ig = grains(ii);

    p_ig = STATE.p(:,:,ig);

    p_on   = 1e-5;
    p_full = 1e-3;

    if isfield(PARAM,'kappa_c_p_on')
        p_on = PARAM.kappa_c_p_on;
    end
    if isfield(PARAM,'kappa_c_p_full')
        p_full = PARAM.kappa_c_p_full;
    end

    if p_full <= p_on
        p_full = p_on + eps;
    end

    w_p = (p_ig - p_on)./(p_full - p_on);
    w_p = min(max(w_p,0),1);
    w_p = w_p.^2.*(3 - 2*w_p);

    p_kap = p_ig .* w_p;

    if max(abs(kappa_eff(:).*p_kap(:))) == 0
        continue
    end


    try
        R = PhaseThermo(MODEL.pars{ig},c_kap{ig});
    catch
        continue
    end

    if isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    H = SymPages(R.H_c);
    J = R.Jac;

    if size(J,1) ~= Ne
        continue
    end

    Nc = size(H,1);

    Huse = Make_KappaC_Hessian(H,PARAM);
    JT   = permute(J,[2 1 3]);

    X = pagemldivide(Huse,JT);

    Cmat = zeros(Nc,N);

    for ic = 1:Nc
        Cmat(ic,:) = reshape(c_kap{ig}{ic},1,N);
    end

    ik = numel(KAPC)+1;

    KAPC(ik).ig = ig;
    KAPC(ik).p  = p_kap;
    KAPC(ik).c  = Cmat;
    KAPC(ik).J  = J;
    KAPC(ik).X  = X;
    KAPC(ik).Nc = Nc;

end

end



function Huse = Make_KappaC_Hessian(H,PARAM)
%MAKE_KAPPAC_HESSIAN Hessian used only for implicit c-kappa response.

[Nc,~,N] = size(H);

h_floor = 1e-8;
h_cap   = inf;
posdef  = 1;

if isfield(PARAM,'kappa_c_H_floor')
    h_floor = PARAM.kappa_c_H_floor;
end
if isfield(PARAM,'kappa_c_H_cap')
    h_cap = PARAM.kappa_c_H_cap;
end
if isfield(PARAM,'kappa_c_posdef')
    posdef = PARAM.kappa_c_posdef;
end

H = 0.5*(H + permute(H,[2 1 3]));

bad = squeeze(any(any(~isfinite(H),1),2));

if any(bad)
    H(:,:,bad) = repmat(eye(Nc),1,1,nnz(bad));
end

[V,D] = pageeig(H);

lam = zeros(Nc,N);

for i = 1:Nc
    lam(i,:) = reshape(D(i,i,:),1,N);
end

scale     = max(1,max(abs(lam),[],1));
floor_val = h_floor .* scale;

if posdef == 1

    lam_use = max(abs(lam),floor_val);

    if isfinite(h_cap)
        lam_use = min(lam_use,h_cap .* scale);
    end

else

    sgn = sign(lam);
    sgn(sgn == 0) = 1;

    lam_abs = max(abs(lam),floor_val);

    if isfinite(h_cap)
        lam_abs = min(lam_abs,h_cap .* scale);
    end

    lam_use = sgn .* lam_abs;

end

Duse = zeros(Nc,Nc,N);

for i = 1:Nc
    Duse(i,i,:) = reshape(lam_use(i,:),1,1,N);
end

Huse = pagemtimes(pagemtimes(V,Duse),permute(V,[2 1 3]));
Huse = 0.5*(Huse + permute(Huse,[2 1 3]));

end


function [R,rows,cols,vals,k] = Add_KappaC_Block( ...
    R,rows,cols,vals,k,row,idMu,l,Ml,kappa_eff,KAPC, ...
    idx_c,idx_L,idx_R,idx_U,idx_D,idx_L2,idx_R2,idx_U2,idx_D2, ...
    idx_UR,idx_DR,idx_UL,idx_DL,dx2,dy2,dx4,dy4)
%ADD_KAPPAC_BLOCK Add composition-space fourth-order contribution.
%
% Faster version with same operator:
%   - vectorizes the sum over phase composition variables ic;
%   - inserts each final neighbour block once per element m;
%   - skips rows where the local kappa stencil coefficient is zero.

Ne = numel(idMu);
[ny,nx] = size(kappa_eff);
Nall    = ny*nx;
Nrow    = numel(idx_c);

idx_off = {idx_c,idx_L,idx_R,idx_U,idx_D, ...
           idx_L2,idx_R2,idx_U2,idx_D2, ...
           idx_UR,idx_DR,idx_UL,idx_DL};

acc = zeros(Nrow,Ne,13);

for ik = 1:numel(KAPC)

    K_l = Ml .* kappa_eff .* KAPC(ik).p;

    if ~any(K_l(:))
        continue
    end

    [q_C,q_L,q_R,q_U,q_D,q_L2,q_R2,q_U2,q_D2,q_UR,q_DR,q_UL,q_DL] = ...
        Kappa4_Coeffs(K_l,idx_c,idx_L,idx_R,idx_U,idx_D, ...
                      idx_L2,idx_R2,idx_U2,idx_D2,idx_UR,idx_DR,idx_UL,idx_DL, ...
                      dx2,dy2,dx4,dy4);

    q_off = {q_C,q_L,q_R,q_U,q_D, ...
             q_L2,q_R2,q_U2,q_D2, ...
             q_UR,q_DR,q_UL,q_DL};

    Nc   = KAPC(ik).Nc;
    Cmat = KAPC(ik).c;              % Nc x Nall
    Xall = KAPC(ik).X;              % Nc x Ne x Nall

    % J_l is Nc x Nall
    J_l  = reshape(KAPC(ik).J(l,:,:),Nc,Nall);
    J_c0 = J_l(:,idx_c);            % Nc x Nrow

    for io = 1:13

        q = q_off{io};

        use = isfinite(q) & q ~= 0;

        if ~any(use)
            continue
        end

        idxo  = idx_off{io};
        idxou = idxo(use);
        qu    = q(use);

        J_c = J_c0(:,use);

        % Explicit old-c contribution:
        %   sum_ic J_l_ic(center) * c_ic(offset)
        JC = sum(J_c .* Cmat(:,idxou),1).';
        R(row(use)) = R(row(use)) - qu .* JC;

        % Implicit dc/dmu contribution:
        %   sum_ic J_l_ic(center) * X_icm(offset)
        Xnei = Xall(:,:,idxou);                         % Nc x Ne x Nuse
        J3   = reshape(J_c,Nc,1,nnz(use));              % Nc x 1  x Nuse
        JX   = reshape(sum(J3 .* Xnei,1),Ne,nnz(use)).'; % Nuse x Ne

        acc(use,:,io) = acc(use,:,io) + qu .* JX;

    end

end

for m = 1:Ne

    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_c),  acc(:,m,1));
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_L),  acc(:,m,2));
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_R),  acc(:,m,3));
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_U),  acc(:,m,4));
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_D),  acc(:,m,5));

    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_L2), acc(:,m,6));
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_R2), acc(:,m,7));
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_U2), acc(:,m,8));
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_D2), acc(:,m,9));

    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_UR), acc(:,m,10));
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_DR), acc(:,m,11));
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_UL), acc(:,m,12));
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_DL), acc(:,m,13));

end

end



function [q_C,q_L,q_R,q_U,q_D,q_L2,q_R2,q_U2,q_D2,q_UR,q_DR,q_UL,q_DL] = ...
    Kappa4_Coeffs(K_l,idx_c,idx_L,idx_R,idx_U,idx_D, ...
                  idx_L2,idx_R2,idx_U2,idx_D2,idx_UR,idx_DR,idx_UL,idx_DL, ...
                  dx2,dy2,dx4,dy4)
%KAPPA4_COEFFS Variable-coefficient lap(K*lap(field)) coefficients.

cx4 = 1/dx4;
cy4 = 1/dy4;
cxy = 1/(dx2*dy2);

K_c = K_l(idx_c);
K_L = K_l(idx_L);
K_R = K_l(idx_R);
K_U = K_l(idx_U);
K_D = K_l(idx_D);

q_L  = -(2*cx4 + 2*cxy) .* (K_L + K_c);
q_R  = -(2*cx4 + 2*cxy) .* (K_R + K_c);
q_U  = -(2*cy4 + 2*cxy) .* (K_U + K_c);
q_D  = -(2*cy4 + 2*cxy) .* (K_D + K_c);

q_L2 = cx4 .* K_L;
q_R2 = cx4 .* K_R;
q_U2 = cy4 .* K_U;
q_D2 = cy4 .* K_D;

q_UR = cxy .* (K_U + K_R);
q_DR = cxy .* (K_D + K_R);
q_UL = cxy .* (K_U + K_L);
q_DL = cxy .* (K_D + K_L);

q_C  = cx4 .* (K_L + K_R + 4*K_c) + ...
       cy4 .* (K_U + K_D + 4*K_c) + ...
       8*cxy .* K_c;

end



function H = SymPages(H)

H = 0.5*(H+permute(H,[2 1 3]));

end



function chi_imp = ConvexifyChi_ForImplicit(chi_raw,PARAM)
%CONVEXIFYCHI_FORIMPLICIT Positive bounded capacity for implicit CH terms.
%
% Faster serial/vectorized version.
%
% Main speedup:
%   pack chi_raw into a numeric Ne x Ne x N array first, so the loop over
%   grid points does not repeatedly access cell arrays.
%
% If pageeig exists in your MATLAB version, it uses page-wise eig.

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

% ------------------------------------------------------------
% Pack cell chi into page array C(:,:,id)
% ------------------------------------------------------------
Cpage = zeros(Ne,Ne,N);

for i = 1:Ne
    for j = 1:Ne
        Cpage(i,j,:) = reshape(chi_raw{i,j},1,1,N);
    end
end

% Symmetrize exactly
Cpage = 0.5*(Cpage + permute(Cpage,[2 1 3]));

% ------------------------------------------------------------
% True page-wise eig path, if available
% ------------------------------------------------------------
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


% ------------------------------------------------------------
% Unpack back to cell
% ------------------------------------------------------------
chi_imp = cell(Ne,Ne);

for i = 1:Ne
    for j = 1:Ne
        chi_imp{i,j} = reshape(Cimp_page(i,j,:),ny,nx);
    end
end

end



function [q,hit,CACHE] = CachedColamdOrdering(CACHE,A,rows,cols,Ntot)
%CACHEDCOLAMDORDERING Reuse COLAMD ordering if sparse pattern is unchanged.
%
% The numerical values of A change every step, but COLAMD depends only on
% the sparsity pattern.  This avoids repeating the symbolic ordering when
% the active phi mask and CH pattern are unchanged.

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
%CALC_JPHI_DIAG_FD
%
% Diagonal real AC-source Jacobian:
%
%   Jdiag{alpha} = dS_alpha/dphi_alpha
%
% omega/mu are frozen. This is not a stabilizer. It is the actual local
% derivative of Calc_S_AllenCahn with respect to the same phi component.
%
% The derivative is computed only on active dphi rows.

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




function tf = Is_Full_M_Local(M,Ne)

tf = iscell(M) && size(M,1) >= Ne && size(M,2) >= Ne;

end


function Mlm = Get_Mobility_Local(PARAM,l,m,ny,nx)
%GET_MOBILITY_LOCAL Return mobility field M_lm.
%
% Supports old diagonal format:
%   PARAM.M{l}
% and new full matrix format:
%   PARAM.M{l,m}

use_full_M = iscell(PARAM.M) && size(PARAM.M,1) >= l && size(PARAM.M,2) >= m && size(PARAM.M,1) > 1 && size(PARAM.M,2) > 1;

if use_full_M

    if size(PARAM.M,1) >= l && size(PARAM.M,2) >= m
        Mlm = PARAM.M{l,m};
    else
        Mlm = [];
    end

else

    if l == m
        Mlm = PARAM.M{l};
    else
        Mlm = [];
    end

end

if isempty(Mlm)
    Mlm = zeros(ny,nx);
elseif isscalar(Mlm)
    Mlm = Mlm*ones(ny,nx);
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