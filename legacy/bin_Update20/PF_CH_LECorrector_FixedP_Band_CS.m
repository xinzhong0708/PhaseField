function [STATE_CORR,DIAG] = PF_CH_LECorrector_FixedP_Band_CS(STATE_TIME,STATE_IT,PARAM,MODEL,GRID,PHYS,NUM)
%PF_CH_LECORRECTOR_FIXEDP_BAND Reduced-band chemical corrector at fixed updated phi and p.
%
% This is a validation/acceleration version of PF_CH_LECorrector_FixedP.
% It preserves exactly the same CH-LE correction equation, but solves
% dmu only in a dilated interface band. Outside the band, dmu_corr = 0.
%
% IMPORTANT:
%   This first reduced-band implementation still assembles the full sparse
%   matrix A and then factors only A(active,active). It is intentionally
%   conservative: first verify it reproduces the full corrector. If it
%   works, assembly can be restricted to the same band in a later cleanup.
%
% PURPOSE
%   Correct the CH/LE chemical field after one coupled AC-CH predictor and
%   one nonlinear LE update, without solving Allen-Cahn again.
%
% USAGE
%   STATE_RAW = PF_Coupled_ACCH_LETangent(...);
%   STATE_IT  = LE_Run(STATE_RAW,PARAM,MODEL);
%   [STATE_CORR,DIAG] = PF_CH_LECorrector_FixedP( ...
%                           STATE_OLD,STATE_IT,PARAM,GRID,PHYS,NUM);
%   STATE_NEW = LE_Run(STATE_CORR,PARAM,MODEL);
%
% INPUT STATE MEANING
%   STATE_TIME : fixed old-time state. Only STATE_TIME.E supplies E^n/dt.
%   STATE_IT   : current nonlinear iterate after the predictor LE_Run.
%                Its phi and p are frozen in this corrector.
%                Its E, mu_e and chi define the CH residual and tangent.
%
% EQUATION
%   For fixed updated p, solve one chemical correction dmu from:
%
%       [ D + A_E * chi_IT ] dmu
%           = E_TIME/dt - A_E*E_IT - D*mu_IT
%
%   then recover:
%
%       mu_CORR = mu_IT + dmu
%       E_CORR  = E_IT  + chi_IT*dmu
%
%   There is no e*dp term because phi and p are held fixed during this
%   corrector. LE_Run must be called afterwards to refresh c, e, chi and
%   omega at corrected E and fixed p.
%
% IMPORTANT
%   This is not a second physical time advance. It corrects the chemical
%   residual of the same timestep using STATE_TIME as the fixed old-time
%   reference. It preserves STATE_IT.phi and STATE_IT.p exactly.
%
% The spatial stencils and reflective boundary convention are copied from
% PF_Coupled_ACCH_LETangent in PFM_Version2.

% ------------------------------------------------------------
% Basic checks and sizes
% ------------------------------------------------------------
if nargin < 7
    error('PF_CH_LECorrector_FixedP_Band_CS requires STATE_TIME, STATE_IT, PARAM, MODEL, GRID, PHYS and NUM.')
end

% Backward compatibility with the temporary signature:
%   (STATE_TIME,STATE_IT,MODEL,PARAM,GRID,PHYS,NUM)
if isfield(PARAM,'pars') && isfield(MODEL,'M')
    tmp   = PARAM;
    PARAM = MODEL;
    MODEL = tmp;
end

E_time = STATE_TIME.E;
E_it   = STATE_IT.E;
mu_it  = STATE_IT.mu_e;
chi_it = STATE_IT.chi;

% Convex-split implicit capacity.
% Raw mu_e is still used as the physical driving force.
if isfield(PARAM,'use_CS_chi') && PARAM.use_CS_chi == 1
    chi_imp = ConvexifyChi_ForImplicit(chi_it,PARAM);
else
    chi_imp = chi_it;
end

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
% Spatially varying kappa: identical choice to coupled solver
% ------------------------------------------------------------
if isfield(PARAM,'kappa_eff') && ~isempty(PARAM.kappa_eff)

    kappa_eff = PARAM.kappa_eff;

    if isscalar(kappa_eff)
        kappa_eff = kappa_eff*ones(ny,nx);
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
% Optional composition-space kappa.
%
% If enabled, the fourth-order regularization is applied to phase
% composition c through dc = inv(H_c)*J'*dmu, instead of to total E.
% This keeps the existing dmu-corrector structure.
% ------------------------------------------------------------
use_kappa_c = isfield(PARAM,'use_kappa_c') && PARAM.use_kappa_c == 1;

if use_kappa_c
    KAPC = Build_KappaC_Data(STATE_IT,PARAM,MODEL,kappa_eff);
else
    KAPC = [];
end

% ------------------------------------------------------------
% Reflective-neighbour indices: identical stencil to coupled solver
% ------------------------------------------------------------
[Igrid,Jgrid] = ndgrid(1:ny,1:nx);

ii       = Igrid(:);
jj       = Jgrid(:);

refI     = @(i,sh) reflect_index(i+sh,ny);
refJ     = @(j,sh) reflect_index(j+sh,nx);

jjL      = refJ(jj,-1);
jjR      = refJ(jj,+1);
iiU      = refI(ii,-1);
iiD      = refI(ii,+1);

jjL2     = refJ(jj,-2);
jjR2     = refJ(jj,+2);
iiU2     = refI(ii,-2);
iiD2     = refI(ii,+2);

iiUR     = refI(ii,-1);
jjUR     = refJ(jj,+1);
iiDR     = refI(ii,+1);
jjDR     = refJ(jj,+1);
iiUL     = refI(ii,-1);
jjUL     = refJ(jj,-1);
iiDL     = refI(ii,+1);
jjDL     = refJ(jj,-1);

idx_c    = sub2ind([ny,nx],ii,jj);
idx_L    = sub2ind([ny,nx],ii,jjL);
idx_R    = sub2ind([ny,nx],ii,jjR);
idx_U    = sub2ind([ny,nx],iiU,jj);
idx_D    = sub2ind([ny,nx],iiD,jj);

idx_L2   = sub2ind([ny,nx],ii,jjL2);
idx_R2   = sub2ind([ny,nx],ii,jjR2);
idx_U2   = sub2ind([ny,nx],iiU2,jj);
idx_D2   = sub2ind([ny,nx],iiD2,jj);

idx_UR   = sub2ind([ny,nx],iiUR,jjUR);
idx_DR   = sub2ind([ny,nx],iiDR,jjDR);
idx_UL   = sub2ind([ny,nx],iiUL,jjUL);
idx_DL   = sub2ind([ny,nx],iiDL,jjDL);

% ------------------------------------------------------------
% Chemical-potential unknown numbering
% ------------------------------------------------------------
Nmu      = Ne*Nnode;
idMu     = cell(1,Ne);

for ie = 1:Ne
    idMu{ie} = ((ie-1)*Nnode+(1:Nnode)).';
end

% ------------------------------------------------------------
% Sparse allocation and residual
% Standard mode has 13 A_E*chi stencil entries. In composition-kappa mode
% A_E = 1/dt, so only the center chi entry is nonzero; the fourth-order
% c-kappa block is assembled separately.
% ------------------------------------------------------------
if use_kappa_c
    max_nnz  = Nnode*Ne*(5+Ne)+1000;
else
    max_nnz  = Nnode*Ne*(5+13*Ne)+1000;
end

% Extra allocation for aggregated composition-space kappa block.
% The old allocation did not include these entries, causing repeated growth
% of rows/cols/vals when use_kappa_c = 1.
if use_kappa_c && ~isempty(KAPC)
    max_nnz = max_nnz + Nnode * Ne * 13 * Ne;
end

rows     = zeros(max_nnz,1);
cols     = zeros(max_nnz,1);
vals     = zeros(max_nnz,1);
R        = zeros(Nmu,1);
k        = 1;

% ============================================================
% CH-only Newton correction block at fixed updated p
% ============================================================
for l = 1:Ne

    row  = idMu{l}(idx_c);

    % --------------------------------------------------------
    % Diffusion operator D on mu_l: copied from coupled solver
    % --------------------------------------------------------
    Ml   = PARAM.M{l};

    M_c  = Ml(idx_c);
    M_L  = Ml(idx_L);
    M_R  = Ml(idx_R);
    M_U  = Ml(idx_U);
    M_D  = Ml(idx_D);

    d_L  = -(M_L+M_c)/2/dx2;
    d_R  = -(M_R+M_c)/2/dx2;
    d_U  = -(M_U+M_c)/2/dy2;
    d_D  = -(M_D+M_c)/2/dy2;
    d_C  = -(d_L+d_R+d_U+d_D);

    % --------------------------------------------------------
    % A_E operator on E_l.
    %
    % Standard mode:
    %   A_E = 1/dt + lap(M*kappa_eff*lap(E))
    %
    % Composition-kappa mode:
    %   A_E = 1/dt only.
    %   The fourth-order term is added separately as a c-space
    %   operator: lap(M*kappa_eff*p_phase*lap(dc_phase)).
    % --------------------------------------------------------
    if use_kappa_c

        z     = zeros(size(idx_c));

        a_C   = 1/dt + z;
        a_L   = z;
        a_R   = z;
        a_U   = z;
        a_D   = z;

        a_L2  = z;
        a_R2  = z;
        a_U2  = z;
        a_D2  = z;

        a_UR  = z;
        a_DR  = z;
        a_UL  = z;
        a_DL  = z;

    else

        [q_C,q_L,q_R,q_U,q_D,q_L2,q_R2,q_U2,q_D2,q_UR,q_DR,q_UL,q_DL] = ...
            Kappa4_Coeffs(Ml.*kappa_eff,idx_c,idx_L,idx_R,idx_U,idx_D, ...
                          idx_L2,idx_R2,idx_U2,idx_D2,idx_UR,idx_DR,idx_UL,idx_DL, ...
                          dx2,dy2,dx4,dy4);

        a_C   = 1/dt + q_C;
        a_L   = q_L;
        a_R   = q_R;
        a_U   = q_U;
        a_D   = q_D;

        a_L2  = q_L2;
        a_R2  = q_R2;
        a_U2  = q_U2;
        a_D2  = q_D2;

        a_UR  = q_UR;
        a_DR  = q_DR;
        a_UL  = q_UL;
        a_DL  = q_DL;

    end

    % --------------------------------------------------------
    % Current nonlinear iterate residual
    %
    % The time RHS must use fixed E^n from STATE_TIME.
    % In composition-kappa mode, A_E = 1/dt and all neighbour A_E
    % coefficients are zero. Avoid gathering 13 E stencils there.
    % --------------------------------------------------------
    E_time_l = E_time{l};
    E_l      = E_it{l};
    mu_l     = mu_it{l};

    Etime_c  = E_time_l(idx_c);
    E_c      = E_l(idx_c);

    mu_c     = mu_l(idx_c);
    mu_L     = mu_l(idx_L);
    mu_R     = mu_l(idx_R);
    mu_U     = mu_l(idx_U);
    mu_D     = mu_l(idx_D);

    D_muit   = ...
        d_C.*mu_c + d_L.*mu_L + d_R.*mu_R + d_U.*mu_U + d_D.*mu_D;

    if use_kappa_c

        R(row) = (Etime_c - E_c)/dt - D_muit;

    else

        E_L      = E_l(idx_L);
        E_R      = E_l(idx_R);
        E_U      = E_l(idx_U);
        E_D      = E_l(idx_D);
        E_L2     = E_l(idx_L2);
        E_R2     = E_l(idx_R2);
        E_U2     = E_l(idx_U2);
        E_D2     = E_l(idx_D2);
        E_UR     = E_l(idx_UR);
        E_DR     = E_l(idx_DR);
        E_UL     = E_l(idx_UL);
        E_DL     = E_l(idx_DL);

        AE_Eit   = ...
            a_C .*E_c  + ...
            a_L .*E_L  + a_R .*E_R  + a_U .*E_U  + a_D .*E_D  + ...
            a_L2.*E_L2 + a_R2.*E_R2 + a_U2.*E_U2 + a_D2.*E_D2 + ...
            a_UR.*E_UR + a_DR.*E_DR + a_UL.*E_UL + a_DL.*E_DL;

        R(row) = Etime_c/dt - AE_Eit - D_muit;

    end

    % Composition-space fourth-order term.
    % This adds the implicit operator from
    %   dc_phase = inv(H_c)*J'*dmu
    % and the explicit current c contribution, while p is kept fixed.
    if use_kappa_c && ~isempty(KAPC)
        [R,rows,cols,vals,k] = Add_KappaC_Block( ...
            R,rows,cols,vals,k,row,idMu,l,Ml,kappa_eff,KAPC, ...
            idx_c,idx_L,idx_R,idx_U,idx_D,idx_L2,idx_R2,idx_U2,idx_D2, ...
            idx_UR,idx_DR,idx_UL,idx_DL,dx2,dy2,dx4,dy4);
    end

    % --------------------------------------------------------
    % D * dmu_l
    % --------------------------------------------------------
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_c),d_C);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_L),d_L);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_R),d_R);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_U),d_U);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{l}(idx_D),d_D);

    % --------------------------------------------------------
    % A_E * chi_IT * dmu
    % --------------------------------------------------------
    if use_kappa_c

        % Fast path: A_E = 1/dt in composition-kappa mode, so only the
        % center chi block is nonzero. The c-gradient fourth-order part is
        % already included by Add_KappaC_Block above.
        for m = 1:Ne

            Chi_lm  = chi_imp{l,m};
            Chi_c   = Chi_lm(idx_c);

            [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{m}(idx_c),Chi_c/dt);

        end

    else

        for m = 1:Ne

            Chi_lm  = chi_imp{l,m};

            Chi_c   = Chi_lm(idx_c);
            Chi_L   = Chi_lm(idx_L);
            Chi_R   = Chi_lm(idx_R);
            Chi_U   = Chi_lm(idx_U);
            Chi_D   = Chi_lm(idx_D);

            Chi_L2  = Chi_lm(idx_L2);
            Chi_R2  = Chi_lm(idx_R2);
            Chi_U2  = Chi_lm(idx_U2);
            Chi_D2  = Chi_lm(idx_D2);

            Chi_UR  = Chi_lm(idx_UR);
            Chi_DR  = Chi_lm(idx_DR);
            Chi_UL  = Chi_lm(idx_UL);
            Chi_DL  = Chi_lm(idx_DL);

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


end

% ------------------------------------------------------------
% Assemble full operator once, then solve only in an interface band.
%
% The fixed-p correction is observed to be localized near interfaces.
% We impose dmu_corr = 0 outside a dilated interface band. The band must
% be wide enough that the correction is already negligible at its edge.
% ------------------------------------------------------------
rows      = rows(1:k-1);
cols      = cols(1:k-1);
vals      = vals(1:k-1);

A         = sparse(rows,cols,vals,Nmu,Nmu);

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

if band_thick < 2
    error('NUM.CHLE_band_thick must be at least 2 because the CH operator contains a +/-2 stencil.')
end

%------------------------------------------------------------
% Correction-band core
%
% IMPORTANT:
% A phase interface is not always represented by cells having 0 < p < 1.
% In a sharp or under-resolved interface, two adjacent cells may be pure
% but their p fields jump from 0 to 1. Therefore the band core must include
% both mixed cells and neighbour jumps in p.
%------------------------------------------------------------
if (isfield(NUM,'CHLE_force_full') && NUM.CHLE_force_full == 1) || use_kappa_c

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

    %Optional residual-based safety extension. Use this if the chemical
    %correction has relevant support away from visible phase interfaces.
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

NUM_L = NUM;
if ~isfield(NUM_L,'linear_cache_id') || isempty(NUM_L.linear_cache_id)
    NUM_L.linear_cache_id = 'CHLE';
end
if isfield(NUM,'ilu_reuse_steps_CHLE') && ~isempty(NUM.ilu_reuse_steps_CHLE)
    NUM_L.ilu_reuse_steps = NUM.ilu_reuse_steps_CHLE;
end
cache_id = NUM_L.linear_cache_id;

if ~any(mask_corr(:))

    %No interface correction required.
    sol             = zeros(size(R));
    relres          = 0;
    active          = zeros(0,1);
    Ared            = sparse(0,0);
    Rred            = zeros(0,1);
    outside_res_ratio = norm(R)/max(norm(R),eps);
    LDIAG = struct();
    LDIAG.linear_solver = 'none';
    LDIAG.flag          = 0;
    LDIAG.relres        = 0;
    LDIAG.iter          = [0 0];
    LDIAG.solve_time    = 0;
    LDIAG.prec_time     = 0;
    LDIAG.fallback_used = 0;
    LDIAG.perm_method   = 'none';
    LDIAG.linear_cache_id = cache_id;
    LDIAG.ilu_cache_hit   = 0;
    LDIAG.ilu_cache_age   = -1;
    LDIAG.ilu_rebuilt     = 0;
else

    active_node = mask_corr(:);
    active      = zeros(Ne*nnz(active_node),1);
    kact        = 1;

    for ie = 1:Ne

        ids = idMu{ie}(active_node);
        active(kact:kact+numel(ids)-1) = ids;
        kact = kact+numel(ids);

    end

    Ared             = A(active,active);
    Rred             = R(active);
    [sol_red,LDIAG]  = PF_LinearSolve(Ared,Rred,NUM_L);
    sol              = zeros(size(R));
    sol(active)      = sol_red;
    relres           = LDIAG.relres;

    outside   = true(Nmu,1);
    outside(active) = false;
    outside_res_ratio = norm(R(outside))/max(norm(R),eps);

end

% ------------------------------------------------------------
% Unpack chemical-potential correction
% ------------------------------------------------------------
dmu       = cell(1,Ne);

for ie = 1:Ne
    dmu{ie} = reshape(sol(idMu{ie}),ny,nx);
end

% ------------------------------------------------------------
% Fixed-p corrected state
% ------------------------------------------------------------
STATE_CORR      = STATE_IT;

%Keep AC geometry exactly fixed during the chemical corrector.
STATE_CORR.phi = STATE_IT.phi;
STATE_CORR.p   = STATE_IT.p;

%Update mu around the current nonlinear iterate.
for ie = 1:Ne
    STATE_CORR.mu_e{ie} = STATE_IT.mu_e{ie}+dmu{ie};
end

%Recover E from the current LE tangent at fixed p.
for ie = 1:Ne

    En = STATE_IT.E{ie};

    for je = 1:Ne
        En = En+chi_imp{ie,je}.*dmu{je};
    end

    STATE_CORR.E{ie} = En;

end

if isfield(NUM,'norm_E') && NUM.norm_E == 1
    STATE_CORR.E = EnforceMeanE_Local(STATE_CORR.E,STATE_TIME.E);
end

%Predict omega only for diagnostic comparison before the next LE_Run.
%LE_Run should overwrite omega, e, chi and c after this function returns.
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

%Retain coefficient fields for diagnostics only; final LE_Run overwrites them.
STATE_CORR.chi = STATE_IT.chi;
STATE_CORR.e   = STATE_IT.e;
STATE_CORR.c   = STATE_IT.c;
if isfield(STATE_IT,'c_ext')
    STATE_CORR.c_ext = STATE_IT.c_ext;
end

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
DIAG.relres           = relres;
DIAG.linear_cache_id  = LDIAG.linear_cache_id;
DIAG.ilu_cache_hit    = LDIAG.ilu_cache_hit;
DIAG.ilu_cache_age    = LDIAG.ilu_cache_age;
DIAG.ilu_rebuilt      = LDIAG.ilu_rebuilt;
DIAG.linear_solver    = LDIAG.linear_solver;
DIAG.solve_flag       = LDIAG.flag;
DIAG.solve_relres     = LDIAG.relres;
DIAG.solve_iter       = LDIAG.iter;
DIAG.solve_time       = LDIAG.solve_time;
DIAG.prec_time        = LDIAG.prec_time;
DIAG.fallback_used    = LDIAG.fallback_used;
DIAG.perm_method      = LDIAG.perm_method;
DIAG.matrix_size_full = Nmu;
DIAG.nnz_full         = nnz(A);
DIAG.matrix_size      = numel(active);
DIAG.nnz              = nnz(Ared);
DIAG.solve_ratio      = numel(active)/max(Nmu,1);
DIAG.outside_res_ratio= outside_res_ratio;
DIAG.mask_core        = mask_core;
DIAG.mask_corr        = mask_corr;
DIAG.Ne               = Ne;
DIAG.band_thick       = band_thick;
DIAG.full_coverage    = all(mask_corr(:));
DIAG.n_core_nodes     = nnz(mask_core);
DIAG.n_corr_nodes     = nnz(mask_corr);

max_dmu          = 0;
max_dE           = 0;

for ie = 1:Ne

    max_dmu = max(max_dmu,max(abs(dmu{ie}(:))));
    max_dE  = max(max_dE,max(abs(STATE_CORR.E{ie}(:)-STATE_IT.E{ie}(:))));

end

DIAG.max_dmu     = max_dmu;
DIAG.max_dE      = max_dE;
DIAG.fixed_phi_error = max(abs(STATE_CORR.phi(:)-STATE_IT.phi(:)));
DIAG.fixed_p_error   = max(abs(STATE_CORR.p(:)-STATE_IT.p(:)));

%Check whether the imposed zero correction boundary is sufficiently far
%away from the interface core. A large value means the band is too narrow.
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

DIAG.edge_dmu_ratio = edge_dmu/max(max_dmu,eps);
DIAG.use_CS_chi     = isfield(PARAM,'use_CS_chi') && PARAM.use_CS_chi == 1;
DIAG.use_kappa_c    = use_kappa_c;
DIAG.fast_kappa_c_ch = use_kappa_c;
DIAG.use_c_ext_kappa = isfield(STATE_IT,'c_ext') && ~isempty(STATE_IT.c_ext);
DIAG.N_kappa_c     = numel(KAPC);
if isempty(KAPC)
    DIAG.kappa_c_grains = [];
else
    DIAG.kappa_c_grains = [KAPC.ig];
end

end


%% ========================================================================
% Local helper functions
% ========================================================================



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

% Use rim/ghost-extended composition for kappa(c) if available.
% Physical STATE.c is still used by LE and thermodynamics outside this block.
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
