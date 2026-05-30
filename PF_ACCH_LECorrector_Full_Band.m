function [STATE_CORR,DIAG] = PF_ACCH_LECorrector_Full_Band(STATE_TIME,STATE_IT,PARAM,MODEL,GRID,PHYS,NUM)
%PF_ACCH_LECORRECTOR_FULL_BAND
%
% Same-timestep full AC-CH-LE tangent corrector around a nonlinear LE state.
%
% INTENDED CALLING SEQUENCE
%
%   [STATE_RAW,DIAG0] = PF_Coupled_ACCH_LETangent( ...
%                            STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);
%
%   STATE_LE0          = LE_Run(STATE_RAW,PARAM,MODEL);
%
%   [STATE_CORR,DIAGC] = PF_ACCH_LECorrector_Full_Band( ...
%                            STATE_OLD,STATE_LE0,PARAM,MODEL,GRID,PHYS,NUM);
%
%   STATE_TRIAL        = LE_Run(STATE_CORR,PARAM,MODEL);
%
% MEANING OF THE TWO INPUT STATES
%
%   STATE_TIME : fixed accepted old-time state. It supplies phi^n and E^n
%                in the time-discrete residual and is never advanced here.
%
%   STATE_IT   : current nonlinear iterate for the same n -> n+1 step.
%                Normally STATE_IT = STATE_LE0, i.e. the LE-evaluated
%                result of the first coupled predictor. It supplies the
%                current phi, p, E, mu_e, chi, e and omega used to
%                linearize the residual.
%
% UNKNOWN VECTOR
%
%   x = [active correction dphi ; chemical correction dmu]
%
% Unlike PF_CH_LECorrector_FixedP, this function corrects both phase
% geometry and chemical redistribution. It is not a second physical
% timestep: the residual always uses STATE_TIME as the old-time reference.
%
% DISCRETE RESIDUALS
%
%   AC residual at STATE_IT:
%
%       R_AC = -(phi_IT-phi_TIME)/dt + LK*Lap(phi_IT) + S_AC(STATE_IT)
%
%   CH residual at STATE_IT:
%
%       R_CH = E_TIME/dt - A_E*E_IT - D*mu_IT
%
%   LE tangent used in the CH correction:
%
%       dE = chi_IT*dmu + B_IT*dphi
%
% with B_IT = dE/dphi at constant mu, evaluated from STATE_IT.
%
% BAND OPTION
%
%   dphi is always solved only around interfaces, as in the original
%   coupled predictor. By default dmu remains global, which is the safest
%   validation mode.
%
%   To restrict dmu to a dilated interface/residual band:
%
%       NUM.FullCorr_mu_band_enable = 1;
%       NUM.FullCorr_mu_band_thick  = 20;
%       NUM.FullCorr_res_rel        = 1e-8;   %or [] for geometry only
%
%   In band mode the code imposes dmu = 0 outside the selected band and
%   reports DIAG.relres_full. Do not use the band result unless it agrees
%   with the full-dmu corrector for representative timesteps.
%
% KAPPA
%
%   This function intentionally uses the same A_E fourth-order stencil as
%   the current PF_Coupled_ACCH_LETangent so the full corrector solves the
%   same discrete equation as your predictor. Since you currently set
%   PARAM.kappa_eff(:) = PHYS.kappa, this is consistent with a constant
%   kappa test. If you later restore spatially varying kappa, patch both
%   predictor and corrector together.

%% ------------------------------------------------------------------------
%  Inputs and sizes
% -------------------------------------------------------------------------
phi_time = STATE_TIME.phi;
E_time   = STATE_TIME.E;

phi_it   = STATE_IT.phi;
p_it     = STATE_IT.p;
E_it     = STATE_IT.E;
mu_it    = STATE_IT.mu_e;
chi0     = STATE_IT.chi;
e0       = STATE_IT.e;

[ny,nx,Ngrain] = size(phi_it);
Ne             = numel(E_time);
Nnode          = nx*ny;

dt             = NUM.dt_phy;
dx             = GRID.dx;
dy             = GRID.dy;

dx2            = dx^2;
dy2            = dy^2;
dx4            = dx2^2;
dy4            = dy2^2;

if dt <= 0 || dx <= 0 || dy <= 0
    error('NUM.dt_phy, GRID.dx and GRID.dy must be positive.')
end

%% ------------------------------------------------------------------------
%  Kappa field: same stencil convention as current predictor
% -------------------------------------------------------------------------
if isfield(PARAM,'kappa_eff') && ~isempty(PARAM.kappa_eff)

    kappa_eff = PARAM.kappa_eff;

    if isscalar(kappa_eff)
        kappa_eff = kappa_eff*ones(ny,nx);
    end

else

    kappa_eff = PHYS.kappa*ones(ny,nx);

end

%% ------------------------------------------------------------------------
%  AC source evaluated at current nonlinear LE iterate
% -------------------------------------------------------------------------
STATE_SRC = STATE_IT;

if isfield(NUM,'use_Aac') && NUM.use_Aac == 1

    if isfield(NUM,'Aac_fac')
        Aac_fac = NUM.Aac_fac;
    else
        Aac_fac = 3;
    end

    PARAM.A_ac = Calc_Aac_FrozenOmega( ...
        STATE_SRC,PARAM,MODEL,Aac_fac,1e-6,0,[]);

else

    PARAM.A_ac = zeros(ny,nx);

end

STATE_SRC = Calc_S_AllenCahn(STATE_SRC,PARAM,MODEL);
S_AC      = STATE_SRC.S_AC;

%% ------------------------------------------------------------------------
%  Grid indexing with reflective boundary
% -------------------------------------------------------------------------
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

iiUR     = refI(ii,-1); jjUR = refJ(jj,+1);
iiDR     = refI(ii,+1); jjDR = refJ(jj,+1);
iiUL     = refI(ii,-1); jjUL = refJ(jj,-1);
iiDL     = refI(ii,+1); jjDL = refJ(jj,-1);

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

%% ------------------------------------------------------------------------
%  Active phi-correction mask, evaluated at STATE_IT
% -------------------------------------------------------------------------
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
    if PHYS.kappa ~= 0
        mask_thick = 2;
    else
        mask_thick = 1;
    end
end

maskPhi = Local_Calc_Interface_Mask(phi_it,phi_cut,pure_cut,mask_thick);

if isfield(NUM,'phi_mask_source_tol') && ~isempty(NUM.phi_mask_source_tol)

    for alpha = 1:Ngrain

        core_source        = abs(S_AC{alpha}) > NUM.phi_mask_source_tol;
        maskPhi(:,:,alpha) = maskPhi(:,:,alpha) | ...
                              Local_Dilate_Mask(core_source,mask_thick);

    end

end

maskPhi        = logical(maskPhi);
active_per_cell= sum(maskPhi,3);

%% ------------------------------------------------------------------------
%  Optional diagonal AC-source Jacobian around STATE_IT
% -------------------------------------------------------------------------
if isfield(NUM,'use_Jphi') && NUM.use_Jphi == 1
    Jphi_diag = Calc_Jphi_Diag_FD(STATE_SRC,PARAM,MODEL,NUM,maskPhi);
else
    Jphi_diag = [];
end

%% ------------------------------------------------------------------------
%  Unknown ids: interface dphi and initially global dmu
% -------------------------------------------------------------------------
idPhiMap = cell(1,Ngrain);
counter  = 0;

for alpha = 1:Ngrain

    idmap       = zeros(ny,nx);
    ids         = find(maskPhi(:,:,alpha));
    nids        = numel(ids);

    idmap(ids)  = counter+(1:nids);

    idPhiMap{alpha} = idmap;
    counter     = counter+nids;

end

Nphi = counter;
Nmu  = Ne*Nnode;
Ntot = Nphi+Nmu;

idMu = cell(1,Ne);

for ie = 1:Ne
    idMu{ie} = Nphi+((ie-1)*Nnode+(1:Nnode)).';
end

%% ------------------------------------------------------------------------
%  LE tangent ingredients evaluated at STATE_IT
% -------------------------------------------------------------------------
eps_phi = 1e-14;
Dphi    = sum(phi_it.^2,3)+eps_phi;
p_tan   = phi_it.^2./Dphi;

facPhi  = cell(1,Ngrain);

for alpha = 1:Ngrain
    facPhi{alpha} = 2.*phi_it(:,:,alpha)./Dphi;
end

e_bar = cell(1,Ne);

for ie = 1:Ne

    tmp = zeros(ny,nx);

    for ig = 1:Ngrain
        tmp = tmp+p_tan(:,:,ig).*e0{ig}{ie};
    end

    e_bar{ie} = tmp;

end

%% ------------------------------------------------------------------------
%  Allocate sparse full coupled correction system
% -------------------------------------------------------------------------
max_nnz = ...
    Nphi*(6+Ne) + ...
    Nnode*Ne*(5+13*Ne) + ...
    13*Ne*max(Nphi,1) + ...
    1000;

rows = zeros(max_nnz,1);
cols = zeros(max_nnz,1);
vals = zeros(max_nnz,1);
R    = zeros(Ntot,1);

k = 1;

%% ========================================================================
%  1. Allen-Cahn correction block
% ========================================================================
for alpha = 1:Ngrain

    mask_a = maskPhi(:,:,alpha);
    idx_a  = find(mask_a);

    if isempty(idx_a)
        continue
    end

    row          = idPhiMap{alpha}(idx_a);
    phi_a_it     = phi_it(:,:,alpha);
    phi_a_time   = phi_time(:,:,alpha);
    lap_phi_it   = laplacian_reflect(phi_a_it,dx,dy);

    %True residual at the current nonlinear iterate.
    rhs_full     = -(phi_a_it-phi_a_time)./dt ...
                   + PARAM.LK.*lap_phi_it + S_AC{alpha};

    R(row)       = rhs_full(idx_a);

    LKc          = PARAM.LK(idx_a);
    Aac          = PARAM.A_ac(idx_a);

    cC           = 1/dt+Aac+2*LKc/dx2+2*LKc/dy2;
    cL           = -LKc/dx2;
    cR           = -LKc/dx2;
    cU           = -LKc/dy2;
    cD           = -LKc/dy2;

    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k, ...
        row,idPhiMap{alpha},idx_a,      cC);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k, ...
        row,idPhiMap{alpha},idx_L(idx_a),cL);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k, ...
        row,idPhiMap{alpha},idx_R(idx_a),cR);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k, ...
        row,idPhiMap{alpha},idx_U(idx_a),cU);
    [rows,cols,vals,k] = add_active_block(rows,cols,vals,k, ...
        row,idPhiMap{alpha},idx_D(idx_a),cD);

    if ~isempty(Jphi_diag)

        Jaa = Jphi_diag{alpha};

        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k, ...
            row,idPhiMap{alpha},idx_a,-Jaa(idx_a));

    end

    %AC coupling to chemical correction dmu.
    for ie = 1:Ne

        B_ie_alpha = facPhi{alpha}.*(e0{alpha}{ie}-e_bar{ie});
        coeff_mu   = -PARAM.L.*B_ie_alpha;

        [rows,cols,vals,k] = add_block(rows,cols,vals,k, ...
            row,idMu{ie}(idx_a),coeff_mu(idx_a));

    end

end

%% ========================================================================
%  2. Cahn-Hilliard correction block
% ========================================================================
for ie = 1:Ne

    row  = idMu{ie}(idx_c);

    Ml   = PARAM.M{ie};
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

    %A_E operator: copied from the present coupled predictor so this
    %corrector targets the same current discrete equation.
    kappa_c = kappa_eff(idx_c);
    MK_c    = M_c.*kappa_c;

    q_L     = MK_c.*(-4/dx4-4/(dx2*dy2));
    q_R     = MK_c.*(-4/dx4-4/(dx2*dy2));
    q_U     = MK_c.*(-4/dy4-4/(dx2*dy2));
    q_D     = MK_c.*(-4/dy4-4/(dx2*dy2));

    q_L2    = MK_c.*(1/dx4);
    q_R2    = MK_c.*(1/dx4);
    q_U2    = MK_c.*(1/dy4);
    q_D2    = MK_c.*(1/dy4);

    q_UR    = MK_c.*(2/(dx2*dy2));
    q_DR    = MK_c.*(2/(dx2*dy2));
    q_UL    = MK_c.*(2/(dx2*dy2));
    q_DL    = MK_c.*(2/(dx2*dy2));

    q_C     = MK_c.*(6/dx4+6/dy4+8/(dx2*dy2));

    a_C     = 1/dt+q_C;
    a_L     = q_L;
    a_R     = q_R;
    a_U     = q_U;
    a_D     = q_D;
    a_L2    = q_L2;
    a_R2    = q_R2;
    a_U2    = q_U2;
    a_D2    = q_D2;
    a_UR    = q_UR;
    a_DR    = q_DR;
    a_UL    = q_UL;
    a_DL    = q_DL;

    Et      = E_time{ie};
    Ei      = E_it{ie};
    mui     = mu_it{ie};

    AE_Eit = ...
        a_C .*Ei(idx_c)  + ...
        a_L .*Ei(idx_L)  + a_R .*Ei(idx_R)  + ...
        a_U .*Ei(idx_U)  + a_D .*Ei(idx_D)  + ...
        a_L2.*Ei(idx_L2) + a_R2.*Ei(idx_R2) + ...
        a_U2.*Ei(idx_U2) + a_D2.*Ei(idx_D2) + ...
        a_UR.*Ei(idx_UR) + a_DR.*Ei(idx_DR) + ...
        a_UL.*Ei(idx_UL) + a_DL.*Ei(idx_DL);

    D_muit = ...
        d_C.*mui(idx_c) + d_L.*mui(idx_L) + d_R.*mui(idx_R) + ...
        d_U.*mui(idx_U) + d_D.*mui(idx_D);

    %True same-timestep CH residual at current nonlinear iterate.
    R(row) = Et(idx_c)/dt-AE_Eit-D_muit;

    %D*dmu diagonal-element block.
    [rows,cols,vals,k] = add_block(rows,cols,vals,k, ...
        row,idMu{ie}(idx_c),d_C);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k, ...
        row,idMu{ie}(idx_L),d_L);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k, ...
        row,idMu{ie}(idx_R),d_R);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k, ...
        row,idMu{ie}(idx_U),d_U);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k, ...
        row,idMu{ie}(idx_D),d_D);

    %A_E*chi*dmu.
    for je = 1:Ne

        C = chi0{ie,je};

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_c), a_C .*C(idx_c));
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_L), a_L .*C(idx_L));
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_R), a_R .*C(idx_R));
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_U), a_U .*C(idx_U));
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_D), a_D .*C(idx_D));

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_L2),a_L2.*C(idx_L2));
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_R2),a_R2.*C(idx_R2));
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_U2),a_U2.*C(idx_U2));
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_D2),a_D2.*C(idx_D2));

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_UR),a_UR.*C(idx_UR));
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_DR),a_DR.*C(idx_DR));
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_UL),a_UL.*C(idx_UL));
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idMu{je}(idx_DL),a_DL.*C(idx_DL));

    end

    %A_E*B*dphi.
    for alpha = 1:Ngrain

        idmap = idPhiMap{alpha};

        if ~any(idmap(:))
            continue
        end

        B = facPhi{alpha}.*(e0{alpha}{ie}-e_bar{ie});

        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_c, a_C .*B(idx_c));
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_L, a_L .*B(idx_L));
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_R, a_R .*B(idx_R));
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_U, a_U .*B(idx_U));
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_D, a_D .*B(idx_D));

        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_L2,a_L2.*B(idx_L2));
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_R2,a_R2.*B(idx_R2));
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_U2,a_U2.*B(idx_U2));
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_D2,a_D2.*B(idx_D2));

        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_UR,a_UR.*B(idx_UR));
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_DR,a_DR.*B(idx_DR));
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_UL,a_UL.*B(idx_UL));
        [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_DL,a_DL.*B(idx_DL));

    end

end

%% ------------------------------------------------------------------------
%  Assemble full correction matrix
% -------------------------------------------------------------------------
rows = rows(1:k-1);
cols = cols(1:k-1);
vals = vals(1:k-1);

A    = sparse(rows,cols,vals,Ntot,Ntot);

%% ------------------------------------------------------------------------
%  Optional dmu correction band
% -------------------------------------------------------------------------
if ~isfield(NUM,'FullCorr_mu_band_enable')
    NUM.FullCorr_mu_band_enable = 0;
end

if ~isfield(NUM,'FullCorr_mu_band_thick')
    NUM.FullCorr_mu_band_thick = 20;
end

if ~isfield(NUM,'FullCorr_res_rel')
    NUM.FullCorr_res_rel = [];
end

if NUM.FullCorr_mu_band_enable == 1

    maskMu = any(maskPhi,3);

    if ~isempty(NUM.FullCorr_res_rel) && NUM.FullCorr_res_rel > 0

        resMu = zeros(ny,nx);

        for ie = 1:Ne
            resMu = max(resMu,abs(reshape(R(idMu{ie}),ny,nx)));
        end

        rmax = max(resMu(:));

        if rmax > 0
            maskMu = maskMu | resMu > NUM.FullCorr_res_rel*rmax;
        end

    end

    maskMu = Local_Dilate_Mask(maskMu,NUM.FullCorr_mu_band_thick);

    idSolve = (1:Nphi).';

    for ie = 1:Ne
        idSolve = [idSolve; idMu{ie}(maskMu(:))]; %#ok<AGROW>
    end

else

    maskMu  = true(ny,nx);
    idSolve = (1:Ntot).';

end

Ared = A(idSolve,idSolve);
Rred = R(idSolve);

%% ------------------------------------------------------------------------
%  Solve correction
% -------------------------------------------------------------------------
q        = colamd(Ared);
y        = Ared(:,q)\Rred;

sol_red  = zeros(size(Rred));
sol_red(q) = y;

sol      = zeros(Ntot,1);
sol(idSolve) = sol_red;

DIAG.relres_reduced = norm(Ared*sol_red-Rred)/max(norm(Rred),eps);
DIAG.relres_full    = norm(A*sol-R)/max(norm(R),eps);

%% ------------------------------------------------------------------------
%  Unpack correction and update current iterate
% -------------------------------------------------------------------------
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

for ie = 1:Ne
    dmu{ie} = reshape(sol(idMu{ie}),ny,nx);
end

STATE_CORR      = STATE_IT;
STATE_CORR.phi = phi_it+dphi;

if isfield(NUM,'norm_phi') && NUM.norm_phi == 1
    STATE_CORR.phi = Normalize_Phi_L2_Local(STATE_CORR.phi);
end

if isfield(NUM,'cut_phi') && NUM.cut_phi == 1
    STATE_CORR.phi = Cut_Phi_Local(STATE_CORR.phi);
end

STATE_CORR.p = Calc_p(MODEL,STATE_CORR.phi);

for ie = 1:Ne
    STATE_CORR.mu_e{ie} = mu_it{ie}+dmu{ie};
end

for ie = 1:Ne

    En = E_it{ie};

    for je = 1:Ne
        En = En+chi0{ie,je}.*dmu{je};
    end

    for ig = 1:Ngrain
        dp = STATE_CORR.p(:,:,ig)-p_it(:,:,ig);
        En = En+e0{ig}{ie}.*dp;
    end

    STATE_CORR.E{ie} = En;

end

if isfield(NUM,'norm_E') && NUM.norm_E == 1
    STATE_CORR.E = EnforceMeanE_Local(STATE_CORR.E,E_time);
end

%Linear omega prediction only; final LE_Run must overwrite it.
STATE_CORR.omg = STATE_IT.omg;

for ig = 1:Ngrain

    domega = zeros(ny,nx);

    for ie = 1:Ne
        domega = domega-e0{ig}{ie}.*dmu{ie};
    end

    STATE_CORR.omg(:,:,ig) = STATE_IT.omg(:,:,ig)+domega;

end

omg_mean = mean(STATE_CORR.omg,3);

for ig = 1:Ngrain
    STATE_CORR.omg(:,:,ig) = STATE_CORR.omg(:,:,ig)-omg_mean;
end

STATE_CORR.chi = chi0;
STATE_CORR.e   = e0;
STATE_CORR.c   = STATE_IT.c;

%% ------------------------------------------------------------------------
%  Diagnostics and safety guard
% -------------------------------------------------------------------------
DIAG.matrix_size      = numel(idSolve);
DIAG.matrix_size_full = Ntot;
DIAG.solve_ratio      = numel(idSolve)/max(Ntot,1);
DIAG.nnz              = nnz(Ared);
DIAG.nnz_full         = nnz(A);
DIAG.Nphi_active      = Nphi;
DIAG.Nmu              = Nmu;
DIAG.maskPhi          = maskPhi;
DIAG.maskMu           = maskMu;
DIAG.active_cell_mean = mean(active_per_cell(:));

DIAG.max_dphi         = max(abs(dphi(:)));
DIAG.max_dmu          = 0;
DIAG.max_dE           = 0;
DIAG.max_E            = 0;
DIAG.max_dp           = max(abs(STATE_CORR.p(:)-p_it(:)));

for ie = 1:Ne

    DIAG.max_dmu = max(DIAG.max_dmu,max(abs(dmu{ie}(:))));
    DIAG.max_dE  = max(DIAG.max_dE,max(abs(STATE_CORR.E{ie}(:)-E_it{ie}(:))));
    DIAG.max_E   = max(DIAG.max_E,max(abs(STATE_CORR.E{ie}(:))));

end

if isfield(NUM,'FullCorr_E_guard') && ~isempty(NUM.FullCorr_E_guard) && ...
   DIAG.max_E > NUM.FullCorr_E_guard

    error(['PF_ACCH_LECorrector_Full_Band generated nonphysical E: ', ...
           'max|E| = %.3e exceeds guard %.3e. Reduce dt or reject step.'], ...
           DIAG.max_E,NUM.FullCorr_E_guard)

end

end


%% ========================================================================
%  Local helper functions
% ========================================================================

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
Jdiag  = cell(1,Ngrain);

for alpha = 1:Ngrain

    Jdiag{alpha} = zeros(ny,nx);
    mask_a       = maskPhi(:,:,alpha);

    if ~any(mask_a(:))
        continue
    end

    STATEP = STATE;
    tmp    = STATEP.phi(:,:,alpha);
    tmp(mask_a) = tmp(mask_a)+eps_fd;
    STATEP.phi(:,:,alpha) = tmp;

    STATEP = Calc_S_AllenCahn(STATEP,PARAM,MODEL);
    Jtmp   = (STATEP.S_AC{alpha}-S0{alpha})/eps_fd;

    Jdiag{alpha}(mask_a) = Jtmp(mask_a);

end

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

L = (AL-2*A+AR)/dx^2+(AU-2*A+AD)/dy^2;

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


function [rows,cols,vals,k] = add_active_block(rows,cols,vals,k,row,idmap,idx_col,coeff)

col   = idmap(idx_col);
row   = row(:);
col   = col(:);
coeff = coeff(:);

keep  = col > 0 & isfinite(coeff) & coeff ~= 0;

if ~any(keep)
    return
end

r = row(keep);
c = col(keep);
v = coeff(keep);
n = numel(r);

if k+n-1 > numel(rows)

    grow = max(numel(rows),n+1000);
    rows = [rows; zeros(grow,1)];
    cols = [cols; zeros(grow,1)];
    vals = [vals; zeros(grow,1)];

end

rows(k:k+n-1) = r;
cols(k:k+n-1) = c;
vals(k:k+n-1) = v;
k             = k+n;

end


function E = EnforceMeanE_Local(E,E_old)

for ie = 1:numel(E)

    target = mean(E_old{ie}(:));
    current= mean(E{ie}(:));

    E{ie} = E{ie}+target-current;

end

end


function phi = Cut_Phi_Local(phi)

phi(phi < 0) = 0;
phi(phi > 1) = 1;

end


function phi = Normalize_Phi_L2_Local(phi)

phi = max(phi,0);
phi = min(phi,1);

s    = sqrt(sum(phi.^2,3));
mask = s > eps;

for ip = 1:size(phi,3)

    tmp       = phi(:,:,ip);
    tmp(mask) = tmp(mask)./s(mask);
    phi(:,:,ip) = tmp;

end

end


function mask = Local_Calc_Interface_Mask(phi,low_cut,pure_cut,thickness)

[ny,nx,Ngrain] = size(phi);
mask            = false(ny,nx,Ngrain);
den             = sum(phi.^2,3)+eps;

for alpha = 1:Ngrain

    q    = phi(:,:,alpha).^2./den;
    core = q > low_cut & q < 1-pure_cut;

    if nx == 1
        qL = q; qR = q;
    else
        qL = q(:,[2,1:nx-1]);
        qR = q(:,[2:nx,nx-1]);
    end

    if ny == 1
        qU = q; qD = q;
    else
        qU = q([2,1:ny-1],:);
        qD = q([2:ny,ny-1],:);
    end

    jump = abs(q-qL) > low_cut | abs(q-qR) > low_cut | ...
           abs(q-qU) > low_cut | abs(q-qD) > low_cut;

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
