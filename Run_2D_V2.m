%Clear and restart
clear;figure(1);clf;addpath([cd,'\bin']);addpath([cd,'\Thermo']);addpath([cd,'\Thermo\Solutions'])

%Load map
load Map2d.mat

%PHYSICS
PHYS.E_sc          =  E_sc;
PHYS.t_sc          =  1;                                                   % Time scale
PHYS.L_sc          =  1;                                                   % Length scale
PHYS.l             =  Lx/150/L_sc;                                         % interface thickness (m)
PHYS.sigma         =  0.3/PHYS.E_sc*PHYS.L_sc^2;                           % surface energy (J/m^2)
PHYS.kappa         =  1e-7/(PHYS.E_sc*PHYS.L_sc^2);                        % 4th order term, can be set to 0 if no solvus
PHYS.D_esti        =  1e-12;
PHYS.chi_ref       =  1e-2;
PHYS.M0            =  PHYS.D_esti*PHYS.t_sc/PHYS.L_sc^2*PHYS.chi_ref;
PHYS.m             =  6*PHYS.sigma/PHYS.l;
PHYS.kap           =  3/4*PHYS.sigma*PHYS.l;
PHYS.dceq          =  0.2;
PHYS.L             =  4*PHYS.m/3/PHYS.kap/(PHYS.dceq^2/PHYS.M0)/2000;
PHYS.eta           =  eta;

%NUMERICS
NUM.dt_phy         =   1e-4/PHYS.t_sc;
NUM.dt_max         =      1/PHYS.t_sc;
NUM.dt_min         =  1e-16/PHYS.t_sc; 
NUM.t_tot          =    1e5/PHYS.t_sc;
NUM.dE_target      =  1e-2;
NUM.dp_target      =  1e-2;
NUM.dmu_target     =  1e5;
NUM.time           =  0;
NUM.dt_good_count  =  0;
NUM.dt_grow_after  =  8;
NUM.dt_grow_fac    =  1.15;
NUM.dt_shrink_fac  =  0.5;
NUM.err_grow       =  0.25;

%GRIDS
GRID.dx            =  dx;
GRID.dy            =  dy;
GRID.nx            =  nx;
GRID.ny            =  ny;

%MODELS
MODEL.dgdphi       =  @(phi) 2*PHYS.m*phi.*(phi - 1).^2 + PHYS.m*phi.^2.*(2.*phi - 2);
MODEL.pars         =  pars;
MODEL.p_fun        =  @(a,phi)   phi(:,:,a).^2./sum(phi.^2,3);
MODEL.dpdphi       =  @(a,b,phi) (a==b)*2*phi(:,:,b)./sum(phi.^2,3) - 2*phi(:,:,a).*phi(:,:,b).^2./sum(phi.^2,3).^2;

%PARAMETERS
PARAM.L            =  PHYS.L*ones(ny,nx);
PARAM.Lm           =  PHYS.L*PHYS.m.*ones(ny,nx);
PARAM.LK           =  PHYS.L*PHYS.kap.*ones(ny,nx);
PARAM.Np           =  length(c);
PARAM.Ne           =  length(E);
PARAM.M            =  repmat({PHYS.M0*ones(ny,nx)},1,PARAM.Ne);

%STATES
STATE.c            =  c;
STATE.e            =  e;
STATE.E            =  E;
STATE.mu_e         =  mu_e;
STATE.chi          =  chi;
STATE.omg          =  zeros(ny,nx,Np);
STATE.phi          =  phi;
STATE.p            =  Calc_p(MODEL,STATE.phi);
STATE.mask         =  ones(ny,nx,Np);
STATE.LE_state     = [   ];

%DISPLAY COMPOSITION

disp([mean(STATE.E{1},'all') mean(STATE.E{2},'all') mean(STATE.E{3},'all') mean(STATE.E{end},'all')])
% NUM.dE_target      =  3e-3;
% NUM.dp_target      =  3e-3;

load 1000
for it = it:1e5
    if mod(it,100)==0
        save(num2str(it))
    end

    tic
    
    %DIRECT COUPLED SOLVER
    %OLD STATE
    STATE_OLD            =    STATE;
    STATE_TRIAL          =    STATE;

    
    % METHOD1
    % Make sure old state is thermodynamically consistent
    PARAM.eta            =    Eta_Damping(STATE_OLD.p,PHYS.eta,0.1*PHYS.eta);
    STATE_OLD            =    LE_Run(STATE_OLD,PARAM,MODEL);

    % One monolithic tangent AC-CH-LE step
    [STATE_TRIAL,DIAG]   =    PF_Coupled_ACCH_LETangent(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);

    % Final nonlinear LE correction only once
    PARAM.eta            =    Eta_Damping(STATE_TRIAL.p,PHYS.eta,0.1*PHYS.eta);
    STATE_TRIAL          =    LE_Run(STATE_TRIAL,PARAM,MODEL);



    
    % % METHOD 2
    % %EXPLICIT AC CH
    % %===========================TRIAL STEPS================================
    % PARAM.A_ac           =    Calc_Aac_FrozenOmega(STATE,PARAM,MODEL,3,1e-6,0,[]);
    % 
    % %ALLEN CAHN
    % STATE_TRIAL          =    Calc_S_AllenCahn(STATE_TRIAL,PARAM,MODEL);
    % STATE_TRIAL          =    PF_IMEX_Solver2D_AllenCahn_Stab(STATE_TRIAL,MODEL,PARAM,GRID,NUM,1);
    % 
    % %INTERFACE PENALTY
    % PARAM.eta            =    Eta_Damping(STATE_TRIAL.p,PHYS.eta,0.5*PHYS.eta); 
    % 
    % %LOCAL EQUILIBRIUM
    % STATE_TRIAL          =    LE_Run(STATE_TRIAL,PARAM,MODEL);
    % 
    % %CAHN HILLARD
    % STATE_TRIAL          =    Calc_S_Diffusion(STATE_TRIAL,STATE_OLD,NUM);
    % STATE_TRIAL          =    PF_IMEX_Solver2D_Diffusion_MuOnly(STATE_TRIAL,PARAM,GRID,PHYS,NUM);


    % % METHOD 3
    % % ============================================================
    % % Fixed old state
    % % ============================================================
    % STATE_OLD = STATE;
    % 
    % % Freeze eta during the whole trial step.
    % % Do not update eta inside Picard loop.
    % PARAM_STEP     = PARAM;
    % PARAM_STEP.eta = Eta_Damping(STATE_OLD.p,PHYS.eta,0.1*PHYS.eta);
    % 
    % % Make old state thermodynamically consistent
    % STATE_OLD = LE_Run(STATE_OLD,PARAM_STEP,MODEL);
    % 
    % % ============================================================
    % % Picard iteration:
    % %   fixed time reference = STATE_OLD
    % %   updated coefficients = STATE_COEF
    % % ============================================================
    % maxPicard = 5;
    % 
    % tol_phi = 1e-4;
    % tol_p   = 1e-4;
    % tol_E   = 0.2 * NUM.dE_target;
    % tol_mu  = 0.2 * NUM.dmu_target;
    % 
    % % Coefficient-state damping.
    % % Set alpha_coef = 1 for no damping.
    % % Use 0.3~0.7 if LE_Run changes mu/omega sharply at interfaces.
    % alpha_coef = 0.5;
    % 
    % STATE_COEF  = STATE_OLD;
    % STATE_TRIAL = STATE_OLD;
    % 
    % for ipic = 1:maxPicard
    % 
    %     STATE_PREV = STATE_TRIAL;
    % 
    %     % Same time step from STATE_OLD.
    %     % STATE_COEF only provides chi/e/omega coefficients.
    %     [STATE_RAW,DIAG_COUPLED] = PF_Coupled_ACCH_LETangent( ...
    %         STATE_OLD,PARAM_STEP,MODEL,GRID,PHYS,NUM,STATE_COEF);
    % 
    %     % Nonlinear LE correction
    %     STATE_TRIAL = LE_Run(STATE_RAW,PARAM_STEP,MODEL);
    % 
    %     % Picard convergence check
    %     dphi_pic = max(abs(STATE_TRIAL.phi(:) - STATE_PREV.phi(:)));
    %     dp_pic   = max(abs(STATE_TRIAL.p(:)   - STATE_PREV.p(:)));
    %     dE_pic   = Max_Cell_Diff(STATE_TRIAL.E,    STATE_PREV.E);
    %     dmu_pic  = Max_Cell_Diff(STATE_TRIAL.mu_e, STATE_PREV.mu_e);
    % 
    %     fprintf('Picard %d: dphi=%e, dp=%e, dE=%e, dmu=%e, linres=%e\n', ...
    %         ipic,dphi_pic,dp_pic,dE_pic,dmu_pic,DIAG_COUPLED.relres);
    % 
    %     if ipic > 1 && ...
    %             dphi_pic < tol_phi && ...
    %             dp_pic   < tol_p   && ...
    %             dE_pic   < tol_E   && ...
    %             dmu_pic  < tol_mu
    %         break
    %     end
    % 
    %     % Update only thermodynamic coefficients for next Picard iteration.
    %     % The time reference remains STATE_OLD.
    %     if alpha_coef >= 1
    %         STATE_COEF = STATE_TRIAL;
    %     else
    %         STATE_COEF = Blend_Tangent_State(STATE_COEF,STATE_TRIAL,alpha_coef);
    %     end
    % end






    %TIME STEP UPDATE
    [STATE,NUM]          =    Update_TimeStep_Soft(STATE,STATE_TRIAL,PARAM,MODEL,NUM);

    toc
    
    %Plotting
    TIME(it)             =    NUM.time;
    DTPHY(it)            =    NUM.dt_phy;
    PHASE1(it)           =    mean(mean(STATE.p(:,:,1)));
    PHASE2(it)           =    mean(mean(STATE.p(:,:,2)));

    if mod(it,5)==0
        disp(NUM.dt_phy)
        disp([mean(mean(STATE.p(:,:,1))),mean(mean(STATE.p(:,:,2))),mean(STATE.E{3}(:)),mean(STATE.E{4}(:))])
        subplot(331);plot(GRID.x,STATE.E{1}(3,:),GRID.x,STATE.E{2}(3,:),GRID.x,STATE.E{3}(3,:),GRID.x,STATE.E{4}(3,:));title('E1')
        subplot(332);plot(STATE.mu_e{1}(3,:));title('mu_e')
        subplot(333);plot(DTPHY,'b.');title('dt')
        subplot(334);plot(GRID.x,STATE.phi(3,:,1),'.-',GRID.x,STATE.phi(3,:,2),'.-',GRID.x,STATE.phi(3,:,3),'.-',GRID.x,STATE.phi(3,:,end),'.-');title('p2')        
        subplot(335);plot(GRID.x,STATE.omg(3,:,1)-STATE.omg(3,:,2),GRID.x,STATE.omg(3,:,2)-STATE.omg(3,:,3),'.-');title('domg12')
        subplot(336);plot(PARAM.eta(3,:));title('eta')
        subplot(337);plot(GRID.x,STATE.c{1}{1}(3,:),GRID.x,STATE.c{1}{2}(3,:),GRID.x,STATE.c{1}{3}(3,:),GRID.x,STATE.c{1}{4}(3,:));title('c11')
        subplot(338);plot(GRID.x,STATE.c{2}{1}(3,:),GRID.x,STATE.c{2}{2}(3,:));title('c21')
        subplot(339);plot(TIME,PHASE2,'.-');title('Phase2')
        drawnow
    end
end







%% FUNCTIONS
function [STATE_NEW,DIAG] = PF_Coupled_ACCH_LETangent(STATE_REF,PARAM,MODEL,GRID,PHYS,NUM,STATE_COEF)
%PF_COUPLED_ACCH_LETANGENT
%
% Coupled AC-CH tangent solve with optional Picard coefficient state.
%
% STATE_REF:
%   fixed old accepted state at time t^n.
%   The unknown increments are:
%       dphi = phi_new - phi_ref
%       dmu  = mu_new  - mu_ref
%
% STATE_COEF:
%   optional thermodynamic coefficient state.
%   It provides:
%       chi, e, omega, mu_e for AC driving/tangent coefficients.
%
% If STATE_COEF is omitted, the function behaves like the original
% one-step coupled tangent solver.
%
% The solved LE tangent is:
%
%   dE = chi_coef*dmu + sum_i e_i_coef*dp_i
%
% with:
%
%   dp_i = sum_alpha dp_i/dphi_alpha at STATE_REF * dphi_alpha
%
% This function does not call LE_Run. Call LE_Run outside after this solve.

if nargin < 7 || isempty(STATE_COEF)
    STATE_COEF = STATE_REF;
end

% ------------------------------------------------------------
% Fixed time-reference state: increments are measured from here
% ------------------------------------------------------------
phi_ref = STATE_REF.phi;
p_ref   = STATE_REF.p;
E_ref   = STATE_REF.E;
mu_ref  = STATE_REF.mu_e;

% ------------------------------------------------------------
% Coefficient state: thermodynamic coefficients only
% ------------------------------------------------------------
chi0 = STATE_COEF.chi;
e0   = STATE_COEF.e;

% Sizes
[ny,nx,Np] = size(phi_ref);
Ne         = numel(E_ref);
Nnode      = nx*ny;

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

PARAM.A_ac = Calc_Aac_FrozenOmega(STATE_SRC,PARAM,MODEL,3,1e-6,0,[]);
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

% ------------------------------------------------------------
% Unknown ids
% ------------------------------------------------------------
Nphi = Np*Nnode;
Nmu  = Ne*Nnode;
Ntot = Nphi + Nmu;

idPhi = cell(1,Np);
for a = 1:Np
    idPhi{a} = ((a-1)*Nnode + (1:Nnode)).';
end

idMu = cell(1,Ne);
for l = 1:Ne
    idMu{l} = Nphi + ((l-1)*Nnode + (1:Nnode)).';
end

% ------------------------------------------------------------
% Tangents dp_i/dphi_alpha at fixed old geometry
% ------------------------------------------------------------
dpdphi = cell(Np,Np);

for alpha = 1:Np
    for ip = 1:Np
        dpdphi{alpha,ip} = MODEL.dpdphi(alpha,ip,phi_ref);
    end
end

% ------------------------------------------------------------
% B{ie,alpha} = dE_ie/dphi_alpha
%              = sum_ip e_ip_ie * dp_ip/dphi_alpha
%
% e is taken from coefficient state.
% dpdphi is taken from fixed old geometry.
% ------------------------------------------------------------
B = cell(Ne,Np);

for ie = 1:Ne
    for alpha = 1:Np

        tmp = zeros(ny,nx);

        for ip = 1:Np
            tmp = tmp + e0{ip}{ie} .* dpdphi{alpha,ip};
        end

        B{ie,alpha} = tmp;
    end
end

% ------------------------------------------------------------
% Sparse matrix allocation
% ------------------------------------------------------------
max_nnz = Nnode * ( ...
    5*Np + Np*Ne + ...
    Ne*(5 + 13*Ne + 13*Np) ...
    ) + 1000;

rows = zeros(max_nnz,1);
cols = zeros(max_nnz,1);
vals = zeros(max_nnz,1);
R    = zeros(Ntot,1);

k = 1;

% ============================================================
% 1. Allen-Cahn block
% ============================================================
for alpha = 1:Np

    row = idPhi{alpha}(idx_c);

    phi_a = phi_ref(:,:,alpha);

    % Old/reference Laplacian for RHS
    lap_phi_ref = laplacian_reflect(phi_a,dx,dy);

    % Original old-to-new RHS.
    % This remains old-time referenced, not correction-form.
    rhs = PARAM.LK .* lap_phi_ref + S_AC{alpha};

    R(row) = rhs(idx_c);

    LKc = PARAM.LK(idx_c);
    Aac = PARAM.A_ac(idx_c);

    % dphi_alpha block:
    % (1/dt)dphi - LK*Lap(dphi) + Aac*dphi
    cC = 1/dt + Aac + 2*LKc/dx2 + 2*LKc/dy2;
    cL = -LKc/dx2;
    cR = -LKc/dx2;
    cU = -LKc/dy2;
    cD = -LKc/dy2;

    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_c),cC);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_L),cL);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_R),cR);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_U),cU);
    [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_D),cD);

    % Coupling to dmu through domega_i = -e_i*dmu.
    %
    % Moving unknown term to left:
    %
    %   -L * sum_i dp_i/dphi_alpha * e_i * dmu
    %
    for ie = 1:Ne

        coeff_mu = zeros(ny,nx);

        for ip = 1:Np
            coeff_mu = coeff_mu - PARAM.L .* dpdphi{alpha,ip} .* e0{ip}{ie};
        end

        [rows,cols,vals,k] = add_block( ...
            rows,cols,vals,k,row,idMu{ie}(idx_c),coeff_mu(idx_c));
    end
end

% ============================================================
% 2. Cahn-Hilliard block with implicit LE closure
% ============================================================
%
% The unknown relation is:
%
%   E_new = E_ref + chi_coef*dmu + B_coef*dphi
%
% The CH equation is assembled as:
%
%   D(dmu) + A_E(chi*dmu + B*dphi)
%       = E_ref/dt - A_E(E_ref) - D(mu_ref)
%
% where:
%
%   A_E = 1/dt + M*kappa*biLaplacian
%   D   = -div(M grad)
%
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
    % RHS uses fixed reference state, not coefficient state
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
    % chi is from coefficient state
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
    % B uses e from coefficient state and dpdphi from reference state
    % --------------------------------------------------------
    for alpha = 1:Np

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

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_c),  a_C  .* B_c);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_L),  a_L  .* B_L);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_R),  a_R  .* B_R);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_U),  a_U  .* B_U);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_D),  a_D  .* B_D);

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_L2), a_L2 .* B_L2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_R2), a_R2 .* B_R2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_U2), a_U2 .* B_U2);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_D2), a_D2 .* B_D2);

        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_UR), a_UR .* B_UR);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_DR), a_DR .* B_DR);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_UL), a_UL .* B_UL);
        [rows,cols,vals,k] = add_block(rows,cols,vals,k,row,idPhi{alpha}(idx_DL), a_DL .* B_DL);
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
use_iter = false;

if isfield(NUM,'use_coupled_iter')
    use_iter = NUM.use_coupled_iter;
end

if use_iter
    tol   = 1e-8;
    maxit = 300;

    setup.type    = 'ilutp';
    setup.droptol = 1e-4;

    try
        [Lilu,Uilu] = ilu(A,setup);
        [sol,flag,relres,iter] = bicgstab(A,R,tol,maxit,Lilu,Uilu);
    catch
        [sol,flag,relres,iter] = bicgstab(A,R,tol,maxit);
    end

    if flag ~= 0 || any(~isfinite(sol))
        warning('Coupled bicgstab failed. Falling back to direct solve.');
        sol    = A \ R;
        flag   = 0;
        relres = norm(A*sol - R)/max(norm(R),eps);
        iter   = 0;
    end
else
    sol    = A \ R;
    flag   = 0;
    relres = norm(A*sol - R)/max(norm(R),eps);
    iter   = 0;
end

% ------------------------------------------------------------
% Unpack dphi and dmu
% ------------------------------------------------------------
dphi = zeros(ny,nx,Np);

for alpha = 1:Np
    dphi(:,:,alpha) = reshape(sol(idPhi{alpha}),ny,nx);
end

dmu = cell(1,Ne);

for ie = 1:Ne
    dmu{ie} = reshape(sol(idMu{ie}),ny,nx);
end

% ------------------------------------------------------------
% Update phi and p from fixed reference state
% ------------------------------------------------------------
STATE_NEW = STATE_REF;

STATE_NEW.phi = phi_ref + dphi;
STATE_NEW.phi = Normalize_Phi_Local(STATE_NEW.phi);
STATE_NEW.p   = Calc_p(MODEL,STATE_NEW.phi);

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
    for ip = 1:Np
        dp_ip = STATE_NEW.p(:,:,ip) - p_ref(:,:,ip);
        En = En + e0{ip}{ie} .* dp_ip;
    end

    STATE_NEW.E{ie} = En;
end

% Keep global mass exactly conserved relative to fixed old state
STATE_NEW.E = EnforceMeanE_Local(STATE_NEW.E,E_ref);

% ------------------------------------------------------------
% Predict omega for diagnostics only.
% Final LE_Run should overwrite this.
% ------------------------------------------------------------
STATE_NEW.omg = STATE_REF.omg;

for ip = 1:Np
    domega = zeros(ny,nx);

    for ie = 1:Ne
        domega = domega - e0{ip}{ie} .* dmu{ie};
    end

    STATE_NEW.omg(:,:,ip) = STATE_REF.omg(:,:,ip) + domega;
end

% Remove common omega offset
omg_mean = mean(STATE_NEW.omg,3);

for ip = 1:Np
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

DIAG.max_dmu     = max_dmu;
DIAG.matrix_size = Ntot;
DIAG.nnz         = nnz(A);

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


function phi = Normalize_Phi_Local(phi)

if exist('Norm_Phi','file') == 2
    phi = Norm_Phi(phi);
    return
end

phi = max(phi,0);

s = sum(phi,3);
s = max(s,eps);

for ip = 1:size(phi,3)
    phi(:,:,ip) = phi(:,:,ip) ./ s;
end

end


function E = EnforceMeanE_Local(E,E_old)

Ne = numel(E);

for ie = 1:Ne
    target_mean = mean(E_old{ie}(:));
    new_mean    = mean(E{ie}(:));
    E{ie}       = E{ie} + target_mean - new_mean;
end

end


function dmax = Max_Cell_Diff(A,B)

dmax = 0;

for i = 1:numel(A)
    if isempty(A{i}) || isempty(B{i})
        continue
    end

    dmax = max(dmax,max(abs(A{i}(:) - B{i}(:))));
end

end


function S = Blend_Tangent_State(Sold,Snew,alpha)

S = Snew;

% Damped fields used as thermodynamic coefficients
S.mu_e = Mix_Cell_Field(Sold.mu_e,Snew.mu_e,alpha);
S.E    = Mix_Cell_Field(Sold.E,   Snew.E,   alpha);
S.chi  = Mix_Cell_Field(Sold.chi, Snew.chi, alpha);

% Nested phase compositions: e{phase}{element}
S.e = Sold.e;

for ip = 1:numel(Sold.e)
    for ie = 1:numel(Sold.e{ip})
        S.e{ip}{ie} = (1-alpha)*Sold.e{ip}{ie} + alpha*Snew.e{ip}{ie};
    end
end

S.omg = (1-alpha)*Sold.omg + alpha*Snew.omg;

end


function C = Mix_Cell_Field(A,B,alpha)

C = A;

for i = 1:numel(A)
    C{i} = (1-alpha)*A{i} + alpha*B{i};
end

end