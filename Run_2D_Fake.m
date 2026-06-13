%Clear and restart
clear;figure(1);clf;addpath('bin');addpath('ThermoData');addpath('Thermo');addpath('Maps');addpath('Thermo\Solutions')

load Map2d.mat

%PHYSICS
PHYS.E_sc           =  E_sc;
PHYS.t_sc           =  1;                                                   % Time scale
PHYS.L_sc           =  1;                                                   % Length scale
PHYS.l              =  3*GRID.dx/PHYS.L_sc;                                 % interface thickness (m)
PHYS.sigma          =  0.50/(PHYS.E_sc*PHYS.L_sc);                          % surface energy (J/m^2)
PHYS.kappa          =  0e-7/(PHYS.E_sc*PHYS.L_sc^2);                        % 4th order term, can be set to 0 if no solvus
PHYS.D_esti         =  1e-12;
PHYS.D_esti         =  1e-6;
PHYS.chi_ref        =  1e-2;
PHYS.M0             =  PHYS.D_esti*PHYS.t_sc/PHYS.L_sc^2*PHYS.chi_ref;
PHYS.M_phs          =  PHYS.M0;
PHYS.m              =  6*PHYS.sigma/PHYS.l;
PHYS.kap            =  3/4*PHYS.sigma*PHYS.l;
PHYS.dceq           =  0.5;
PHYS.L              =  4*PHYS.m/3/PHYS.kap/(PHYS.dceq^2/PHYS.M0);
PHYS.eta            =  eta;

%NUMERICS
NUM.dt_phy          =   1e-2/PHYS.t_sc;
NUM.dt_max          =    1e3/PHYS.t_sc;
NUM.dt_min          =  1e-16/PHYS.t_sc; 
NUM.t_tot           =    1e5/PHYS.t_sc;
NUM.dE_target       =  1e-2;
NUM.dp_target       =  1e-2;
NUM.dmu_target      =  1e2;
NUM.time            =  0;
NUM.dt_good_count   =  0;
NUM.dt_grow_after   =  8;
NUM.dt_grow_fac     =  1.15;
NUM.dt_shrink_fac   =  0.5;
NUM.err_grow        =  0.25;
NUM.phi_mask_cut    =  1e-8;
NUM.phi_mask_thick  =  3;
NUM.norm_phi        =  1;
NUM.cut_phi         =  0;
NUM.norm_E          =  1;
NUM.int_damp        =  0.5;
NUM.use_Jphi        =  0;
NUM.CHLE_p_cut      =  1e-8;
NUM.CHLE_band_thick =  16;
NUM.CHLE_res_rel    =  [];
NUM.phi_mask_cut    =  1e-8;

NUM.linear_solver   = 'bicgstab_ilu';
% NUM.linear_solver   = 'direct';
NUM.linear_tol      =  1e-9;
NUM.linear_maxit    =  500;
NUM.ilu_reuse       =  1;
NUM.ilu_reuse_steps =  2;
NUM.ilu_rebuild     =  1;
NUM.direct_fallback =  1;

NUM.ilu_reuse_ACCH          =  1;
NUM.ilu_reuse_steps_ACCH    =  10;
NUM.ACCH_mask_update        =  NUM.ilu_reuse_steps_ACCH;
NUM.ilu_cache_check_pattern =  0;

%GRIDS
GRID.dx             =  GRID.dx/PHYS.L_sc;
GRID.dy             =  GRID.dy/PHYS.L_sc;
GRID.x              =  GRID.x /PHYS.L_sc;
GRID.y              =  GRID.y /PHYS.L_sc;
GRID.nx             =  GRID.nx;
GRID.ny             =  GRID.ny;

%MODELS
MODEL.dgdphi        =  @(phi) 2*PHYS.m*phi.*(phi - 1).^2 + PHYS.m*phi.^2.*(2.*phi - 2);
MODEL.pars          =  pars;
MODEL.p_fun         =  @(a,phi)   phi(:,:,a).^2./sum(phi.^2,3);
MODEL.dpdphi        =  @(a,b,phi) (a==b)*2*phi(:,:,b)./sum(phi.^2,3) - 2*phi(:,:,a).*phi(:,:,b).^2./sum(phi.^2,3).^2;

%PARAMETERS
PARAM.L             =  PHYS.L*ones(GRID.ny,GRID.nx);
PARAM.Lm            =  PHYS.L*PHYS.m.*ones(GRID.ny,GRID.nx);
PARAM.LK            =  PHYS.L*PHYS.kap.*ones(GRID.ny,GRID.nx);
PARAM.Np            =  length(STATE.c);
PARAM.Ne            =  length(STATE.E);
PARAM.M             =  repmat({PHYS.M0*ones(GRID.ny,GRID.nx)},1,PARAM.Ne);

PARAM.kappa_phase   =  PHYS.kappa .* cellfun(@(x) size(x.n,1) > 1, pars);   % kappa nonzero automatically for phase with >1 endmembers
PARAM.use_CS_chi    =  1;
PARAM.CS_chi_floor  =  1e-9;
PARAM.CS_chi_cap    =  1e-2;

PARAM.use_kappa_c   =  1; 
PARAM.L_fac         =  0.5;
PARAM.LE_mode       =  'LE';

%STATES
STATE.c             =  STATE.c;
STATE.e             =  STATE.e;
STATE.E             =  STATE.E;
STATE.mu_e          =  STATE.mu_e;
STATE.chi           =  STATE.chi;
STATE.omg           =  zeros(GRID.ny,GRID.nx,Np);
STATE.phi           =  STATE.phi;
STATE.p             =  Calc_p(MODEL,STATE.phi);
STATE.mask          =  ones(GRID.ny,GRID.nx,Np);
STATE.LE_state      = [   ];

%DISPLAY COMPOSITION
% disp([mean(STATE.E{1},'all') mean(STATE.E{2},'all') mean(STATE.E{3},'all') mean(STATE.E{4},'all') mean(STATE.E{end},'all')])

% load test
PARAM.aniso_delta = 0.15;
PARAM.aniso_nfold = 4;

PARAM.theta_grain = zeros(1,size(STATE.p,3));

for it = 1:50000
    
    PARAM.kappa_eff      =    PHYS.kappa*ones(GRID.ny,GRID.nx);

    % DIRECT COUPLED SOLVER WITH MU-MATCHED E RECOVERY
    STATE_OLD            =    STATE;
    t_step               =    tic;

    % DAMPED ETA
    PARAM.eta            =    Eta_Damping(STATE_OLD.p,PHYS.eta,NUM.int_damp*PHYS.eta);

    % FIRST LOCAL EQUILIBRIUM
    t                    =    tic;
    STATE_OLD            =    LE_Run_Mode_New(STATE_OLD,PARAM,MODEL);
    STATE_OLD            =    Extend_AbsentPhaseC_Rim(STATE_OLD,PARAM);
    t_LE1                =    toc(t);

    % Full AC + CH COUPLED PREDICTOR
    t                    =    tic;
    PARAM                =    Update_PF_SolverMasks(PARAM,STATE_OLD,MODEL,GRID,PHYS,NUM,'ACCH');
    PARAM                =    Calc_AC_Anisotropy_Simple(STATE_OLD,PARAM,GRID);
    [STATE_RAW,DIAG_RAW] =    PF_Coupled_ACCH_LETangent_CS_offdiagM(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);
    t_ACCH               =    toc(t);

    % SECOND LOCAL EQUILIBRIUM
    t                    =    tic;
    STATE_LE0            =    LE_Run_Mode_New(STATE_RAW,PARAM,MODEL);
    STATE_LE0            =    Extend_AbsentPhaseC_Rim(STATE_LE0,PARAM);
    t_LE2                =    toc(t);

    % Fixed-p chemical corrector
    t                    =    tic;
    STATE_TRIAL          =    PF_CH_LECorrector_FixedP_Band_CS_offdiagM(STATE_OLD,STATE_LE0,PARAM,MODEL,GRID,PHYS,NUM);
    t_CHLE               =    toc(t);

    % TIME STEP UPDATE
    dt_try               =    NUM.dt_phy;
    time_old             =    NUM.time;
    [STATE,NUM]          =    Update_TimeStep_Soft(STATE,STATE_TRIAL,PARAM,MODEL,NUM);
    
    % UPDATE M
    % PARAM                =    Compute_M(STATE,PARAM,MODEL,PHYS);

    % UPDATE L
    % PARAM                =    Compute_L(STATE,PARAM,PHYS);

    % PRINT TIME
    t_total              =    toc(t_step);
    disp(['Total time:',num2str(t_total),' LE1:',num2str(t_LE1),' ACCH:',num2str(t_ACCH),' LE2:',num2str(t_LE2),' CH:',num2str(t_CHLE)])


    %Plotting
    TIME(it)             =    NUM.time;
    DTPHY(it)            =    NUM.dt_phy;
    phase_ids            =    unique(MODEL.phase_index,'stable');
    PHASE(it,:)          =    zeros(1,numel(phase_ids));
    for iph = 1:numel(phase_ids)
        grains = find(MODEL.phase_index == phase_ids(iph));
        PHASE(it,iph) = mean(sum(STATE.p(:,:,grains),3),'all');
    end


    if mod(it,10)==0
        disp(PHASE(end,end))
        PF_Plot([2,3,1],'E1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([2,3,2],'mu_e1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([2,3,3],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([2,3,4],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([2,3,5],'omg12',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        subplot(236);plot(TIME,PHASE(:,1))
        drawnow
    end

end









%Analytical
xl_eq  = 0.5;
x0     = 0.3;
xs_eq  = 0;
fun    = @(lam) (xl_eq - xs_eq).*lam.*sqrt(pi).*exp(lam.^2).*erfc(lam) - (x0 - xl_eq);

lambda = fzero(fun,[-10,10]);
t      = linspace(0,100000,1000);
s      = 2*lambda*sqrt(PHYS.D_esti*t);
clf
plot(t,abs(s)+0.1)
hold on
plot(TIME,PHASE(:,1),'k--')


