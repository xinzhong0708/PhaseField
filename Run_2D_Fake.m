%Clear and restart
clear;figure(3);clf;addpath('bin');addpath('ThermoData');addpath('Thermo');addpath('Maps');addpath('Thermo\Solutions')

load Map2d.mat

%PHYSICS
PHYS.E_sc          =  E_sc;
PHYS.t_sc          =  1;                                                   % Time scale
PHYS.L_sc          =  1;                                                   % Length scale
PHYS.l             =  5*GRID.dx/PHYS.L_sc;                                 % interface thickness (m)
PHYS.sigma         =  0.4/(PHYS.E_sc*PHYS.L_sc);                           % surface energy (J/m^2)
PHYS.kappa         =  0e-6/(PHYS.E_sc*PHYS.L_sc^2);                        % 4th order term, can be set to 0 if no solvus
PHYS.D_esti        =  1e-6;
PHYS.chi_ref       =  1e-2;
PHYS.M0            =  PHYS.D_esti*PHYS.t_sc/PHYS.L_sc^2*PHYS.chi_ref;
PHYS.m             =  6*PHYS.sigma/PHYS.l;
PHYS.kap           =  3/4*PHYS.sigma*PHYS.l;
PHYS.dceq          =  0.5;
PHYS.L             =  0.01*4*PHYS.m/3/PHYS.kap/(PHYS.dceq^2/PHYS.M0)/1;
PHYS.eta           =  eta;

%NUMERICS
NUM.dt_phy         =   1e-3/PHYS.t_sc;
NUM.dt_max         =    1e5/PHYS.t_sc;
NUM.dt_min         =  1e-16/PHYS.t_sc; 
NUM.t_tot          =    1e5/PHYS.t_sc;
NUM.dE_target      =  2e-2;
NUM.dp_target      =  2e-2;
NUM.dmu_target     =  1e2;
NUM.time           =  0;
NUM.dt_good_count  =  0;
NUM.dt_grow_after  =  8;
NUM.dt_grow_fac    =  1.15;
NUM.dt_shrink_fac  =  0.5;
NUM.err_grow       =  0.25;
NUM.phi_mask_cut   =  1e-8;
NUM.phi_mask_thick =  2;
NUM.norm_phi       =  1;
NUM.cut_phi        =  0;
NUM.norm_E         =  1;
NUM.int_damp       =  0.5;
NUM.kappa_p_cut    =  1e-6;
NUM.use_Jphi       =  0;
NUM.CHLE_p_cut     =  1e-8;
NUM.CHLE_band_thick=  5;
NUM.CHLE_res_rel   =  [];
PARAM.use_CS_chi   =  1;
PARAM.CS_chi_floor =  1e-9;
PARAM.CS_chi_cap   =  1e-2;
PARAM.use_kappa_c  =  1; 

%GRIDS
GRID.dx            =  GRID.dx/PHYS.L_sc;
GRID.dy            =  GRID.dy/PHYS.L_sc;
GRID.x             =  GRID.x /PHYS.L_sc;
GRID.y             =  GRID.y /PHYS.L_sc;
GRID.nx            =  GRID.nx;
GRID.ny            =  GRID.ny;

%MODELS
MODEL.dgdphi       =  @(phi) 2*PHYS.m*phi.*(phi - 1).^2 + PHYS.m*phi.^2.*(2.*phi - 2);
MODEL.pars         =  pars;
MODEL.p_fun        =  @(a,phi)   phi(:,:,a).^2./sum(phi.^2,3);
MODEL.dpdphi       =  @(a,b,phi) (a==b)*2*phi(:,:,b)./sum(phi.^2,3) - 2*phi(:,:,a).*phi(:,:,b).^2./sum(phi.^2,3).^2;

%PARAMETERS
PARAM.L            =  PHYS.L*ones(GRID.ny,GRID.nx);
PARAM.Lm           =  PHYS.L*PHYS.m.*ones(GRID.ny,GRID.nx);
PARAM.LK           =  PHYS.L*PHYS.kap.*ones(GRID.ny,GRID.nx);
PARAM.Np           =  length(STATE.c);
PARAM.Ne           =  length(STATE.E);
PARAM.M            =  repmat({PHYS.M0*ones(GRID.ny,GRID.nx)},1,PARAM.Ne);
PARAM.kappa_phase  =  PHYS.kappa .* cellfun(@(x) size(x.n,1) > 1, pars);   % kappa nonzero automatically for phase with >1 endmembers

%STATES
STATE.c            =  STATE.c;
STATE.e            =  STATE.e;
STATE.E            =  STATE.E;
STATE.mu_e         =  STATE.mu_e;
STATE.chi          =  STATE.chi;
STATE.omg          =  zeros(GRID.ny,GRID.nx,Np);
STATE.phi          =  STATE.phi;
STATE.p            =  Calc_p(MODEL,STATE.phi);
STATE.mask         =  ones(GRID.ny,GRID.nx,Np);
STATE.LE_state     = [   ];

% load 10000
% NUM.dE_target      =  1e-2;
% NUM.dp_target      =  1e-2;
% load Run_L0_1

for it = 1:1e5
    if mod(it,100)==0
        save(num2str(it))
    end

    tic


    % DIRECT COUPLED SOLVER WITH MU-MATCHED E RECOVERY
    STATE_OLD            =    STATE;

    % DAMPED ETA
    PARAM.eta            =    Eta_Damping(STATE_OLD.p,PHYS.eta,NUM.int_damp*PHYS.eta);

    % Old-state local equilibrium
    PARAM.LE_mode        =   'LE';
    STATE_OLD            =    LE_Run_Mode_New(STATE_OLD,PARAM,MODEL);
    STATE_OLD            =    Extend_AbsentPhaseC_Rim(STATE_OLD,PARAM);

    % Full AC + CH coupled predictor
    STATE_RAW            =    PF_Coupled_ACCH_LETangent_CS(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);

    % Re-run fixed-E LE on corrected E
    % PARAM.LE_mode        =   'GP';
    STATE_LE0            =    LE_Run_Mode_New(STATE_RAW,PARAM,MODEL);
    STATE_LE0            =    Extend_AbsentPhaseC_Rim(STATE_LE0,PARAM);

    % Fixed-p chemical corrector
    STATE_CORR           =    PF_CH_LECorrector_FixedP_Band_CS(STATE_OLD,STATE_LE0,MODEL,PARAM,GRID,PHYS,NUM);

    % Final nonlinear fixed-E state
    STATE_TRIAL          =    LE_Run_Mode_New(STATE_CORR,PARAM,MODEL);
    STATE_TRIAL          =    Extend_AbsentPhaseC_Rim(STATE_TRIAL,PARAM);

    % TIME STEP UPDATE
    [STATE,NUM]          =    Update_TimeStep_Soft(STATE,STATE_TRIAL,PARAM,MODEL,NUM);

    
    % Parallel
    for ip = 1:size(STATE.p,3)
        STATE.p(:,:,ip)  =    repmat(mean(STATE.p(:,:,ip)) , GRID.ny , 1);
    end
    for ie = 1:length(STATE.E)
        STATE.E{ie}      =    repmat(mean(STATE.E{ie}) , GRID.ny , 1);
    end

    % CHECK
    dE{1}                =    STATE_TRIAL.E{1}-STATE.E{1};
  

    toc

    %Check E c consistency
    e_check = Calc_e(MODEL.pars,STATE_TRIAL.c);
    E_check = Calc_E_Tot(e_check,STATE_TRIAL.p);
    for ie = 1:numel(STATE_TRIAL.E)
        diff_E{ie} = E_check{ie} - STATE_TRIAL.E{ie};
        cons_E{ie} = E_check{ie} - STATE_TRIAL.E{ie} + STATE_TRIAL.mu_e{ie}./PARAM.eta;
    end

    
    %Plotting
    TIME(it)             =    NUM.time;
    DTPHY(it)            =    NUM.dt_phy;
    phase_ids            =    unique(MODEL.phase_index,'stable');
    PHASE(it,:)          =    zeros(1,numel(phase_ids));
    for iph = 1:numel(phase_ids)
        grains = find(MODEL.phase_index == phase_ids(iph));
        PHASE(it,iph) = mean(sum(STATE.p(:,:,grains),3),'all');
    end




    if mod(it,5)==0
        %Calculate G profile
        NUM.dt_phy
        %Plot
        cc       = linspace(0,1,1000);
        R1       = PhaseThermo(MODEL.pars{1},{cc});
        R2       = PhaseThermo(MODEL.pars{2},{cc});
        c1_pt    = STATE.c{1}{1}(STATE.p(:,:,1)>1e-6);
        gpt1     = PhaseThermo(MODEL.pars{1},{c1_pt});
        c2_pt    = STATE.c{2}{1}(STATE.p(:,:,2)>1e-6);
        gpt2     = PhaseThermo(MODEL.pars{2},{c2_pt});
        mid      = GRID.nx/2;
        mid      = 2;
        subplot(331);plot(GRID.x,STATE.E{1}(mid,:));title('E1')
        subplot(332);plot(STATE.mu_e{1}(mid,:));title('mu_e')
        subplot(333);plot(DTPHY,'b.');title('dt')
        subplot(336);plot(cc,R1.g,cc,R2.g,c1_pt,gpt1.g,'k.' ,c2_pt,gpt2.g,'m.')
        ylim([-5,10])
        subplot(334);plot(GRID.x,STATE.p(mid,:,1),'.-',GRID.x,STATE.p(mid,:,2),'.-');title('p2')        
        subplot(335);plot(GRID.x,STATE.omg(mid,:,1)-STATE.omg(mid,:,2));title('domg12')
        subplot(337);plot(GRID.x,STATE.c{1}{1}(mid,:));title('c11')
        subplot(338);plot(STATE.E{1}(1,:)-STATE.E{1}(2,:))
        % plot(GRID.x,STATE.c{2}{1}(mid,:));title('c11')
        subplot(339);plot(TIME,PHASE(:,1));
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



