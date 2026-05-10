%Clear and restart
clear;figure(1);clf;addpath([cd,'\bin']);addpath([cd,'\Thermo']);addpath([cd,'\Thermo\Solutions'])

%Load map
load Map2d.mat

%PHYSICS
PHYS.E_sc          =  E_sc;
PHYS.t_sc          =  1;                                                   % Time scale
PHYS.L_sc          =  1;                                                   % Length scale
PHYS.l             =  Lx/200/L_sc;                                         % interface thickness (m)
PHYS.sigma         =  1.0/PHYS.E_sc*PHYS.L_sc^2;                           % surface energy (J/m^2)
PHYS.kappa         =  0e-7/(PHYS.E_sc*PHYS.L_sc^2);                        % 4th order term, can be set to 0 if no solvus
PHYS.D_esti        =  1e-12;
PHYS.chi_ref       =  1e-2;
PHYS.M0            =  PHYS.D_esti*PHYS.t_sc/PHYS.L_sc^2*PHYS.chi_ref;
PHYS.m             =  6*PHYS.sigma/PHYS.l;
PHYS.kap           =  3/4*PHYS.sigma*PHYS.l;
PHYS.dceq          =  0.2;
PHYS.L             =  4*PHYS.m/3/PHYS.kap/(PHYS.dceq^2/PHYS.M0)/50;
PHYS.eta           =  eta;

%NUMERICS
NUM.dt_phy         =   1e-5/PHYS.t_sc;
NUM.dt_max         =    1e5/PHYS.t_sc;
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

%TIME STEPPING
for it = 1:1e5
    if mod(it,100)==0
        save(num2str(it))
    end
    
    %OLD STATE
    STATE_OLD            =    STATE;
    STATE_TRIAL          =    STATE;

    %===========================TRIAL STEPS================================
    PARAM.A_ac           =    Calc_Aac_FrozenOmega(STATE,PARAM,MODEL,3,1e-6,0,[]);

    %ALLEN CAHN
    STATE_TRIAL          =    Calc_S_AllenCahn(STATE_TRIAL,PARAM,MODEL);
    STATE_TRIAL          =    PF_IMEX_Solver2D_AllenCahn_Stab(STATE_TRIAL,MODEL,PARAM,GRID,NUM,1);

    %INTERFACE PENALTY
    PARAM.eta            =    Eta_Damping(STATE_TRIAL.p,PHYS.eta,0.5*PHYS.eta); 

    %LOCAL EQUILIBRIUM
    STATE_TRIAL          =    LE_Run(STATE_TRIAL,PARAM,MODEL);

    %CAHN HILLARD
    STATE_TRIAL          =    Calc_S_Diffusion(STATE_TRIAL,STATE_OLD,NUM);
    STATE_TRIAL          =    PF_IMEX_Solver2D_Diffusion_MuOnly(STATE_TRIAL,PARAM,GRID,PHYS,NUM);

    %TIME STEP UPDATE
    [STATE,NUM,DIAG]     =    Update_TimeStep_Soft(STATE,STATE_TRIAL,PARAM,MODEL,NUM);


    %Plotting
    TIME(it)  = NUM.time;
    PHASE1(it)= mean(mean(STATE.p(:,:,1)));
    PHASE2(it)= mean(mean(STATE.p(:,:,2)));

    %1D
    if mod(it,10)==0
        disp(NUM.dt_phy)
        disp([mean(mean(STATE.p(:,:,1))),mean(mean(STATE.p(:,:,2))),mean(STATE.E{1}(:))])
        subplot(331);plot(GRID.x,STATE.E{1}(3,:),GRID.x,STATE.E{2}(3,:),GRID.x,STATE.E{3}(3,:));title('E1')
        subplot(332);plot(STATE.mu_e{1}(3,:));title('mu_e')
        subplot(333);plot(it,NUM.dt_phy,'b.');hold on;title('dt')
        subplot(334);plot(GRID.x,STATE.phi(3,:,2),'.-',GRID.x,STATE.phi(3,:,3),'.-',GRID.x,STATE.phi(3,:,4),'.-');title('p2')        
        subplot(335);plot(STATE.omg(3,:,1)-STATE.omg(3,:,2));title('domg12')
        subplot(336);plot(PARAM.eta(3,:));title('eta')
        subplot(337);plot(STATE.c{1}{1}(3,:));title('c11')
        subplot(338);plot(STATE.c{2}{1}(3,:));title('c21')
        subplot(339);plot(TIME,PHASE2,'.-');title('Phase2')
        drawnow
    end

end

