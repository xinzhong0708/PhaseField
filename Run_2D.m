%Clear and restart
clear;figure(1);clf;addpath('bin');addpath('ThermoData');addpath('Thermo');addpath('Maps');addpath('Thermo\Solutions')

load Map2d.mat
rng(1)
for ii = 1:5
% STATE.c{1}{ii}     =  STATE.c{1}{ii}+0.001*(rand(size(STATE.c{1}{ii}))-0.5);
% STATE.c{2}{ii}     =  STATE.c{2}{ii}+ repmat(0.002*(rand(1,length(STATE.c{2}{ii}))-0.5) , GRID.ny , 1);
end


%PHYSICS
PHYS.E_sc           =  E_sc;
PHYS.t_sc           =  1;                                                   % Time scale
PHYS.L_sc           =  1e-6;                                                % Length scale
PHYS.l              =  4*GRID.dx/PHYS.L_sc;                                 % interface thickness (m)
PHYS.sigma          =  0.4/(PHYS.E_sc*PHYS.L_sc);                           % surface energy (J/m^2)
PHYS.kappa          =  0e-7/(PHYS.E_sc*PHYS.L_sc^2);                        % 4th order term, can be set to 0 if no solvus
PHYS.D_esti         =  1e-12;
PHYS.chi_ref        =  1e-2;
PHYS.M0             =  PHYS.D_esti*PHYS.t_sc/PHYS.L_sc^2*PHYS.chi_ref;
PHYS.M_phs          = [PHYS.M0/100 PHYS.M0/100 PHYS.M0/100 PHYS.M0/100 PHYS.M0/100
                       PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0
                       PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0
                       PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0
                       PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0];
% PHYS.M_phs          = [PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0
%                        PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0
%                        PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0     PHYS.M0];
PHYS.m              =  6*PHYS.sigma/PHYS.l;
PHYS.kap            =  3/4*PHYS.sigma*PHYS.l;
PHYS.dceq           =  0.5;
PHYS.L              =  4*PHYS.m/3/PHYS.kap/(PHYS.dceq^2/PHYS.M0);
PHYS.eta            =  eta;

%NUMERICS
NUM.dt_phy          =   2e-4/PHYS.t_sc;
NUM.dt_max          =    1e1/PHYS.t_sc;
NUM.dt_min          =  1e-16/PHYS.t_sc; 
NUM.t_tot           =    1e5/PHYS.t_sc;
NUM.dE_target       =  3e-2;
NUM.dp_target       =  3e-2;
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
NUM.int_damp        =  0.1;
NUM.use_Jphi        =  1;
NUM.CHLE_p_cut      =  1e-8;
NUM.CHLE_band_thick =  16;
NUM.CHLE_res_rel    =  [];
NUM.phi_mask_cut    =  1e-8;

NUM.linear_solver   = 'bicgstab_ilu';
NUM.linear_tol      =  1e-9;
NUM.linear_maxit    =  500;
NUM.ilu_reuse       =  1;
NUM.ilu_reuse_steps =  10;
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
PARAM.CS_chi_floor  =  1e-7;
PARAM.CS_chi_cap    =  1e-2;

PARAM.use_kappa_c   =  1; 
PARAM.L_fac         =  1;
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
disp([mean(STATE.E{1},'all') mean(STATE.E{2},'all') mean(STATE.E{3},'all') mean(STATE.E{4},'all') mean(STATE.E{end},'all')])

% load 180

% NUM.linear_solver   = 'direct';
% NUM.phi_mask_thick  =  5;
%TIME LOOP
% PARAM.L(:)  = mean(PARAM.L(:));
% PARAM.Lm(:) = mean(PARAM.Lm(:));
% PARAM.LK(:) = mean(PARAM.LK(:));

% load 2800
% NUM.linear_solver   = 'bicgstab_ilu';


for it = 1:20000
    if mod(it,100)==0
        save(num2str(it))
    end


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
    [STATE_RAW,DIAG_RAW] =    PF_Coupled_ACCH_LETangent_CS(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);
    t_ACCH               =    toc(t);

    % SECOND LOCAL EQUILIBRIUM
    t                    =    tic;
    STATE_LE0            =    LE_Run_Mode_New(STATE_RAW,PARAM,MODEL);
    STATE_LE0            =    Extend_AbsentPhaseC_Rim(STATE_LE0,PARAM);
    t_LE2                =    toc(t);

    % Fixed-p chemical corrector
    t                    =    tic;
    STATE_TRIAL          =    PF_CH_LECorrector_FixedP_Band_CS(STATE_OLD,STATE_LE0,PARAM,MODEL,GRID,PHYS,NUM);
    t_CHLE               =    toc(t);

    % TIME STEP UPDATE
    dt_try               =    NUM.dt_phy;
    time_old             =    NUM.time;
    [STATE,NUM]          =    Update_TimeStep_Soft(STATE,STATE_TRIAL,PARAM,MODEL,NUM);
    
    % UPDATE M
    PARAM                =    Compute_M(STATE,PARAM,MODEL,PHYS);

    % UPDATE L
    PARAM                =    Compute_L(STATE,PARAM,PHYS);

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

    if mod(it,2)==0
        disp(PHASE(end,end))
        PF_Plot([3,3,1],'E3',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,2],'mu_e1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,3],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,4],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,5],'omg12',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,6],'omg23',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,7],'omg34',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,8],'c51',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,9],'PhaseStack',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        drawnow
    end


    % if mod(it,2)==0
    %     disp(NUM.dt_phy)
    %     disp([mean(STATE.E{1},'all') mean(STATE.E{2},'all') mean(STATE.E{3},'all') ])
    %     disp(PHASE(end,:))
    %     % subplot(331);plot(GRID.x,STATE.E{1}(3,:),GRID.x,STATE.E{2}(3,:),GRID.x,STATE.E{3}(3,:),GRID.x,STATE.E{4}(3,:),GRID.x,STATE.E{end}(3,:));title('E1')
    %     % subplot(331);plot(GRID.x,STATE.E{1}(3,:),GRID.x,STATE.E{2}(3,:),GRID.x,STATE.E{3}(3,:));title('E1')
    %     subplot(331);plot(GRID.x,STATE.E{3}(3,:));title('E1')
    %     subplot(332);plot(STATE.mu_e{1}(3,:));title('mu_e')
    %     subplot(333);plot(DTPHY,'b.');title('dt')
    %     % subplot(336);plot(GRID.x,CHI');title('chi')
    %     subplot(334);plot(GRID.x,STATE.p(3,:,1),'.-',GRID.x,STATE.p(3,:,2),'.-',GRID.x,STATE.p(3,:,end-1),'.-',GRID.x,STATE.p(3,:,end),'.-');title('p2')        
    %     subplot(335);plot(GRID.x,STATE.omg(3,:,1)-STATE.omg(3,:,2),GRID.x,STATE.omg(3,:,end-1)-STATE.omg(3,:,end),'.-');title('domg12')
    %     subplot(337);plot(GRID.x,STATE.c{1}{1}(3,:),GRID.x,STATE.c{2}{1}(3,:),GRID.x,STATE.c{3}{1}(3,:));title('c11')
    %     % subplot(339);plot(GRID.x,diff_E{1}(2,:) , GRID.x,diff_E{2}(2,:) , GRID.x,diff_E{3}(2,:));
    %     subplot(338);plot(TIME,PHASE);
    %     drawnow
    % end

    % if mod(it,2)==0
    %     disp(NUM.dt_phy)
    %     disp(PHASE(end,:))
    %     PF_Plot([2,2,1],'E1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([2,2,2],'mu_e1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([2,2,3],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     subplot(224);pcolor(GRID.x,GRID.y,CHI);colorbar;shading interp;axis equal
    %     drawnow
    % end

end









