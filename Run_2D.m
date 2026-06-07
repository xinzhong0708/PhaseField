%Clear and restart
clear;figure(1);clf;addpath('bin');addpath('ThermoData');addpath('Thermo');addpath('Maps');addpath('Thermo\Solutions')

load Map2d.mat
rng(1)
for ii = 1:5
% STATE.c{1}{ii}     =  STATE.c{1}{ii}+0.005*(rand(size(STATE.c{1}{ii}))-0.5);
STATE.c{2}{ii}     =  STATE.c{2}{ii}+ repmat(0.002*(rand(1,length(STATE.c{2}{ii}))-0.5) , GRID.ny , 1);
end

STATE.e            =  Calc_e(MODEL.pars,STATE.c);
STATE.E            =  Calc_E_Tot(STATE.e,STATE.p);

%PHYSICS
PHYS.E_sc          =  E_sc;
PHYS.t_sc          =  1;                                                   % Time scale
PHYS.L_sc          =  1e-6;                                                % Length scale
PHYS.l             =  3*GRID.dx/PHYS.L_sc;                                 % interface thickness (m)
PHYS.sigma         =  0.4/(PHYS.E_sc*PHYS.L_sc);                           % surface energy (J/m^2)
PHYS.kappa         =  1e-6/(PHYS.E_sc*PHYS.L_sc^2);                        % 4th order term, can be set to 0 if no solvus
PHYS.D_esti        =  1e-12;
PHYS.chi_ref       =  1e-2;
PHYS.M0            =  PHYS.D_esti*PHYS.t_sc/PHYS.L_sc^2*PHYS.chi_ref;
PHYS.m             =  6*PHYS.sigma/PHYS.l;
PHYS.kap           =  3/4*PHYS.sigma*PHYS.l;
PHYS.dceq          =  1;
PHYS.L             =  4*PHYS.m/3/PHYS.kap/(PHYS.dceq^2/PHYS.M0)/20;
PHYS.eta           =  eta;

%NUMERICS
NUM.dt_phy         =   1e-3/PHYS.t_sc;
NUM.dt_max         =    1e1/PHYS.t_sc;
NUM.dt_min         =  1e-16/PHYS.t_sc; 
NUM.t_tot          =    1e5/PHYS.t_sc;
NUM.dE_target      =  3e-2;
NUM.dp_target      =  3e-2;
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
NUM.use_Jphi       =  1;
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

%DISPLAY COMPOSITION
disp([mean(STATE.E{1},'all') mean(STATE.E{2},'all') mean(STATE.E{3},'all') mean(STATE.E{4},'all') mean(STATE.E{end},'all')])

% load temp
% NUM.dE_target      =  2e-2;
% NUM.dp_target      =  2e-2;
% NUM.int_damp       =  0.35;
% NUM.use_order_cache= 1;

for it = 1:1e5
    if mod(it,50)==0
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
    PARAM                =    Calc_Kappa_WScale_InterfaceRamp(STATE_OLD,PARAM,MODEL,PHYS,NUM);
    STATE_RAW            =    PF_Coupled_ACCH_LETangent_CS(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);

    % Re-run fixed-E LE on corrected E
    PARAM.LE_mode        =   'LE';
    STATE_LE0            =    LE_Run_Mode_New(STATE_RAW,PARAM,MODEL);
    STATE_LE0            =    Extend_AbsentPhaseC_Rim(STATE_LE0,PARAM);

    % Fixed-p chemical corrector
    STATE_CORR           =    PF_CH_LECorrector_FixedP_Band_CS(STATE_OLD,STATE_LE0,MODEL,PARAM,GRID,PHYS,NUM);

    % Final nonlinear fixed-E state
    STATE_TRIAL          =    LE_Run_Mode_New(STATE_CORR,PARAM,MODEL);
    STATE_TRIAL          =    Extend_AbsentPhaseC_Rim(STATE_TRIAL,PARAM);

    % TIME STEP UPDATE
    [STATE,NUM]          =    Update_TimeStep_Soft(STATE,STATE_TRIAL,PARAM,MODEL,NUM);



    % CHECK
    dE{1}                =    STATE_TRIAL.E{1}-STATE.E{1};
    dE{2}                =    STATE_TRIAL.E{2}-STATE.E{2};
    dE{3}                =    STATE_TRIAL.E{3}-STATE.E{3};
    dE{4}                =    STATE_TRIAL.E{4}-STATE.E{4};
  

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


    %Calculate chi
    for i = 1:GRID.nx
        for j = 1:GRID.ny
            chi = [];
            for ii = 1:5
                for jj = 1:5
                    chi(ii,jj) = STATE.chi{ii,jj}(j,i);
                end
            end
            CHI(j,i) = det(chi);
        end
    end

    % if mod(it,2)==0
    %     disp(NUM.dt_phy)
    %     disp(PHASE(end,:))
    %     PF_Plot([3,3,1],'E3',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,2],'mu_e1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,3],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,4],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,5],'omg12',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,6],'omg13',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,7],'omg23',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     subplot(337);pcolor(GRID.x,GRID.y,dE{1});colorbar;shading interp;axis equal
    %     subplot(338);pcolor(GRID.x,GRID.y,CHI);colorbar;shading interp;axis equal
    %     PF_Plot([3,3,9],'PhaseStack',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     drawnow
    % end


    if mod(it,2)==0
        disp(NUM.dt_phy)
        disp([mean(STATE.E{1},'all') mean(STATE.E{2},'all') mean(STATE.E{3},'all') ])
        disp(PHASE(end,:))
        % subplot(331);plot(GRID.x,STATE.E{1}(3,:),GRID.x,STATE.E{2}(3,:),GRID.x,STATE.E{3}(3,:),GRID.x,STATE.E{4}(3,:),GRID.x,STATE.E{end}(3,:));title('E1')
        % subplot(331);plot(GRID.x,STATE.E{1}(3,:),GRID.x,STATE.E{2}(3,:),GRID.x,STATE.E{3}(3,:));title('E1')
        subplot(331);plot(GRID.x,STATE.E{3}(3,:));title('E1')
        subplot(332);plot(STATE.mu_e{1}(3,:));title('mu_e')
        subplot(333);plot(DTPHY,'b.');title('dt')
        subplot(336);plot(GRID.x,CHI');title('chi')
        subplot(334);plot(GRID.x,STATE.p(3,:,1),'.-',GRID.x,STATE.p(3,:,2),'.-',GRID.x,STATE.p(3,:,end-1),'.-',GRID.x,STATE.p(3,:,end),'.-');title('p2')        
        subplot(335);plot(GRID.x,STATE.omg(3,:,1)-STATE.omg(3,:,2),GRID.x,STATE.omg(3,:,end-1)-STATE.omg(3,:,end),'.-');title('domg12')
        subplot(337);plot(GRID.x,STATE.c{1}{1}(3,:),GRID.x,STATE.c{2}{1}(3,:),GRID.x,STATE.c{3}{1}(3,:));title('c11')
        subplot(339);plot(GRID.x,dE{1}(2,:)',GRID.x,dE{2}(2,:)',GRID.x,dE{3}(2,:)',GRID.x,dE{4}(2,:)');
        % subplot(339);plot(GRID.x,diff_E{1}(2,:) , GRID.x,diff_E{2}(2,:) , GRID.x,diff_E{3}(2,:));
        subplot(338);plot(TIME,PHASE);
        drawnow
    end

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









