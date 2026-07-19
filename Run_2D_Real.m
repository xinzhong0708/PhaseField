% Clear and restart
clear; figure(3); clf

% Repository root
repo_root = fileparts(mfilename('fullpath'));

% Required code paths
addpath(fullfile(repo_root,'config'), genpath(fullfile(repo_root,'src')), genpath(fullfile(repo_root,'Thermo')))

% Load map or checkpoint
load('Map2d.mat')

% Load current metadata after STATE, GRID, MODEL and PARAM are known
[PHYS,NUM,PARAM,MODEL,GRID] = Read_PFM_Metadata('MetaData.xlsx',GRID,MODEL,STATE,PARAM.eta,PARAM);

%DISPLAY ELEMENT
disp([mean(STATE.E{1},'all') mean(STATE.E{2},'all') mean(STATE.E{3},'all')  ])

% load 1300
% NUM.dt_phy    = NUM.dt_phy/4;
% NUM.dp_target = 0.01;
% NUM.dE_target = 0.01;

for it = 1:2e4

    % SAVE CHECKPOINT
    if mod(it,100)==0
        save(num2str(it))
    end
    
    % DIRECT COUPLED SOLVER WITH MU-MATCHED E RECOVERY
    STATE_OLD            =    STATE;
    t_step               =    tic;

    % UPDATE P-T THERMODYNAMICS
    [Tcur,Pcur]          =    PT_Path(NUM.t_phy,PARAM.PT);
    [MODEL,PARAM]        =    Update_Model_PT(MODEL,PARAM,PHYS,Tcur,Pcur);

    % DAMPED ETA
    p_eta                =    CollapsePForEta(STATE_OLD.p,MODEL.phase_index);
    PARAM.eta            =    Eta_Damping_SmoothHalo(p_eta, PHYS.eta, NUM.int_damp*PHYS.eta, 4, 3e-3, 4, 3e-3, NUM.int_damp*PHYS.eta, 3, NUM.int_damp*PHYS.eta,1e-3);
    
    % W scaling
    PARAM                =    Calc_Kappa_WScale_InterfaceRamp(STATE_OLD,PARAM,MODEL,PHYS);
    
    % FIRST LOCAL EQUILIBRIUM
    t                    =    tic;
    STATE_OLD            =    LE_Run_Mode(STATE_OLD,PARAM,MODEL);
    t_LE1                =    toc(t);

    % BUILD THE MASK FOR ACCH SOLVER
    t                    =    tic;
    PARAM                =    Update_PF_SolverMasks(PARAM,STATE_OLD,MODEL,GRID,PHYS,NUM,'ACCH');

    % UPDATE M and L
    PARAM                =    Compute_M_And_L(STATE_OLD,PARAM,MODEL,PHYS);  

    % Full AC + CH COUPLED PREDICTOR
    [STATE_RAW,DIAG_RAW] =    PF_Coupled_ACCH_LETangent_CS_offdiagM(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);
    t_ACCH               =    toc(t);
    STATE_RAW            =    Snap_PureP(STATE_RAW,PARAM,MODEL);
    
    % TEST ONLY FOR PSEUDO 1D
    % STATE_RAW            =    Horizon_Ave(STATE_RAW,MODEL);

    % EXTEND PHASE
    STATE_RAW            =    Extend_AbsentPhaseC_Rim(STATE_RAW,PARAM);

    % SECOND LOCAL EQUILIBRIUM
    t                    =    tic;
    STATE_LE0            =    LE_Run_Mode(STATE_RAW,PARAM,MODEL);
    t_LE2                =    toc(t);

    % Fixed-p chemical corrector
    t                    =    tic;
    [STATE_TRIAL,DIAG_CHLE] = PF_CH_LECorrector_FixedP_Band_CS_offdiagM(STATE_OLD,STATE_LE0,PARAM,MODEL,GRID,PHYS,NUM);
    t_CHLE               =    toc(t);



    % EXTEND PHASE
    STATE_TRIAL          =    Extend_AbsentPhaseC_Rim(STATE_TRIAL,PARAM);

    % FINAL LOCAL EQUILIBRIUM AFTER FIXED-P CHEMICAL CORRECTOR
    t                    =    tic;
    STATE_TRIAL          =    LE_Run_Mode(STATE_TRIAL,PARAM,MODEL);
    t_LE3                =    toc(t);
    STATE_TRIAL.CHLE_diag=    DIAG_CHLE;
    
    % TIME STEP UPDATE
    dt_try               =    NUM.dt_phy;
    time_old             =    NUM.t_phy;
    [STATE,NUM,DIAG_TS]  =    Update_TimeStep_Soft(STATE,STATE_TRIAL,PARAM,MODEL,NUM);
    
    % PRINT TIME
    t_total              =    toc(t_step);
    disp(['IT: ',num2str(it),'. Total time:',num2str(t_total),' LE1:',num2str(t_LE1),' ACCH:',num2str(t_ACCH),' LE2:',num2str(t_LE2),' LE3:',num2str(t_LE3),' CH:',num2str(t_CHLE),' accept:',num2str(DIAG_TS.accept)])


    %Plotting
    TIME(it)             =    NUM.t_phy;
    DTPHY(it)            =    NUM.dt_phy;
    phase_ids            =    1:numel(MODEL.phs_name);
    PHASE(it,:)          =    zeros(1,numel(phase_ids));
    for iph = 1:numel(phase_ids)
        grains = find(MODEL.phase_index == phase_ids(iph));
        PHASE(it,iph) = mean(sum(STATE.p(:,:,grains),3),'all');
    end

    if mod(it,10)==0
        disp(PHASE(end,:))
        PF_Plot([3,3,1],'E3',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,2],'mu_e1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,3],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,4],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        % subplot(334);plot(GRID.x,STATE.chi{3,3}(2,:),'.-')
        PF_Plot([3,3,5],'omg12',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        % PF_Plot([3,3,6],'c31',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,6],'omg23',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        % subplot(337);plot(GRID.x,STATE.p(2,:,2),'.-',GRID.x,STATE.p(2,:,3),'.-',GRID.x,STATE.p(2,:,4),'.-',GRID.x,STATE.p(2,:,5),'.-',GRID.x,STATE.p(2,:,6),'.-',GRID.x,STATE.p(2,:,7),'.-',GRID.x,STATE.p(2,:,8),'.-')
        subplot(337);plot(GRID.x,STATE.p(2,:,2),'.-',GRID.x,STATE.p(2,:,3),'.-',GRID.x,STATE.p(2,:,4),'.-')
        subplot(338);plot(GRID.x,STATE.E{1}(2,:), GRID.x,STATE.E{2}(2,:)  , GRID.x,STATE.E{3}(2,:) )
        subplot(339);plot(TIME,PHASE(:,:))
        drawnow
    end

    % if mod(it,2)==0
    %     disp(PHASE(end,:))
    %     subplot(231); plot(STATE.E{3}(2,:));title('E')
    %     subplot(232); plot(STATE.mu_e{1}(2,:));title('mu_e')
    %     subplot(233); plot(STATE.p(2,:,3),'.-');title('p')
    %     hold on;plot(STATE.p(2,:,4),'.-');title('p')
    %     hold on;plot(STATE.p(2,:,5),'.-');title('p');hold off
    %     subplot(234); plot(STATE.omg(2,:,3)-STATE.omg(2,:,4),'.-');title('\delta\omega')
    %     subplot(235); plot(STATE.omg(2,:,2)-STATE.omg(2,:,3),'.-');title('\delta\omega')
    %     PF_Plot([2,3,6],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)        
    %     drawnow
    % end

end






