%Clear and restart
clear;figure(1);clf;addpath('bin');addpath('ThermoData');addpath('Thermo');addpath('Maps');addpath('Thermo\Solutions')

%LOAD MAP
load Map2d.mat

%LOAD METADATA
[PHYS,NUM,PARAM,MODEL]   =    Read_PFM_Metadata('Metadata.xlsx',GRID,MODEL,STATE);

%START TIME STEPS
for it = 1:50000
    % SAVE CHECKPOINT
    if mod(it,200)==0
        save(num2str(it))
    end
    
    PARAM.kappa_eff      =    PHYS.kappa*ones(GRID.ny,GRID.nx);

    % DIRECT COUPLED SOLVER WITH MU-MATCHED E RECOVERY
    STATE_OLD            =    STATE;
    t_step               =    tic;

    % UPDATE P-T THERMODYNAMICS
    [Tcur,Pcur]          =    PT_Path(NUM.t_phy,PARAM.PT);
    [MODEL,PARAM]        =    Update_Model_PT(MODEL,PARAM,PHYS,Tcur,Pcur);

    % DAMPED ETA
    PARAM.eta            =    Eta_Damping(STATE_OLD.p,PHYS.eta,NUM.int_damp*PHYS.eta);

    % FIRST LOCAL EQUILIBRIUM
    t                    =    tic;
    STATE_OLD            =    LE_Run_Mode_New(STATE_OLD,PARAM,MODEL);
    STATE_OLD            =    Extend_AbsentPhaseC_Rim(STATE_OLD,PARAM);
    t_LE1                =    toc(t);

    % BUILD THE MASK FOR ACCH SOLVER
    t                    =    tic;
    PARAM                =    Update_PF_SolverMasks(PARAM,STATE_OLD,MODEL,GRID,PHYS,NUM,'ACCH');
    
    % UPDATE M and L
    PARAM                =    Compute_M_And_L(STATE,PARAM,MODEL,PHYS);

    % CALCULATE FACET
    PARAM                =    Calc_AC_Anisotropy_FacetedStiffness(STATE_OLD,PARAM,MODEL,GRID);

    % Full AC + CH COUPLED PREDICTOR
    [STATE_RAW,DIAG_RAW] =    PF_Coupled_ACCH_LETangent_CS_offdiagM(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);
    t_ACCH               =    toc(t);

    % SECOND LOCAL EQUILIBRIUM
    t                    =    tic;
    STATE_LE0            =    LE_Run_Mode_New(STATE_RAW,PARAM,MODEL);
    STATE_LE0            =    Extend_AbsentPhaseC_Rim(STATE_LE0,PARAM);
    t_LE2                =    toc(t);

    % FIXED-p CHEMISTRY CORRECTOR
    t                    =    tic;
    STATE_TRIAL          =    PF_CH_LECorrector_FixedP_Band_CS_offdiagM(STATE_OLD,STATE_LE0,PARAM,MODEL,GRID,PHYS,NUM);
    t_CHLE               =    toc(t);

    % TIME STEP UPDATE
    [STATE,NUM]          =    Update_TimeStep_Soft(STATE,STATE_TRIAL,PARAM,MODEL,NUM);
    
    % PRINT TIME
    t_total              =    toc(t_step);
    disp(['Total time:',num2str(t_total),' LE1:',num2str(t_LE1),' ACCH:',num2str(t_ACCH),' LE2:',num2str(t_LE2),' CH:',num2str(t_CHLE)])


    % PLOTTING
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
        % PF_Plot([3,3,6],'omg23',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,7],'c11',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,8],'c21',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,9],'PhaseStack',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        drawnow
    end

end

