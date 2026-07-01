%Clear and restart
clear;figure(2);clf;addpath('bin');addpath('ThermoData');addpath('Thermo');addpath('Maps');addpath('Thermo\Solutions')

%LOAD MAP OR CHECKPOINT
metadata_file = 'Metadata.xlsx';
map_file      = fullfile('Maps','Map2d.mat');
restart_file  = '';      

if isempty(restart_file)
    load(map_file)
    NUM_RESTART = [];
    it_offset = 0;
else
    load(restart_file)
    NUM_RESTART = NUM;
    [~,restart_base] = fileparts(restart_file);
    it_offset = str2double(restart_base);
    if isnan(it_offset)
        it_offset = 0;
    end
end

%LOAD CURRENT METADATA AFTER THE STATE IS KNOWN
[PHYS,NUM,PARAM,MODEL,GRID] = Read_PFM_Metadata(metadata_file,GRID,MODEL,STATE,PARAM.eta,PARAM);

%If restarting, keep the accepted physical time but refresh metadata-controlled
%thermodynamics, mobility, and timestep limits from metadata_file.
if ~isempty(NUM_RESTART)
    t_restart     = NUM_RESTART.t_phy;
    dt_restart    = NUM_RESTART.dt_phy;
    NUM.t_phy     = t_restart;
    NUM.time      = t_restart;
    NUM.dt_phy    = min(max(dt_restart,NUM.dt_min),NUM.dt_max);
    NUM.dt_good_count = 0;
end

%DISPLAY ELEMENT
disp([mean(STATE.E{1},'all') mean(STATE.E{2},'all') mean(STATE.E{3},'all')  mean(STATE.E{4},'all')])

for it_local = 1:1e4
    it = it_offset + it_local;
    if isfield(NUM,'t_tot') && NUM.t_phy >= NUM.t_tot
        disp('Reached end of metadata P-T path time.')
        break
    end
    if isfield(NUM,'t_tot') && NUM.t_phy + NUM.dt_phy > NUM.t_tot
        NUM.dt_phy = max(NUM.t_tot - NUM.t_phy,NUM.dt_min);
    end
    NUM.istep = it;

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
    p_eta                =    CollapsePForEta(STATE_OLD.p,MODEL.phase_index);
    PARAM.eta            =    Eta_Damping_SmoothHalo(p_eta, PHYS.eta, NUM.int_damp*PHYS.eta, 4, 3e-3, 4, 3e-3, NUM.int_damp*PHYS.eta, 2, NUM.int_damp*PHYS.eta,1e-3);

    %W scaling
    PARAM                =    Calc_WScale(STATE_OLD,PARAM,MODEL,PHYS,NUM);

    % FIRST LOCAL EQUILIBRIUM
    t                    =    tic;
    STATE_OLD            =    LE_Run_Mode(STATE_OLD,PARAM,MODEL);
    % STATE_OLD            =    Extend_AbsentPhaseC_Rim(STATE_OLD,PARAM);
    t_LE1                =    toc(t);

    % BUILD THE MASK FOR ACCH SOLVER
    t                    =    tic;
    PARAM                =    Update_PF_SolverMasks(PARAM,STATE_OLD,MODEL,GRID,PHYS,NUM,'ACCH');
    
    % UPDATE M and L
    PARAM                =    Compute_M_And_L(STATE_OLD,PARAM,MODEL,PHYS);

    % CALCULATE FACET
    % PARAM.aniso_scale_chemical = 0;   % physical capillary mode, do not scale chemical driving force
    % PARAM                =    Calc_Anisotropy_FacetAngular(STATE_OLD,PARAM,GRID,PHYS,MODEL);
    % PARAM                =    Calc_Anisotropy_Facet(STATE_OLD,PARAM,GRID,PHYS);

    % Full AC + CH COUPLED PREDICTOR
    [STATE_RAW,DIAG_RAW] =    PF_Coupled_ACCH_LETangent_CS_offdiagM(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);
    t_ACCH               =    toc(t);

    % SECOND LOCAL EQUILIBRIUM
    t                    =    tic;
    STATE_LE0            =    LE_Run_Mode(STATE_RAW,PARAM,MODEL);
    t_LE2                =    toc(t);

    % Fixed-p chemical corrector
    t                    =    tic;
    [STATE_TRIAL,DIAG_CHLE] = PF_CH_LECorrector_FixedP_Band_CS_offdiagM(STATE_OLD,STATE_LE0,PARAM,MODEL,GRID,PHYS,NUM);
    t_CHLE               =    toc(t);

    % FINAL LOCAL EQUILIBRIUM AFTER FIXED-P CHEMICAL CORRECTOR
    t                    =    tic;
    STATE_TRIAL          =    LE_Run_Mode(STATE_TRIAL,PARAM,MODEL);
    t_LE3                =    toc(t);
    STATE_TRIAL.CHLE_diag =   DIAG_CHLE;

    % TIME STEP UPDATE
    dt_try               =    NUM.dt_phy;
    time_old             =    NUM.time;
    [STATE,NUM,DIAG_TS]  =    Update_TimeStep_Soft(STATE,STATE_TRIAL,PARAM,MODEL,NUM);
    
    % PRINT TIME
    t_total              =    toc(t_step);
    disp(['IT: ',num2str(it),'. Total time:',num2str(t_total),' LE1:',num2str(t_LE1),' ACCH:',num2str(t_ACCH),' LE2:',num2str(t_LE2),' LE3:',num2str(t_LE3),' CH:',num2str(t_CHLE),' accept:',num2str(DIAG_TS.accept)])


    %Plotting
    TIME(it)             =    NUM.time;
    DTPHY(it)            =    NUM.dt_phy;
    phase_ids            =    1:numel(MODEL.phs_name);
    PHASE(it,:)          =    zeros(1,numel(phase_ids));
    for iph = 1:numel(phase_ids)
        grains = find(MODEL.phase_index == phase_ids(iph));
        PHASE(it,iph) = mean(sum(STATE.p(:,:,grains),3),'all');
    end

    if mod(it,20)==0
        disp(PHASE(end,:))
        PF_Plot([3,3,1],'E3',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,2],'mu_e1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,3],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,4],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,5],'omg23',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,6],'c11',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,7],'c21',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,8],'p3',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        subplot(339); plot(TIME,PHASE(:,:))
        drawnow
    end

end
