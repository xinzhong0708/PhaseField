%Clear and restart
clear;figure(2);clf;addpath('bin');addpath('ThermoData');addpath('Thermo');addpath('Maps');addpath('Thermo\Solutions')

%LOAD MAP OR CHECKPOINT
metadata_file = 'Metadata.xlsx';
map_file      = fullfile('Maps','Map2d.mat');
load(map_file)

%LOAD CURRENT METADATA AFTER THE STATE IS KNOWN
[PHYS,NUM,PARAM,MODEL,GRID] = Read_PFM_Metadata(metadata_file,GRID,MODEL,STATE,PARAM.eta,PARAM);

%DISPLAY ELEMENT
disp([mean(STATE.E{1},'all') mean(STATE.E{2},'all') mean(STATE.E{3},'all')  mean(STATE.E{4},'all')])

% load 800
% NUM.dt_max=2e5;
for it = 1:2e4

    % SAVE CHECKPOINT
    if mod(it,200)==0
        save(num2str(it))
    end
    
    % PARAM.kappa_eff      =    PHYS.kappa*ones(GRID.ny,GRID.nx);

    % DIRECT COUPLED SOLVER WITH MU-MATCHED E RECOVERY
    STATE_OLD            =    STATE;
    t_step               =    tic;

    % UPDATE P-T THERMODYNAMICS
    [Tcur,Pcur]          =    PT_Path(NUM.t_phy,PARAM.PT);
    [MODEL,PARAM]        =    Update_Model_PT(MODEL,PARAM,PHYS,Tcur,Pcur);

    % DAMPED ETA
    p_eta                =    CollapsePForEta(STATE_OLD.p,MODEL.phase_index);
    PARAM.eta            =    Eta_Damping_SmoothHalo(p_eta, PHYS.eta, NUM.int_damp*PHYS.eta, 4, 3e-3, 4, 3e-3, NUM.int_damp*PHYS.eta, 3, NUM.int_damp*PHYS.eta,1e-3);

    %W scaling
    PARAM                =    Calc_Kappa_WScale_InterfaceRamp(STATE_OLD,PARAM,MODEL,PHYS,NUM);
    


    % % FIRST LOCAL EQUILIBRIUM
    % t                    =    tic;
    % PARAM.LE_mode        =   'GP';
    % STATE_OLD            =    LE_Run_Mode(STATE_OLD,PARAM,MODEL);
    % t_LE1                =    toc(t);
    % 
    % % ACCH predictor
    % t                    =    tic;
    % PARAM                =    Update_PF_SolverMasks(PARAM,STATE_OLD,MODEL,GRID,PHYS,NUM,'ACCH');
    % PARAM                =    Compute_M_And_L(STATE_OLD,PARAM,MODEL,PHYS);
    % [STATE_RAW,DIAG_RAW] =    PF_Coupled_ACCH_LETangent_CS_offdiagM(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);
    % t_ACCH               =    toc(t);
    % 
    % % SECOND LOCAL EQUILIBRIUM
    % t                    =    tic;
    % PARAM.LE_mode        =   'GP';
    % STATE_LE0            =    LE_Run_Mode(STATE_RAW,PARAM,MODEL);
    % STATE_TRIAL          =    STATE_LE0;
    % t_LE2                =    toc(t);
    % 
    % t_LE3                =    0;
    % t_CHLE               =    0;








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

    % if mod(it,20)==0
    %     disp(PHASE(end,:))
    %     PF_Plot([3,3,1],'E3',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,2],'mu_e1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,3],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,4],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,5],'omg12',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,6],'c11',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,7],'c21',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,8],'p1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     subplot(339); plot(TIME,PHASE(:,:))
    %     drawnow
    % end

    if mod(it,20)==0
        disp(PHASE(end,:))
        subplot(231); plot(STATE.E{3}(2,:));title('E')
        subplot(232); plot(STATE.mu_e{3}(2,:));title('mu_e')
        subplot(233); plot(STATE.p(2,:,4),'.-');title('p')
        subplot(234); plot(STATE.omg(2,:,3)-STATE.omg(2,:,4),'.-');title('\delta\omega')
        subplot(235); plot(STATE.omg(2,:,2)-STATE.omg(2,:,3),'.-');title('\delta\omega')
        PF_Plot([2,3,6],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)        
        drawnow
    end

end




function kappa_eff = Kappa_InteriorOnly(p_phase,kappa0,PARAM)
%KAPPA_INTERIORONLY Keep kappa inside phase interiors only.
%
% This prevents the E-kappa term from acting directly across a moving
% phase interface, but still keeps spinodal stabilization inside pure phases.

[ny,nx,Np] = size(p_phase);

p0    = 0.90;
p1    = 0.995;
ngrow = 4;

if isfield(PARAM,'kappa_p0'),    p0    = PARAM.kappa_p0;    end
if isfield(PARAM,'kappa_p1'),    p1    = PARAM.kappa_p1;    end
if isfield(PARAM,'kappa_ngrow'), ngrow = PARAM.kappa_ngrow; end

pmax = max(p_phase,[],3);

% Smooth interior weight
w = (pmax - p0)./(p1 - p0);
w = min(max(w,0),1);
w = w.^2.*(3 - 2*w);

% Remove a halo around phase interfaces.
interface = pmax < p0;

mask = interface;
K    = ones(3,3);

for it = 1:ngrow
    mask = conv2(double(mask),K,'same') > 0;
end

w(mask) = 0;

kappa_eff = kappa0 .* w;

end