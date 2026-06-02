%Clear and restart
clear;figure(3);clf;addpath('bin');addpath('ThermoData');addpath('Thermo');addpath('Maps');addpath('Thermo\Solutions')

load Map2d.mat
rng(1)
md                 =  GRID.nx/2;
for ii = 1:5
STATE.c{end}{ii}       =  STATE.c{end}{ii}+0.001*(rand(size(STATE.c{end}{ii}))-0.5);
end
STATE.e            =  Calc_e(MODEL.pars,STATE.c);
STATE.E            =  Calc_E_Tot(STATE.e,STATE.p);

% PARAM.eta          =  eta .* ones(GRID.ny,GRID.nx);

% % Use fixed-E LE once to initialize thermodynamic mu from E,c,p
% STATE              =  LE_Run(STATE,PARAM,MODEL);
% 
% % Then, if using GP-like evolution, suppress negative chi if desired
% [STATE,DIAG_CHI_INIT] = Regularize_Chi_Positive_GP(STATE,PARAM);


%PHYSICS
PHYS.E_sc          =  E_sc;
PHYS.t_sc          =  1;                                                   % Time scale
PHYS.L_sc          =  1e-6;                                                % Length scale
PHYS.l             =  2*GRID.dx/PHYS.L_sc;                                 % interface thickness (m)
PHYS.sigma         =  0.3/(PHYS.E_sc*PHYS.L_sc);                           % surface energy (J/m^2)
PHYS.kappa         =  1e-6/(PHYS.E_sc*PHYS.L_sc^2);                        % 4th order term, can be set to 0 if no solvus
PHYS.D_esti        =  1e-12;
PHYS.chi_ref       =  1e-2;
PHYS.M0            =  PHYS.D_esti*PHYS.t_sc/PHYS.L_sc^2*PHYS.chi_ref;
PHYS.m             =  6*PHYS.sigma/PHYS.l;
PHYS.kap           =  3/4*PHYS.sigma*PHYS.l;
PHYS.dceq          =  0.5;
PHYS.L             =  4*PHYS.m/3/PHYS.kap/(PHYS.dceq^2/PHYS.M0)/10;
PHYS.eta           =  eta;

%NUMERICS
NUM.dt_phy         =   2e-3/PHYS.t_sc;
NUM.dt_max         =    1e3/PHYS.t_sc;
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
NUM.int_damp       =  1;
NUM.kappa_p_cut    =  1e-6;
NUM.use_Jphi       =  1;
NUM.CHLE_p_cut     =  1e-8;
NUM.CHLE_band_thick=  50;
NUM.CHLE_res_rel   =  [];

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
PARAM.use_WScale   =  0;                                                   % 1-make excess energy small at interface; 2-do nothing

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

% Use fixed-E LE once to initialize thermodynamic mu from E,c,p
STATE              =  LE_Run_Mode(STATE,PARAM,MODEL);

% load 100

for it = 1:1e6
    if mod(it,100)==0
        save(num2str(it))
    end

    tic

    %DIRECT COUPLED SOLVER WITH UNCONDITIONAL FIXED-P CH-LE CORRECTOR
    STATE_OLD            =    STATE;

    %DAMPED ETA
    PARAM.eta            =    Eta_Damping(STATE_OLD.p,PHYS.eta,NUM.int_damp*PHYS.eta);

    %Old-state local equilibrium
    PARAM.LE_mode        =   'LE';
    STATE_OLD            =    LE_Run_Mode_New(STATE_OLD,PARAM,MODEL);
    
    %Kappa for fourth-order term
    PARAM.kappa_phase_name = 'Cpx';
    PARAM.kappa_pmin       = 0.999;
    PARAM.kappa_buffer     = 2;
    PARAM                  = Calc_Kappa_PhaseCoreRim(STATE_OLD,PARAM,MODEL,PHYS);


    % Full AC + CH coupled predictor
    STATE_RAW            =    PF_Coupled_ACCH_LETangent(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);

    % Nonlinear thermodynamic response of predictor
    PARAM.LE_mode        =   'GP';
    STATE_LE0            =    LE_Run_Mode_New(STATE_RAW,PARAM,MODEL);

    % Fixed-p chemical corrector
    STATE_CORR           =    PF_CH_LECorrector_FixedP_Band(STATE_OLD,STATE_LE0,PARAM,GRID,PHYS,NUM);

    % Final nonlinear LE-consistent state
    PARAM.LE_mode        =   'GP';
    STATE_TRIAL          =    LE_Run_Mode_New(STATE_CORR,PARAM,MODEL);

    % CHECK
    dE{1} = STATE_TRIAL.E{1}-STATE.E{1};
    dE{2} = STATE_TRIAL.E{2}-STATE.E{2};
    dE{3} = STATE_TRIAL.E{3}-STATE.E{3};
    dE{4} = STATE_TRIAL.E{4}-STATE.E{4};
    dE{5} = STATE_TRIAL.E{5}-STATE.E{5};
  
    % TIME STEP UPDATE
    [STATE,NUM]          =    Update_TimeStep_Soft(STATE,STATE_TRIAL,PARAM,MODEL,NUM);

    toc
    
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

    % if mod(it,1)==0
    %     % disp(NUM.dt_phy)
    %     disp(PHASE(it,:))
    %     % disp(mean(STATE.E{1},'all'))
    %     subplot(331);plot(GRID.x,STATE.E{1}(3,:),GRID.x,STATE.E{2}(3,:),GRID.x,STATE.E{3}(3,:),GRID.x,STATE.E{4}(3,:),GRID.x,STATE.E{end}(3,:));title('E1')
    %     subplot(331);plot(GRID.x,STATE.E{1}(3,:));title('E1')
    %     subplot(332);plot(STATE.mu_e{1}(3,:));title('mu_e')
    %     subplot(333);plot(DTPHY,'b.');title('dt')
    %     subplot(336);plot(CHI)
    %     subplot(334);plot(GRID.x,STATE.p(3,:,1),'.-',GRID.x,STATE.p(3,:,2),'.-',GRID.x,STATE.p(3,:,end-1),'.-',GRID.x,STATE.p(3,:,end),'.-');title('p2')        
    %     subplot(335);plot(GRID.x,STATE.omg(3,:,1)-STATE.omg(3,:,2),GRID.x,STATE.omg(3,:,end-1)-STATE.omg(3,:,end),'.-');title('domg12')
    %     subplot(337);plot(GRID.x,STATE.c{1}{1}(3,:),GRID.x,STATE.c{1}{end}(3,:));title('c11')
    %     % subplot(338);plot(GRID.x,STATE.c{2}{1}(3,:),GRID.x,STATE.c{2}{end}(3,:));title('c21')
    %     subplot(338);plot(GRID.x,dE');
    %     subplot(339);plot(GRID.x,dp');
    % 
    %     % subplot(339);plot(TIME,PHASE,'.-');title('Phase2')
    %     drawnow
    % end


    % 
    % if mod(it,2)==0
    %     disp(NUM.dt_phy)
    %     % disp(PHASE(it,:))
    %     PF_Plot([2,2,1],'E1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([2,2,2],'mu_e1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([2,2,3],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     subplot(224);pcolor(CHI);shading interp;colorbar
    %     drawnow
    % end

    if mod(it,2)==0
        disp(NUM.dt_phy)
        disp(PHASE(end,:))
        PF_Plot([3,3,1],'E1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,2],'mu_e1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,3],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,4],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        % PF_Plot([3,3,5],'omg12',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        % PF_Plot([3,3,6],'omg13',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        % PF_Plot([3,3,7],'omg23',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        subplot(337);pcolor(dE{1});colorbar;shading interp
        subplot(338);pcolor(CHI);colorbar;shading interp
        % PF_Plot([3,3,9],'Phase%',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        PF_Plot([3,3,9],'PhaseStack',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
        drawnow
    end

end



function PARAM = Calc_Kappa_PhaseCoreRim(STATE,PARAM,MODEL,PHYS)
%CALC_KAPPA_PHASECORERIM Switch on kappa in selected phase with weak rim.
%
% Rule:
%   selected phase core -> full kappa
%   selected phase rim  -> weak kappa
%   outside phase       -> zero kappa

ny = size(STATE.p,1);
nx = size(STATE.p,2);

%Default parameters
phase_name = 'Cpx';
pmin       = 0.999;
nbuf       = 2;
rim_frac   = 0.05;
npass      = 1;

if isfield(PARAM,'kappa_phase_name')
    phase_name = PARAM.kappa_phase_name;
end

if isfield(PARAM,'kappa_pmin')
    pmin = PARAM.kappa_pmin;
end

if isfield(PARAM,'kappa_buffer')
    nbuf = PARAM.kappa_buffer;
end

if isfield(PARAM,'kappa_rim_frac')
    rim_frac = PARAM.kappa_rim_frac;
end

if isfield(PARAM,'kappa_smooth_npass')
    npass = PARAM.kappa_smooth_npass;
end

%Allow one phase name or several phase names
if ischar(phase_name) || isstring(phase_name)
    phase_name = cellstr(phase_name);
end

%Find selected phase ids
iph_list = [];

for i = 1:numel(phase_name)

    iph = find(strcmpi(MODEL.phs_name,phase_name{i}),1);

    if isempty(iph)
        error('Calc_Kappa_PhaseCoreRim: cannot find phase %s.',phase_name{i})
    end

    iph_list(end+1) = iph;
end

%Collapse selected phase fraction
p_phase = zeros(ny,nx);

for ig = 1:size(STATE.p,3)
    if any(MODEL.phase_index(ig) == iph_list)
        p_phase = p_phase + STATE.p(:,:,ig);
    end
end

%Selected phase mask
mask_phase = p_phase > pmin;

%Eroded core mask
mask_core = ErodeMask(mask_phase,nbuf);

%Rim is selected phase but not core
mask_rim = mask_phase & ~mask_core;

%Build kappa
kappa_eff = zeros(ny,nx);
kappa_eff(mask_core) = PHYS.kappa;
kappa_eff(mask_rim)  = rim_frac*PHYS.kappa;

%Smooth only inside selected phase
if npass > 0
    kappa_eff = SmoothInsideMask(kappa_eff,mask_phase,npass);
end

%Never leak outside selected phase
kappa_eff(~mask_phase) = 0;

PARAM.kappa_eff = kappa_eff;

%Diagnostics
PARAM.kappa_phase_name = phase_name;
PARAM.kappa_phase_id   = iph_list;
PARAM.kappa_phase_p    = p_phase;
PARAM.kappa_phase_mask = mask_phase;
PARAM.kappa_core_mask  = mask_core;
PARAM.kappa_rim_mask   = mask_rim;

end


function mask = ErodeMask(mask,nbuf)

if nbuf <= 0
    return
end

for k = 1:nbuf

    [ny,nx] = size(mask);

    L = mask(:,[1 1:nx-1]);
    R = mask(:,[2:nx nx]);

    if ny > 1
        U = mask([1 1:ny-1],:);
        D = mask([2:ny ny],:);
        mask = mask & L & R & U & D;
    else
        mask = mask & L & R;
    end
end

end


function A = SmoothInsideMask(A,mask,npass)

for ipass = 1:npass

    [ny,nx] = size(A);

    Lm = mask(:,[1 1:nx-1]);
    Rm = mask(:,[2:nx nx]);

    LA = A(:,[1 1:nx-1]);
    RA = A(:,[2:nx nx]);

    num = A;
    den = double(mask);

    num = num + LA.*Lm + RA.*Rm;
    den = den + double(Lm) + double(Rm);

    if ny > 1
        Um = mask([1 1:ny-1],:);
        Dm = mask([2:ny ny],:);

        UA = A([1 1:ny-1],:);
        DA = A([2:ny ny],:);

        num = num + UA.*Um + DA.*Dm;
        den = den + double(Um) + double(Dm);
    end

    A(mask)  = num(mask)./den(mask);
    A(~mask) = 0;
end

end