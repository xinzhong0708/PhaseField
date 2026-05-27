%Clear and restart
clear;figure(3);clf;addpath('bin');addpath('ThermoData');addpath('Thermo');addpath('Maps');addpath('Thermo\Solutions')

%Load map
load Map2d.mat

%PHYSICS
PHYS.E_sc          =  E_sc;
PHYS.t_sc          =  1;                                                   % Time scale
PHYS.L_sc          =  1;                                                   % Length scale
PHYS.l             =  GRID.dx*4/L_sc;                                      % interface thickness (m)
PHYS.sigma         =  0.02/PHYS.E_sc*PHYS.L_sc^2;                           % surface energy (J/m^2)
PHYS.kappa         =  0e-10/(PHYS.E_sc*PHYS.L_sc^2);                        % 4th order term, can be set to 0 if no solvus
PHYS.D_esti        =  1e-12;
PHYS.chi_ref       =  1e-2;
PHYS.M0            =  PHYS.D_esti*PHYS.t_sc/PHYS.L_sc^2*PHYS.chi_ref;
PHYS.m             =  6*PHYS.sigma/PHYS.l;
PHYS.kap           =  3/4*PHYS.sigma*PHYS.l;
PHYS.dceq          =  0.5;
PHYS.L             =  4*PHYS.m/3/PHYS.kap/(PHYS.dceq^2/PHYS.M0)/100;
PHYS.eta           =  eta;

%NUMERICS
NUM.dt_phy         =   1e-3/PHYS.t_sc;
NUM.dt_max         =    1e1/PHYS.t_sc;
NUM.dt_min         =  1e-16/PHYS.t_sc; 
NUM.t_tot          =    1e5/PHYS.t_sc;
NUM.dE_target      =  2e-2;
NUM.dp_target      =  2e-2;
NUM.dmu_target     =  5e-2;
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
NUM.int_damp       =  0.05;
NUM.kappa_p_cut    =  1e-6;
NUM.use_Jphi       =  1;
NUM.CHLE_p_cut     =  1e-8;
NUM.CHLE_band_thick=  10;
NUM.CHLE_res_rel   =  [];

%GRIDS
GRID.dx            =  GRID.dx;
GRID.dy            =  GRID.dy;
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
% 
% load 1300
% NUM.dE_target      =  1e-2;
% NUM.dp_target      =  1e-2;
% NUM.dmu_target     =  1e-2;

for it = 1:1e5
    if mod(it,100)==0
        save(num2str(it))
    end

    tic


    %DIRECT COUPLED SOLVER WITH UNCONDITIONAL FIXED-P CH-LE CORRECTOR
    STATE_OLD            =    STATE;

    %DAMPED ETA
    PARAM.eta            =    Eta_Damping(STATE_OLD.p,PHYS.eta,NUM.int_damp*PHYS.eta);

    %Old-state local equilibrium
    STATE_OLD            =    LE_Run(STATE_OLD,PARAM,MODEL);

    %Kappa for fourth-order term
    PARAM                =    Calc_Kappa(STATE_OLD,MODEL,PARAM,NUM);
    PARAM.kappa_eff(:)   =    PHYS.kappa;

    % Full AC + CH coupled predictor
    [STATE_RAW,DIAG]     =    PF_Coupled_ACCH_LETangent(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);

    % Nonlinear thermodynamic response of predictor
    STATE_LE0            =    LE_Run(STATE_RAW,PARAM,MODEL);

    % Fixed-p chemical corrector
    STATE_CORR           =    PF_CH_LECorrector_FixedP_Band(STATE_OLD,STATE_LE0,PARAM,GRID,PHYS,NUM);

    % Final nonlinear LE-consistent state
    STATE_TRIAL          =    LE_Run(STATE_CORR,PARAM,MODEL);





    %TIME STEP UPDATE
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

    if mod(it,5)==0
        disp(NUM.dt_phy)
        disp(PHASE(it,:))
        subplot(331);plot(GRID.x,STATE.E{1}(3,:),GRID.x,STATE.E{2}(3,:),GRID.x,STATE.E{3}(3,:),GRID.x,STATE.E{4}(3,:),GRID.x,STATE.E{end}(3,:));title('E1')
        subplot(332);plot(STATE.mu_e{1}(3,:));title('mu_e')
        subplot(333);plot(DTPHY,'b.');title('dt')
        % subplot(334);plot(GRID.x,STATE.phi(3,:,1),'.-',GRID.x,STATE.phi(3,:,2),'.-',GRID.x,STATE.phi(3,:,3),'.-',GRID.x,STATE.phi(3,:,4),'.-',GRID.x,STATE.phi(3,:,5),'.-',GRID.x,STATE.phi(3,:,end-1),'.-',GRID.x,STATE.phi(3,:,end),'.-');title('p2')        
        subplot(334);plot(GRID.x,STATE.phi(3,:,1),'.-',GRID.x,STATE.phi(3,:,2),'.-',GRID.x,STATE.phi(3,:,end-1),'.-',GRID.x,STATE.phi(3,:,end),'.-');title('p2')        
        subplot(335);plot(GRID.x,STATE.omg(3,:,1)-STATE.omg(3,:,2),GRID.x,STATE.omg(3,:,end-1)-STATE.omg(3,:,end),'.-');title('domg12')
        % subplot(336);plot(GRID.x,STATE.omg(3,:,4)-STATE.omg(3,:,5),GRID.x,STATE.omg(3,:,6)-STATE.omg(3,:,7),'.-');title('domg12')
        subplot(337);plot(GRID.x,STATE.c{1}{1}(3,:),GRID.x,STATE.c{1}{end}(3,:));title('c11')
        subplot(338);plot(GRID.x,STATE.c{2}{1}(3,:),GRID.x,STATE.c{2}{end}(3,:));title('c21')
        subplot(339);plot(TIME,PHASE,'.-');title('Phase2')
        % subplot(339);plot(TIME,PHASE,'k--');title('Phase2');hold on
        drawnow
    end


    % if mod(it,2)==0
    %     disp(NUM.dt_phy)
    %     disp(PHASE(it,:))
    %     PF_Plot([3,3,1],'E1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,2],'mu_e1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,3],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,4],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,5],'omg12',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,6],'omg13',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,7],'omg23',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,8],'c31',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     % PF_Plot([3,3,9],'Phase%',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     PF_Plot([3,3,9],'PhaseStack',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    %     drawnow
    % end


end



function PARAM = Calc_Kappa(STATE,MODEL,PARAM,NUM)
%CALC_KAPPA_AUTO
%
% Simple phase-based spatial kappa.
%
% PARAM.kappa_phase:
%   1 x Nthermo vector.
%
% Example:
%   PARAM.kappa_phase = PARAM.kappa * [1 1 1 0 0 0 1];
%
% This gives kappa only to selected thermodynamic phases.
% The returned kappa_eff is ny x nx.

[ny,nx,~]  = size(STATE.p);

phase_index = MODEL.phase_index(:).';
phase_ids   = unique(phase_index,'stable');
Nphase      = numel(phase_ids);

if isfield(NUM,'kappa_p_cut')
    p_cut = NUM.kappa_p_cut;
else
    p_cut = 1e-5;
end

kappa_eff = zeros(ny,nx);

for iph = 1:Nphase

    pid = phase_ids(iph);

    if pid > numel(PARAM.kappa_phase)
        continue
    end

    kp = PARAM.kappa_phase(pid);

    if kp == 0
        continue
    end

    grains = find(phase_index == pid);
    p_sum  = sum(STATE.p(:,:,grains),3);

    % Smooth weighted version
    p_sum(p_sum < p_cut) = 0;

    kappa_eff = kappa_eff + kp .* p_sum;

end
PARAM.kappa_eff = kappa_eff;
end




function STATE_COEF = Blend_Thermo_Coefficients(STATE_OLD,STATE_NEW,beta)
%BLEND_THERMO_COEFFICIENTS
%
% Build updated thermodynamic coefficients for one corrector solve.
% STATE_REF remains STATE_OLD in the coupled solver, so this does not
% advance the physical state twice.
%
% beta = 1.0 : full updated LE coefficients
% beta = 0.5 : midpoint/secant-like thermodynamic coefficients

if beta < 0 || beta > 1
    error('beta must satisfy 0 <= beta <= 1.')
end

STATE_COEF = STATE_OLD;

%Elemental chemical potentials
for ie = 1:length(STATE_OLD.mu_e)

    STATE_COEF.mu_e{ie} = STATE_OLD.mu_e{ie} ...
                        + beta.*(STATE_NEW.mu_e{ie}-STATE_OLD.mu_e{ie});

end

%Susceptibility matrix
for ie = 1:size(STATE_OLD.chi,1)
    for je = 1:size(STATE_OLD.chi,2)

        STATE_COEF.chi{ie,je} = STATE_OLD.chi{ie,je} ...
                               + beta.*(STATE_NEW.chi{ie,je}-STATE_OLD.chi{ie,je});

    end
end

%Phase elemental compositions
for ig = 1:length(STATE_OLD.e)
    for ie = 1:length(STATE_OLD.e{ig})

        STATE_COEF.e{ig}{ie} = STATE_OLD.e{ig}{ie} ...
                              + beta.*(STATE_NEW.e{ig}{ie}-STATE_OLD.e{ig}{ie});

    end
end

%Grand potentials used by the Allen-Cahn source
STATE_COEF.omg = STATE_OLD.omg ...
               + beta.*(STATE_NEW.omg-STATE_OLD.omg);

end


function val = Max_Cell_Diff_Local(A,B)
%Maximum absolute difference between two cell-array fields.

val = 0;

for i = 1:numel(A)
    val = max(val,max(abs(A{i}(:)-B{i}(:))));
end

end


function val = Mu_Kink_1D(mu_e)
%Maximum first spatial jump of elemental chemical potentials along x.

val = 0;

for ie = 1:length(mu_e)
    dmu = diff(mu_e{ie},1,2);
    val = max(val,max(abs(dmu(:))));
end

end


function val = Omega_Kink_1D(omg,phase_index)
%Maximum first spatial jump of pairwise phase grand-potential differences.

phase_id = unique(phase_index,'stable');
Nphase   = length(phase_id);
val      = 0;

for ip = 1:Nphase-1

    ig = find(phase_index == phase_id(ip),1,'first');

    for jp = ip+1:Nphase

        jg   = find(phase_index == phase_id(jp),1,'first');
        domg = omg(:,:,ig)-omg(:,:,jg);
        dd   = diff(domg,1,2);

        val  = max(val,max(abs(dd(:))));

    end

end

end
