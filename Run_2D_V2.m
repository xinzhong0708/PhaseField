%Clear and restart
clear;figure(1);addpath([cd,'\bin']);addpath([cd,'\Thermo']);addpath([cd,'\Thermo\Solutions'])

%Load map
load Map2d.mat

%PHYSICS
PHYS.E_sc          =  E_sc;
PHYS.t_sc          =  1;                                                   % Time scale
PHYS.L_sc          =  1;                                                   % Length scale
PHYS.l             =  Lx/90/L_sc;                                          % interface thickness (m)
PHYS.sigma         =  0.5/PHYS.E_sc*PHYS.L_sc^2;                           % surface energy (J/m^2)
PHYS.kappa         =  2e-7/(PHYS.E_sc*PHYS.L_sc^2);                        % 4th order term, can be set to 0 if no solvus
PHYS.D_esti        =  1e-12;
PHYS.chi_ref       =  1e-2;
PHYS.M0            =  PHYS.D_esti*PHYS.t_sc/PHYS.L_sc^2*PHYS.chi_ref;
PHYS.m             =  6*PHYS.sigma/PHYS.l;
PHYS.kap           =  3/4*PHYS.sigma*PHYS.l;
PHYS.dceq          =  0.5;
PHYS.L             =  4*PHYS.m/3/PHYS.kap/(PHYS.dceq^2/PHYS.M0)/50;
PHYS.eta           =  eta;

%NUMERICS
NUM.dt_phy         =   1e-3/PHYS.t_sc;
NUM.dt_max         =    1e3/PHYS.t_sc;
NUM.dt_min         =  1e-16/PHYS.t_sc; 
NUM.t_tot          =    1e5/PHYS.t_sc;
NUM.dE_target      =  5e-2;
NUM.dp_target      =  5e-2;
NUM.dmu_target     =  1e5;
NUM.time           =  0;
NUM.dt_good_count  =  0;
NUM.dt_grow_after  =  8;
NUM.dt_grow_fac    =  1.15;
NUM.dt_shrink_fac  =  0.5;
NUM.err_grow       =  0.25;
NUM.phi_mask_cut   =  1e-6;
NUM.phi_mask_thick =  2;
NUM.norm_phi       =  1;
NUM.cut_phi        =  0;
NUM.norm_E         =  1;
NUM.int_damp       =  0.5;

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

%DISPLAY COMPOSITION
disp([mean(STATE.E{1},'all') mean(STATE.E{2},'all') mean(STATE.E{3},'all') mean(STATE.E{end},'all')])

% load temp
% NUM.dE_target      =  4e-2;
% NUM.dp_target      =  4e-2;

for it = 1:1e5
    if mod(it,100)==0
        save(num2str(it))
    end

    tic
    
    %DIRECT COUPLED SOLVER
    %OLD STATE
    STATE_OLD            =    STATE;
    
    % METHOD1
    % Make sure old state is thermodynamically consistent
    PARAM.eta            =    Eta_Damping(STATE_OLD.p,PHYS.eta,NUM.int_damp*PHYS.eta);
    STATE_OLD            =    LE_Run(STATE_OLD,PARAM,MODEL);

    % One monolithic tangent AC-CH-LE step
    [STATE_TRIAL,DIAG]   =    PF_Coupled_ACCH_LETangent(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);

    % Final nonlinear LE correction only once
    STATE_TRIAL          =    LE_Run(STATE_TRIAL,PARAM,MODEL);


    %TIME STEP UPDATE
    [STATE,NUM]          =    Update_TimeStep_Soft(STATE,STATE_TRIAL,PARAM,MODEL,NUM);

    toc
    
    %Plotting
    TIME(it)             =    NUM.time;
    DTPHY(it)            =    NUM.dt_phy;
    PHASE(it,:)          =    squeeze(mean(STATE.p,[1 2]));
    
    % if mod(it,5)==0
    %     disp(NUM.dt_phy)
    %     disp([mean(mean(STATE.p(:,:,1))),mean(mean(STATE.p(:,:,2))),mean(mean(STATE.p(:,:,3))),mean(mean(STATE.p(:,:,4))),mean(mean(STATE.p(:,:,5))),mean(mean(STATE.p(:,:,6))),mean(mean(STATE.p(:,:,7)))])
    %     subplot(331);plot(GRID.x,STATE.E{1}(3,:),GRID.x,STATE.E{2}(3,:),GRID.x,STATE.E{3}(3,:),GRID.x,STATE.E{4}(3,:));title('E1')
    %     subplot(332);plot(STATE.mu_e{1}(3,:));title('mu_e')
    %     subplot(333);plot(DTPHY,'b.');title('dt')
    %     subplot(334);plot(GRID.x,STATE.phi(3,:,1),'.-',GRID.x,STATE.phi(3,:,2),'.-',GRID.x,STATE.phi(3,:,3),'.-',GRID.x,STATE.phi(3,:,4),'.-',GRID.x,STATE.phi(3,:,5),'.-',GRID.x,STATE.phi(3,:,end-1),'.-',GRID.x,STATE.phi(3,:,end),'.-');title('p2')        
    %     % subplot(334);plot(GRID.x,STATE.phi(3,:,1),'.-',GRID.x,STATE.phi(3,:,2),'.-');title('p2')        
    %     subplot(335);plot(GRID.x,STATE.omg(3,:,1)-STATE.omg(3,:,2),GRID.x,STATE.omg(3,:,2)-STATE.omg(3,:,3),'.-');title('domg12')
    %     subplot(336);plot(PARAM.eta(3,:));title('eta')
    %     subplot(337);plot(GRID.x,STATE.c{2}{1}(3,:),GRID.x,STATE.c{2}{2}(3,:),GRID.x,STATE.c{2}{3}(3,:),GRID.x,STATE.c{2}{4}(3,:));title('c11')
    %     subplot(338);plot(GRID.x,STATE.c{3}{1}(3,:),GRID.x,STATE.c{3}{2}(3,:));title('c21')
    %     subplot(339);plot(TIME,PHASE,'.-');title('Phase2')
    %     % subplot(339);plot(TIME,PHASE,'k--');title('Phase2');hold on
    %     drawnow
    % end


    if mod(it,2)==0
        disp(NUM.dt_phy)
        disp([mean(mean(STATE.p(:,:,1))),mean(mean(STATE.p(:,:,2))),mean(mean(STATE.p(:,:,3))),mean(mean(STATE.p(:,:,4))),mean(mean(STATE.p(:,:,5))),mean(mean(STATE.p(:,:,6))),mean(mean(STATE.p(:,:,7)))])
        [~,phase_ID] = max(STATE.p,[],3);
        subplot(331);pcolor(STATE.E{1});colorbar;shading interp;title('E1')
        subplot(332);pcolor(STATE.mu_e{1});colorbar;shading interp;title('mu_e')
        subplot(333);plot(DTPHY,'b.');title('dt')
        subplot(334);pcolor(phase_ID);colorbar;shading interp;title('p2')
        subplot(335);pcolor(STATE.omg(:,:,1)-STATE.omg(:,:,2));colorbar;shading interp;title('domg12')
        subplot(336);pcolor(STATE.p(:,:,2));colorbar;shading interp;title('p2')
        subplot(337);plot(GRID.y,STATE.phi(10,:,2),'.-',GRID.y,STATE.phi(20,:,2),'.-',GRID.y,STATE.phi(30,:,2),'.-',GRID.y,STATE.phi(40,:,2),'.-',GRID.y,STATE.phi(50,:,2),'.-')
        subplot(338);pcolor(STATE.c{2}{1});colorbar;shading interp;title('c21')
        subplot(339);plot(TIME,PHASE,'.-');title('Phase2')
        drawnow
    end


end


