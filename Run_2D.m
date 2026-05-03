%Clear and restart
clear;figure(1);clf;addpath([cd,'\bin']);addpath([cd,'\Thermo']);addpath([cd,'\Thermo\Solutions'])

%Load map
load Map2d.mat

%Scales (E_sc from map)
t_sc            =  1;                         % Time scale vary depending on D
L_sc            =  1;                         % Length scale, fixed to 1m (due to g=J/m3)

%Physical parameters
l               =  Lx/60/L_sc;                % interface thickness (m)
sigma           =   1.0/E_sc*L_sc^2;          % surface energy (J/m^2)
kappa           =  0e-8/E_sc*L_sc;            % 4th order term, can be set to 0 if no solvus
M0              =  1e-16*E_sc/L_sc^5/t_sc;    % Diffusion mobility

%Dependent variables
m               =  6*sigma/l;
kap             =  3/4*sigma*l;
dceq            =  0.2;
L               =  4*m/3/kap/(4*dceq^2/M0)/1e2;
disp('Time scale of interface')
disp(l^2/L/kap)
disp(1/L/m)
pause(1)

%Time step
dt_phy          =  1e-12/t_sc;
dt_max          =    1e2/t_sc;
dt_min          =  1e-16/t_sc; 
t_tot           =    1e5/t_sc;
dE_target       =  6e-3;
dp_target       =  6e-3;
dmu_target      =  1e5;
time            =  0;

%Load interface 
F.dgdphi        =  @(phi) 2*m*phi.*(phi - 1).^2 + m*phi.^2.*(2.*phi - 2);

%Phase transition L*kappa matrix
LL              =   L*ones(ny,nx);
Lm              =  LL.*m;
LK              =   L*kap*ones(ny,nx);

%Initialize Source germ
omg             =  zeros(ny,nx,Np);
omg_t           =  omg;
s               =  Calc_S_AllenCahn(phi,p,LL,F,omg);
src             =  Calc_S_Diffusion(p,p,e,dt_phy);

%Get the mask for each phase. Ndial is the dilation pixel
ndial           =  nx+1;
mask            =  Calc_Mask(phi,ndial);

%Mobility
M{1}            =  M0*ones(ny,nx);
M{2}            =  M0*ones(ny,nx);
M{3}            =  M0*ones(ny,nx);
M{4}            =  M0*ones(ny,nx);

%Time loop
dt_phy_0        =  dt_phy;
LE_state        =  [];

%Display E
disp('Composition')

%Time iteration
disp([mean(E{1}(:)) mean(E{2}(:)) mean(E{3}(:)) mean(E{4}(:)) ])

load 14200
dE_target       =  3e-3;
dp_target       =  3e-3;

for it = it:1e5
    
    if mod(it,100)==0
        save(num2str(it))
    end

    tic
    
    %Accepted old state
    Eo                                 =  E;
    co                                 =  c;
    mu_eo                              =  mu_e;
    chio                               =  chi;
    po                                 =  p;
    phio                               =  phi;
    omgo                               =  omg;
    dt_try                             =  dt_phy;

    %TRIAL
    %Build S
    s                                  =  Calc_S_AllenCahn(phio, po, LL, F, omgo);

    %Allen Cahn trial
    phi_t                              =  PF_IMEX_Solver2D_AllenCahn(LK,dx,dy,nx,ny,dt_phy,phio,s,mask);
    
    %Clim and normalize phi
    phi_t                              =  Norm_Phi(phi_t);
    p_t                                =  Calc_p(F,phi_t);

    %Variable penalty
    eta_vec                            =  kappa_from_p_smooth_full(p_t, eta, eta*0.1);

    %LE at trial
    [c_t,mu_e_t,chi_t,omg_t,LE_state]  =  LE_Run(pars, p_t, co, Eo, mu_eo, chio, eta_vec, [0.5 1000], [1 1000], LE_state);

    %Diffusion source using trial p and LE-updated c
    e_t                                =  Calc_e(pars,c_t);
    src                                =  Calc_S_Diffusion(p_t, po, e_t, dt_phy);

    %Solve Cahn Hilliard trial
    [E_t,mu_t]                         =  PF_IMEX_Solver2D_Diffusion_MuOnly(M, kappa, dx, dy, nx, ny, Eo, mu_e_t, dt_phy, chi_t, src);

    %Diagnostics on time increment
    dE                                 =  max(abs(cell2mat(E_t)  - cell2mat(Eo)),   [], 'all');
    dmu                                =  max(abs(cell2mat(mu_t) - cell2mat(mu_eo)), [], 'all');
    dp                                 =  max(abs(p_t(:) - po(:)));

    %Timestep for NEXT step only
    scale                              =  min([dE_target  / max(dE,  eps), dp_target  / max(dp,  eps), dmu_target / max(dmu, eps) ]);
    scale                              =  min(max(scale, 0.5), 1.5);
    dt_next                            =  min(max(dt_try * scale, dt_min), dt_max);

    %Accept or Reject
    if dE > dE_target || dp > dp_target || dmu > dmu_target
        %Reject: keep accepted state unchanged
        dt_phy   =  max(dt_try/5, dt_min);
    else
        %Optional final LE correction at trial (p_t, E_t)
        [c_corr,mu_corr,chi_corr,omg_corr,LE_state] =  LE_Run(pars, p_t, c_t, E_t, mu_t, chi_t, eta_vec, [0.5 1000], [1 1000], LE_state);
        %Accept: now overwrite accepted state
        phi      =  phi_t;
        p        =  p_t;
        c        =  c_corr;
        mu_e     =  mu_corr;
        chi      =  chi_corr;
        omg      =  omg_corr;
        E        =  E_t;
        e        =  Calc_e(pars,c);
        time     =  time + dt_try;
        dt_phy   =  dt_next;
    end
    toc
    
    %Plotting
    TIME(it)  = time;
    PHASE1(it)= mean(mean(p(:,:,1)));
    PHASE2(it)= mean(mean(p(:,:,2)));

    %2D
    if mod(it,5)==0
        [~,phase_ID] = max(p,[],3);
        disp([dt_phy])
        disp([mean(mean(p(:,:,1))),mean(mean(p(:,:,2))),mean(E{1}(:))])
        subplot(331);pcolor(E{1});colorbar;shading interp;title('E1')
        subplot(332);pcolor(mu_e{1});colorbar;shading interp;title('mu_e')
        subplot(333);pcolor(phi(:,:,1));colorbar;shading interp;title('p1')
        subplot(334);pcolor(phase_ID);colorbar;shading interp;title('p2')
        subplot(335);pcolor(omg(:,:,1)-omg(:,:,2));colorbar;shading interp;title('domg12')
        subplot(336);pcolor(eta_vec);colorbar;shading interp;title('eta')
        subplot(337);pcolor(c{1}{1});colorbar;shading interp;title('c11')
        subplot(338);plot(x,phi(20,:,1),'.-',y,phi(:,20,2),'.-');title(num2str([mean(p(:,:,1),'all'),mean(p(:,:,2),'all')]))
        subplot(339);plot(TIME,PHASE2,'.-');title('Phase2')
        drawnow
    end

end

