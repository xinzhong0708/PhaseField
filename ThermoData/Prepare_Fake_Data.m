clear; clc; clf
addpath('..\bin')
addpath('..\Thermo')
addpath('..\Thermo\Solutions')

% ------------------------------------------------------------
% Fake pressure-temperature and scaling
% ------------------------------------------------------------
T              =  1;          % fake K
P              =  0;          % fake Pa
E_sc           =  1;          % fake energy scale
vref           =  1;          % fake volume scale

% ------------------------------------------------------------
% Fake elements
Cname          =  {'Fe','Mg'};

% ------------------------------------------------------------
% Fake phases
phase_all      =  {'A','B'};
phase_short    =  {'A','B'};

% ------------------------------------------------------------
% Endmember compositions
n              =  [1 0 ; 0 1 ];

% ------------------------------------------------------------
% Ideal mixing setup
% ------------------------------------------------------------
zt             =  eye(2);
mtpl           =  ones(1,2);

% ------------------------------------------------------------
% Alpha and W
alp            =  ones(3,2);
W0             =  zeros(2,2);

% ------------------------------------------------------------
% Mechanical endmember reference energies
g0_A           =  [6; 0]+1;
g0_B           =  [5; 6];

g0_all         =  {g0_A,g0_B};

% ------------------------------------------------------------
% Save fake phase parameter files
% ------------------------------------------------------------
for ip = 1:numel(phase_all)
    if ip == 1
        w              =  zeros(2,2,3);
        w(:,:,1)       =  W0;
        w(1,1,1)       =  1000;
    end
    if ip == 2
        w              =  zeros(2,2,3);
        w(:,:,1)       =  W0;
    end

    phs_name = phase_all{ip};

    pars               =  struct();
    pars.P             =  P;
    pars.T             =  T;
    pars.E_sc          =  E_sc;
    pars.vref          =  vref;

    pars.phase_name    =  phs_name;
    pars.Cname         =  Cname;

    pars.g0            =  g0_all{ip};
    pars.n             =  n;

    pars.mtpl          =  mtpl;
    pars.zt            =  zt;
    pars.alp           =  alp;
    pars.w             =  w;

    % No penalty endmembers for the first simple benchmark
    pars.gN            =  [];
    pars.nN            =  [];
    pars.penalty       =  1e6;

    % Optional names, useful for debugging
    pars.p_name        =  { [phs_name '_Fe']; [phs_name '_Mg'];  };

    par                =  pars;
    save(['Data_',phase_short{ip}],'par')

end


%Plot
clear
c{1}{1}  = linspace(0,0.10,10000);
c{2}{1}  = linspace(0,1   ,10000);
load Data_A
R1       = PhaseThermo(par,c{1});
load Data_B
R2       = PhaseThermo(par,c{2});

plot(c{1}{1},R1.g , c{2}{1},R2.g)
ylim([-2,3])