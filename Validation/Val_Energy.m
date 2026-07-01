clear;figure(1);clf;addpath('../bin')
addpath ../Thermo/Utilities/ ../Thermo/Solutions/

%% General
% Pressure-temperature initial condition
T              = 500 + 273.15;       % K
P              = 0.5e9;              % Pa
E_sc           = 1;
vref           = 1;
Cname          = {'Fe' 'Mg' 'Ca' 'Al' 'Na' 'Si' 'O'};
solmod         = 'solution_models_PFM';

%% Thermolab
%Choose phases considered in Gibbs minimization
phs_name      = {'Melt(H18)'};
td            =  init_thermo(phs_name,Cname,solmod);
p             =  props_generate(td);     % generate endmember proportions

% Minimization refinement
[g0,v0]       =  tl_g0(T,P,td);
[g,Npc,pc_id] =  tl_gibbs_energy(T,P,phs_name,td,p,g0,v0);

%Normalize
g             =  g./sum(Npc(1:end-1,:))';



%% PF code
c_value{1}     = num2cell(p{1}(:,1:end-1)',2)'; 

pars           = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref);
tic
R              = PhaseThermo(pars{1},c_value{1});
toc
plot(g - R.g,'ro')



function pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref)

%BUILD_PARS_PHASES Build pars for present phases at given T and P.
%
% Full-name-only version. No Data_*.mat files are loaded.
%
% Recomputes:
%   td = init_thermo(...)
%   g0 = tl_g0(T,P,td)
%   n  = td.n_em(:,1:end-1)

Nphase     = numel(phs_name);
pars_phase = cell(1,Nphase);

for ip = 1:Nphase

    phase_name = phs_name(ip);

    td = init_thermo(phase_name,Cname,solmod);
    g0 = cell2mat(tl_g0(T,P,td));
    n  = td.n_em(:,1:end-1);

    td.n_em(:,1:end-1) = n;

    pars            = td;
    pars.n          = n;
    pars.P          = P;
    pars.T          = T;
    pars.g0         = g0;
    pars.E_sc       = E_sc;
    pars.vref       = vref;
    pars.phase_name = phase_name;

    pars_phase{ip}  = pars;

end

end