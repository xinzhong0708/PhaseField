clear;figure(1);clf;addpath('..\bin');addpath('..\Thermo');addpath('..\Thermo\Solutions')

%Pressure Temperature
T              =  700 + 273.15;    % K
P              =  1.0*1e9;         % Pa
E_sc           =  1e9;             % J
vref           =  2e-5;            % m3/mol

%Elements
solmod         = 'solution_models_PFM';
Cname          = {'Fe' 'Mg' 'Ca' 'Al' 'Na' 'Si' 'O'};

%Phases
phase_all      = {'Olivine'    ,'Clinopyroxene','Orthopyroxene','Garnet','Corundum','Quartz','Corundum','Spinel','Kyanite','Andalusite','Sillimanite','Feldspar'};
% phase_all      = {'Olivine'    ,'Clinopyroxene','Orthopyroxene','Garnet','Corundum','Quartz','Corundum','Spinel','Kyanite','Andalusite','Sillimanite','Feldspar_Pen'};
phase_short    = {'Olv'        ,'Cpx'          ,'Opx'          ,'Grt'   ,'Cor'     ,'Qtz'   ,'Cor'     ,'Spl'   ,'Kya'    ,'And'       ,'Sil'        ,'Fel'     };

%Noise
eps            =  1e-4;
for ip = 1:length(phase_all)

    %Phase
    phs_name       =  phase_all(ip);
    td             =  init_thermo(phs_name,Cname,solmod);
    g0             =  cell2mat(tl_g0(T,P,td));
    n              =  td.n_em(:,1:end-1);
    if size(n)>1
    n              =  Perturb_n_LE(n,eps);
    end
    disp(rank(n))
    
    %Thermodynamic values
    td.n_em(:,1:end-1) =  n;
    pars               =  td;
    pars.n             =  n;
    pars.P             =  P;
    pars.T             =  T;
    pars.g0            =  g0;
    pars.E_sc          =  E_sc;
    pars.vref          =  vref;
    pars.phase_name    =  phs_name;

    %Saver
    par                =  pars;
    save(['Data_',phase_short{ip}],'par')
end


