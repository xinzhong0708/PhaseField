clear;figure(1);clf;addpath([cd,'\bin']);addpath([cd,'\Thermo']);addpath([cd,'\Thermo\Solutions'])

%Pressure Temperature
T              =  950 + 273.15;    % K
P              =  1e9;             % Pa
E_sc           =  1e8;             % J
vref           =  2e-5;            % m3/mol

%Elements
solmod         = 'solution_models_H18';
Cname          = {'Fe' 'Mg' 'Ca' 'Al' 'Si' 'O'};

%Phases
phase_all      = {'Olivine'    ,'Clinopyroxene','Orthopyroxene','Garnet','Corundum','Quartz'};
phase_short    = {'Olv'        ,'Cpx'          ,'Opx'          ,'Grt'   ,'Cor'     ,'Qtz'   };
scale          =  0.1;
for ip = 1:6

    %Phase
    phs_name       =  phase_all(ip);
    td             =  init_thermo(phs_name,Cname,solmod);
    g0             =  cell2mat(tl_g0(T,P,td));
    n              =  td.n_em(:,1:end-1);

    %Exceptions
    if strcmp(phase_short{ip},'Olv')==1
        n(2  ,end)   =  n(2  ,end)+2e-3*scale;
        n(3  ,end)   =  n(3  ,end)+2e-3*scale;
        n(4  ,end)   =  n(4  ,end)-4e-3*scale;
        n(4  ,3  )   =  n(4  ,3  )+1e-3*scale;
        gN           = [];
        nN           = [];
        penalty      =  0;
        rank([n ; nN])
    end
    if strcmp(phase_short{ip},'Cpx')==1
        n(1  ,end)   =  n(1  ,end)+2e-3*scale;
        n(2  ,end)   =  n(2  ,end)+3e-3*scale;
        n(3  ,end)   =  n(3  ,end)-2e-3*scale;
        n(4  ,end)   =  n(4  ,end)+3e-3*scale;
        gN           = [];
        nN           = [];
        penalty      =  0;
        rank([n ; nN])
    end
    if strcmp(phase_short{ip},'Opx')==1
        gN           = [];
        nN           = [];
        penalty      =  0;
        rank([n ; nN])
    end
    if strcmp(phase_short{ip},'Grt')==1
        n(1  ,end)   =  n(1  ,end)+2e-3*scale;
        n(2  ,end)   =  n(2  ,end)+1e-3*scale;
        n(3  ,end)   =  n(3  ,end)+2e-3*scale;
        gN           = [];
        nN           = [];
        penalty      =  0;
        rank([n ; nN])
    end
    if strcmp(phase_short{ip},'Cor')==1
        gN           = [];
        nN           = [];
        penalty      =  0;
        rank([n ; nN])
    end
    if strcmp(phase_short{ip},'Qtz')==1
        gN           = [];
        nN           = [];
        penalty      =  0;
        rank([n ; nN])
    end

    td.n_em(:,1:end-1) =  n;
    pars               =  td;
    pars.n             =  n;
    pars.P             =  P;
    pars.T             =  T;
    pars.g0            =  g0;
    pars.E_sc          =  E_sc;
    pars.vref          =  vref;
    pars.phase_name    =  phs_name;
    pars.gN            =  gN;
    pars.nN            =  nN;
    pars.penalty       =  penalty;

    %Pseudocompound
    pp                 =  props_generate(td);
    c_psc              =  pp{1};
    e_psc              = (c_psc*n)./sum(c_psc*n,2);
    pars.c_psc         =  c_psc;
    pars.e_psc         =  e_psc;
    if isempty(num2cell(pp{1}(:,1:end-1)',2))
        R              =  PhaseThermo(pars,{1});
    else
        R              =  PhaseThermo(pars,num2cell(pp{1}(:,1:end-1)',2)');
    end
    pars.g_psc         =  R.g;

    %Saver
    par                =  pars;
    save(['Data_',phase_short{ip}],'par')
end


