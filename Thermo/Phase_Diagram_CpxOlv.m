clear,clf,addpath ./ ./Utilities/ ./Solutions/ ./EOS ./bin
run_name = 'Mantle';
T        = linspace(1500,1501, 2) + 273.15;
P        = linspace(1.0 ,1.01, 2) * 1e9;
solmod   = 'solution_models_H18';
dz       = 1/4;
eps_solv = 2;
Cname    = {'Si'  ,'Al' , 'Cr',    'Ti'    ,'Fe'   ,'Mn',    'Mg',   'Ca',   'Na',    'K',    'H'   ,'O'};
Oxname   = {'SiO2','Al2O3','Cr2O3','TiO2'  ,'FeO'  ,'MnO',   'MgO',  'CaO',  'Na2O',  'K2O',  'H2O'     };
wtOx     = [50     0      0        0        20      0         20       0       0        0       0       ];
noxy     = [2      3      3        2        1       1         1        1       1        1       1       ];
ncat     = [1      2      2        1        1       1         1        1       2        2       2       ];
molmOx   = [60.084  101.961 151.9904 79.8658  71.844  70.93744  40.304  56.077 61.97894 94.196  18.01528];
NsysOx   = wtOx./molmOx;
Nsys     = NsysOx.*ncat;Nsys = [Nsys Nsys*(noxy./ncat)']; Nsys = Nsys/sum(Nsys);% someone measured rock composition of a pelite

% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Melt(H18)','Clinopyroxene','Olivine'};
td       = init_thermo(phs_name,Cname,solmod);
c_exc    = find(Nsys==0);
for i_sol = 1:length(phs_name)
    exc_sol(i_sol) = sum(sum((td(i_sol).n_em(:,c_exc))>0,2)>0)==size(td(i_sol).n_em,1);
end
Cname(Nsys==0) = [];%
Nsys(Nsys==0)  = [];
phs_name(exc_sol) = [];%
td       = init_thermo(phs_name,Cname,solmod);
% for i_sol = 1:length(phs_name)
%     td(i_sol).dz(:) = dz;
% end
options.nref     = 20; % max number of iterations
options.eps_dg   = 1e-12; % tolerance to stop iterations when difference between global gibbs minimimum is below this
options.dz_tol   = 1e-14; % tolerance to stop iterations when z window becomes below this
options.z_window = ones(size(phs_name))*0.001; % the window over which the refined grid is generated
options.dz_fact  = ones(size(phs_name))*3.0;   % the factor to determine new dz spacing, the larger, the more pseudocompounds
options.ref_fact = 8; % the factor to control how the z_window is narrowed each iteration, the larger, the smaller the z window over which new grid is generated
options.disp_ref = 1; % show refinement graphically
options.solver   = 0;
p         = props_generate(td);     % generate endmember proportions
refine_id = ones(length(T)*length(P),1);
% Minimization refinement
[T2d,P2d] = ndgrid(T,P);
for iPT = 1:length(T2d(:))
    if refine_id(iPT) == 1
        [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min{iPT}] = tl_minimizer(T2d(iPT),P2d(iPT),Nsys,phs_name,p,td);
    end
end
id = find(strcmp(phs_name,'Melt(H18)')==1);

% Postprocessing
assemblage_id = zeros(length(T)*length(P),length(phs_name));
for iPT = 1:length(T2d(:))
    [alph_all{iPT},Npc_all{iPT},p_ref{iPT},pc_id_ref{iPT}] = cluster_p(alph_all{iPT},Npc_all{iPT},p_ref{iPT},pc_id_ref{iPT},eps_solv,phs_name);
    assemblage_id(iPT,1:length(alph_all{iPT})) = pc_id_ref{iPT};
end
figure(1);clf
tl_psection(T-273.15,P/1e9,Cname,assemblage_id,phs_name,0,[0,0],8);

%Set values
clearvars -except Npc_all alph_all Nsys
c_melt_Fe = Npc_all{1}(2,1)/sum(Npc_all{1}(1:3,1));
c_melt_Mg = Npc_all{1}(3,1)/sum(Npc_all{1}(1:3,1));
c_cpx_Fe  = Npc_all{1}(2,2)/sum(Npc_all{1}(1:3,2));
c_cpx_Mg  = Npc_all{1}(3,2)/sum(Npc_all{1}(1:3,2));
p_melt    = alph_all{1}(1)/sum(alph_all{1});
p_cpx     = alph_all{1}(2)/sum(alph_all{1});
c         = Nsys(1:3)./sum(Nsys(1:3));
save temp



disp('Cpx')
disp(alph_all{1}(2)/sum(alph_all{1}))
disp(Npc_all{1}(1:end-1,2)'/sum(Npc_all{1}(1:end-1,2)))
disp('Melt')
disp(alph_all{1}(1)/sum(alph_all{1}))
disp(Npc_all{1}(1:end-1,1)'/sum(Npc_all{1}(1:end-1,1)))

