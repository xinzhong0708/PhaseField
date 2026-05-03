clear,clf,addpath ./ ./Utilities/ ./Solutions/ ./EOS ./bin
run_name = 'Mantle';
T        = linspace(2000,2000, 10);
P        = linspace(1.0 ,2.0 , 10) * 1e9;
solmod   = 'solution_models_H18';
dz       = 1/4;
eps_solv = 2;
Cname    = {'Si'  ,'Al' , 'Cr',    'Ti'    ,'Fe
    '   ,'Mn',    'Mg',   'Ca',   'Na',    'K',    'H'   ,'O'};
Oxname   = {'SiO2','Al2O3','Cr2O3','TiO2'  ,'FeO'  ,'MnO',   'MgO',  'CaO',  'Na2O',  'K2O',  'H2O'     };
wtOx     = [30     0      0        0        0       0         30       0       0        0       0       ];
noxy     = [2      3      3        2        1       1         1        1       1        1       1       ];
ncat     = [1      2      2        1        1       1         1        1       2        2       2       ];
molmOx   = [60.084  101.961 151.9904 79.8658  71.844  70.93744  40.304  56.077 61.97894 94.196  18.01528];
NsysOx   = wtOx./molmOx;
Nsys     = NsysOx.*ncat;Nsys = [Nsys Nsys*(noxy./ncat)']; Nsys = Nsys/sum(Nsys);% someone measured rock composition of a pelite

% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Melt(H18)','Clinopyroxene'};
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
p         = props_generate(td);     % generate endmember proportions
refine_id = ones(length(T)*length(P),1);
% Minimization refinement
[T2d,P2d] = ndgrid(T,P);
for iPT = 1:length(T2d(:))
    if refine_id(iPT) == 1
        [alph_all{iPT},Npc_all{iPT},pc_id_ref{iPT},p_ref{iPT},g_min{iPT}] = tl_minimizer(T2d(iPT),P2d(iPT),Nsys,phs_name,p,td);
        disp(iPT/length(T2d(:)))
    end
end
id = find(strcmp(phs_name,'Melt(H18)')==1);
for iPT = 1:length(T2d(:))
    %Save melt
    b         = pc_id_ref{iPT};
    melt(iPT) = sum(alph_all{iPT}(b==id));
end

% Postprocessing
assemblage_id = zeros(length(T)*length(P),length(phs_name));
for iPT = 1:length(T2d(:))
    [alph_all{iPT},Npc_all{iPT},p_ref{iPT},pc_id_ref{iPT}] = cluster_p(alph_all{iPT},Npc_all{iPT},p_ref{iPT},pc_id_ref{iPT},eps_solv,phs_name);
    assemblage_id(iPT,1:length(alph_all{iPT})) = pc_id_ref{iPT};
end
figure(1);clf
tl_psection(T-273.15,P/1e9,Cname,assemblage_id,phs_name,0,[0,0],8);
figure(2);clf
melt2d = reshape(melt,length(T),[])';
save temp
pcolor(T-273.15,P/1e9,melt2d/max(melt2d(:)));colormap jet;shading interp;colorbar

