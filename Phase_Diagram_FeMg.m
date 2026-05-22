clear,addpath ./ ./Thermo/Utilities/ ./Thermo/Solutions/ ./EOS
T        = linspace(1000,1000 ,1) + 273.15;
P        = linspace( 2.5, 2.5 ,1) * 1e9;
solmod   = 'solution_models_PFM';
% solmod   = 'solution_models_H18';
Cname    = {'Si'  'Fe'      'Mg'   'Ca' 'Al'    'Na'  'O'};
Nsys     = [ 0    0.130955415951136   0.200127032038209   0.056218857518135   0.190079083285261   0.026587583381404];
Nsys(1)  =  1-sum(Nsys);
Nsys     = [Nsys Nsys(1)*2+Nsys(2)+Nsys(3)+Nsys(4)+Nsys(5)/2*3];
% Nsys     = [Nsys Nsys(1)*2+Nsys(2)+Nsys(3)];
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name      = {'Olivine','Corundum','Orthopyroxene','Clinopyroxene','Garnet','Quartz','Kyanite'};
phs_name      = {'Olivine','Garnet','Feldspar'};
% phs_name      = {'Olivine','Clinopyroxene'};
td            =  init_thermo(phs_name,Cname,solmod);
p             =  props_generate(td);     % generate endmember proportions

% Minimization refinement
[g0,v0]       =  tl_g0(T,P,td);
[g,Npc,pc_id] =  tl_gibbs_energy(T,P,phs_name,td,p,g0,v0);

% Normalize
g             =  g./sum(Npc(1:end-1,:))';
Npc           =  Npc(1:end-1,:)./repmat(sum(Npc(1:end-1,:)),size(Npc,1)-1,1);

%Normalize based on Npc
LB            =  zeros(1,length(g));
UB            =   ones(1,length(g));
[alph,gmin]   =  linprog(real(g),[],[],Npc,Nsys(1:end-1),LB,UB);     
alph          =  alph/sum(alph);
id            =  find(alph>1e-10);
disp('Phase proportion')
disp(num2str(alph(id)))
disp('Phase ID')
disp(num2str(pc_id(id)))
disp('Phase comp')
disp(num2str(Npc(:,id)))

% g'*alph

% pp = [p{1} ; p{2} ; p{3}];
% pp(id,:)

% p{1}(id(1),:)

% ncol  = 5;
% p_glo = cellfun(@(A) [A, nan(size(A,1), ncol-size(A,2))],p, 'UniformOutput', false);
% p_glo = vertcat(p_glo{:});
% 
% p_glo(id,:)