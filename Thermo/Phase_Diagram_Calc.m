clear,clf,addpath ./ ./Utilities/ ./Solutions/ ./EOS
T        = linspace( 800, 800 ,1) + 273.15;
P        = linspace( 1.0, 1.0 ,1) * 1e9;
solmod   = 'solution_models_H18';
Cname    = {'Si'  'Fe'      'Mg'      'Al'      'Ca'     'O'};
Nsys     = [ 2     1         2        0           1        ];
% Nsys(1)  =  1-sum(Nsys);
Nsys     = [Nsys Nsys(1)*2+Nsys(2)+Nsys(3)+Nsys(4)*1.5+Nsys(5)];

% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name      = {'Melt(H18)'};
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
[alph,gmin]   =  linprog(real(g),[],[],Npc,Nsys(1:end-1),LB,[]);     
alph          =  alph/sum(alph);
id            =  find(alph>1e-10);
disp('Phase proportion')
disp(num2str(alph(id)))
disp('Phase ID')
disp(num2str(pc_id(id)))
disp('Phase comp')
disp(num2str(Npc(:,id)))


sum(g.*alph)

