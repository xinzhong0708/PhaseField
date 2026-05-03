clear,clf,addpath ./ ./Utilities/ ./Solutions/ ./EOS ./bin
run_name = 'Mantle';
T        = 1500 + 273.15;
P        = 1    * 1e9;
solmod   = 'solution_models_H18';
Cname    = {'Si'  ,'Al' , 'Cr',    'Ti'    ,'Fe'   ,'Mn',    'Mg',   'Ca',   'Na',    'K',    'H'   ,'O'};
Nsys     = [ 0      0      0        0       0.2    0        0.2     0      0        0       0        ];
Nsys(1)  =  1-sum(Nsys);
Nsys     = [Nsys Nsys(1)*2+Nsys(2)/2*3+Nsys(5)+Nsys(7)];

% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Clinopyroxene','Olivine','Orthopyroxene','Melt(H18)'};
% phs_name = {'Orthopyroxene'};
td       = init_thermo(phs_name,Cname,solmod);
c_exc    = find(Nsys==0);
for i_sol = 1:length(phs_name)
    exc_sol(i_sol) = sum(sum((td(i_sol).n_em(:,c_exc))>0,2)>0)==size(td(i_sol).n_em,1);
end
Cname(Nsys==0)    = [];%
Nsys(Nsys==0)     = [];
phs_name(exc_sol) = [];%
td                = init_thermo(phs_name,Cname,solmod);
p                 = props_generate(td);     % generate endmember proportions
[g0,v0]           = tl_g0(T,P,td);
[g,Npc,pc_id]     = tl_gibbs_energy(T,P,phs_name,td,p,g0,v0);
LB                = zeros(1,size(g,1));

%Convert to concentration of Si Fe Mg  
g                 = g./sum(Npc(1:end-1,:))';
X_Si              = Npc(1,:)./sum(Npc(1:end-1,:));
X_Fe              = Npc(2,:)./sum(Npc(1:end-1,:));
X_Mg              = Npc(3,:)./sum(Npc(1:end-1,:));
X_Ca              = Npc(3,:)./sum(Npc(1:end-1,:));
Npc               = [X_Si ; X_Fe ; X_Mg ];
Nsys              = Nsys(1:end-1)/sum(Nsys(1:end-1));
[sol,gmin]        = linprog(g,[],[],Npc,Nsys,LB,[]);

%Assign
a                 = find(sol>1e-5);
npc               = Npc(:,a);
factor            = sum(npc);
phs               = pc_id;

%Normalize to 1 mole of cation atom
disp('Phase1')
disp(sum(sol(a(phs(a)==1)))/sum(sol))
% % disp(sol(a(phs(a)==1))' * Npc(:,a(phs(a)==1))' / sum( sol(a(phs(a)==1))' * Npc(:,a(phs(a)==1))' ) );

disp('Phase2')
disp(sum(sol(a(phs(a)==2)))/sum(sol))
% disp(sol(a(phs(a)==2))' * Npc(:,a(phs(a)==2))' / sum( sol(a(phs(a)==2))' * Npc(:,a(phs(a)==2))' ) );

disp('Phase3')
disp(sum(sol(a(phs(a)==3)))/sum(sol))

disp('Phase4')
disp(sum(sol(a(phs(a)==4)))/sum(sol))



