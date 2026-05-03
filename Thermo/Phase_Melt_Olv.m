clear,clf,addpath ./ ./Utilities/ ./Solutions/ ./EOS ./bin
run_name = 'Mantle';
T        = 1500 + 273.15;
P        = 1    * 1e9;
solmod   = 'solution_models_H18';
Cname    = {'Si'  ,'Al' , 'Cr',    'Ti'    ,'Fe'   ,'Mn',    'Mg',   'Ca',   'Na',    'K',    'H'   ,'O'};
Nsys     = [ 0.5    0       0        0       0.3     0        0.2      0      0        0       0        ];
Nsys     = [Nsys Nsys(1)*2+Nsys(5)+Nsys(7)];

% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Melt(H18)','Clinopyroxene'};
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
Npc               = [X_Si ; X_Fe ; X_Mg ];
Nsys              = Nsys(1:end-1)/sum(Nsys(1:end-1));
[sol,gmin]        = linprog(g,[],[],Npc,Nsys,LB,[]);


% Normal linprog
% [sol,gmin]        = linprog(g,[],[],Npc(1:end-1,:),Nsys(1:end-1),LB,[]);


%Assign
a                 = find(sol>1e-5);
npc               = Npc(:,a);
factor            = sum(npc);
phs               = pc_id;
% sol(a)            = sol(a).*factor';

%Normalize to 1 mole of cation atom
disp('Melt')
disp(sum(sol(a(phs(a)==1)))/sum(sol))
% % disp(sol(a(phs(a)==1))' * Npc(:,a(phs(a)==1))' / sum( sol(a(phs(a)==1))' * Npc(:,a(phs(a)==1))' ) );

disp('Olv')
disp(sum(sol(a(phs(a)==2)))/sum(sol))
% disp(sol(a(phs(a)==2))' * Npc(:,a(phs(a)==2))' / sum( sol(a(phs(a)==2))' * Npc(:,a(phs(a)==2))' ) );


g1=g;
sol1=sol;
% plot(g1,'.')

plot(g,'o')









%Second one
load Data_Olv
load Data_Melt
c_olv_Si     = Npc(1,phs==2);
c_olv_Fe     = Npc(2,phs==2);
c_olv_Mg     = Npc(3,phs==2);

method       = 'cubic';
a1_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{1},c_olv_Fe,c_olv_Mg,method);
a2_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{2},c_olv_Fe,c_olv_Mg,method);
a3_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{3},c_olv_Fe,c_olv_Mg,method);
a4_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{4},c_olv_Fe,c_olv_Mg,method);
a5_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{5},c_olv_Fe,c_olv_Mg,method);
a6_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{6},c_olv_Fe,c_olv_Mg,method);

c_melt_Si    = Npc(1,phs==1);
c_melt_Fe    = Npc(2,phs==1);
c_melt_Mg    = Npc(3,phs==1);
a1_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{1},c_melt_Fe,c_melt_Mg);
a2_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{2},c_melt_Fe,c_melt_Mg);
a3_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{3},c_melt_Fe,c_melt_Mg);
a4_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{4},c_melt_Fe,c_melt_Mg);
a5_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{5},c_melt_Fe,c_melt_Mg);
a6_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{6},c_melt_Fe,c_melt_Mg);

%Energy
F_olv        = a1_olv +a2_olv .*c_olv_Fe +a3_olv .*c_olv_Mg +a4_olv .*c_olv_Fe.^2 +a5_olv .*c_olv_Fe .*c_olv_Mg +a6_olv .*c_olv_Mg.^2;
F_melt       = a1_melt+a2_melt.*c_melt_Fe+a3_melt.*c_melt_Mg+a4_melt.*c_melt_Fe.^2+a5_melt.*c_melt_Fe.*c_melt_Mg+a6_melt.*c_melt_Mg.^2;


g            = [F_melt(:) ; F_olv(:)];
% g            = g1;

g2           = g;
hold on;plot(g,'.')
b            = isfinite(g);

[sol,gmin]   = linprog(g,[],[],Npc,Nsys,LB,[]);


%Assign
a                 = find(sol>1e-5);
npc               = Npc(:,a);
factor            = sum(npc);
npc               = npc./repmat(factor,size(npc,1),1);
phs               = pc_id;
% sol(a)            = sol(a).*factor';

%Normalize to 1 mole of cation atom
disp('Melt')
disp(sum(sol(a(phs(a)==1)))/sum(sol))
% % disp(sol(a(phs(a)==1))' * Npc(:,a(phs(a)==1))' / sum( sol(a(phs(a)==1))' * Npc(:,a(phs(a)==1))' ) );

disp('Olv')
disp(sum(sol(a(phs(a)==2)))/sum(sol))
% disp(sol(a(phs(a)==2))' * Npc(:,a(phs(a)==2))' / sum( sol(a(phs(a)==2))' * Npc(:,a(phs(a)==2))' ) );


