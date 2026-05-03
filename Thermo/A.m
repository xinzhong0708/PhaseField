clear;figure(1);clf
%Load data
load Data_Olv
load Data_Melt
%Bulk composition
c_Fe         = 0.3;
c_Mg         = 0.2;
c_Si         = 1-c_Fe-c_Mg;

%Pseudocompound
[x_Fe,x_Mg]  = meshgrid(0:0.005:0.7,0:0.005:0.7);

%Interpolate
x_Si         = 1-x_Fe-x_Mg;
c_olv_Fe     = x_Fe;
c_olv_Mg     = x_Mg;
c_melt_Fe    = x_Fe;
c_melt_Mg    = x_Mg;

method       ='cubic';
a1_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{1},c_olv_Fe,c_olv_Mg,method);
a2_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{2},c_olv_Fe,c_olv_Mg,method);
a3_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{3},c_olv_Fe,c_olv_Mg,method);
a4_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{4},c_olv_Fe,c_olv_Mg,method);
a5_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{5},c_olv_Fe,c_olv_Mg,method);
a6_olv       = griddata(c_Fe_Olv,c_Mg_Olv,pval_Olv{6},c_olv_Fe,c_olv_Mg,method);

a1_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{1},c_melt_Fe,c_melt_Mg,method);
a2_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{2},c_melt_Fe,c_melt_Mg,method);
a3_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{3},c_melt_Fe,c_melt_Mg,method);
a4_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{4},c_melt_Fe,c_melt_Mg,method);
a5_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{5},c_melt_Fe,c_melt_Mg,method);
a6_melt      = griddata(c_Fe_Melt,c_Mg_Melt,pval_Melt{6},c_melt_Fe,c_melt_Mg,method);

%Energy
F_olv        = a1_olv +a2_olv .*c_olv_Fe +a3_olv .*c_olv_Mg +a4_olv .*c_olv_Fe.^2 +a5_olv .*c_olv_Fe .*c_olv_Mg +a6_olv .*c_olv_Mg.^2;
F_melt       = a1_melt+a2_melt.*c_melt_Fe+a3_melt.*c_melt_Mg+a4_melt.*c_melt_Fe.^2+a5_melt.*c_melt_Fe.*c_melt_Mg+a6_melt.*c_melt_Mg.^2;

g            = [F_olv(:) ; F_melt(:)];
Nsys         = [c_Si c_Fe c_Mg];
Npc          = [ [x_Si(:);x_Si(:)]  , [x_Fe(:);x_Fe(:)] , [x_Mg(:);x_Mg(:)]   ]';
phs          = [ones(1,length(F_olv(:))) , 2*ones(1,length(F_olv(:)))];

%Remove zeros
id           = isfinite(g);
g            = g(id);
Npc          = Npc(:,id);
LB           = zeros(1,length(g));
phs          = phs(id);

F_olv = griddata(c_Fe_Olv ,c_Mg_Olv ,f_Olv ,c_olv_Fe ,c_olv_Mg);
F_melt= griddata(c_Fe_Melt,c_Mg_Melt,f_Melt,c_melt_Fe,c_melt_Mg);

sol          = linprog(g,[],[],Npc,Nsys,LB,[]);
a            = find(sol>1e-5);


plot(sol)

disp('olv')
disp(sum(sol(a(phs(a)==1)))/sum(sol))
disp(sol(a(phs(a)==1))' * Npc(:,a(phs(a)==1))' / sum( sol(a(phs(a)==1))' * Npc(:,a(phs(a)==1))' ) );

disp('Melt')
disp(sum(sol(a(phs(a)==2)))/sum(sol))
disp(sol(a(phs(a)==2))' * Npc(:,a(phs(a)==2))' / sum( sol(a(phs(a)==2))' * Npc(:,a(phs(a)==2))' ) );
Npc(:,a)