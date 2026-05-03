clear;figure(1);clf
%Load data
Phase = {'Cpx','Olv','Melt'};
for i = 1:length(Phase)
    load(['Data_',Phase{i}])
end

%Bulk composition
C1           = 0.2;
C2           = 0.2;
C3           = 1-C1-C2;

%Pseudocompound
[x1,x2]      = meshgrid(0.05:0.002:0.55 , 0.05:0.002:0.55);
x1           = x1(:);
x2           = x2(:);
x3           = 1-x1-x2; 

%Asign phases
c1_int{1}    = c_Fe_Cpx;
c2_int{1}    = c_Mg_Cpx;
pv_int{1}    = pval_Cpx;
c1_int{2}    = c_Fe_Olv;
c2_int{2}    = c_Mg_Olv;
pv_int{2}    = pval_Olv;
c1_int{3}    = c_Fe_Melt;
c2_int{3}    = c_Mg_Melt;
pv_int{3}    = pval_Melt;

%Interpolate p value
method       ='cubic';
for i = 1:length(c1_int)
    a1{i}    = griddata(c1_int{i},c2_int{i},pv_int{i}{1},x1(:),x2(:),method);
    a2{i}    = griddata(c1_int{i},c2_int{i},pv_int{i}{2},x1(:),x2(:),method);
    a3{i}    = griddata(c1_int{i},c2_int{i},pv_int{i}{3},x1(:),x2(:),method);
    a4{i}    = griddata(c1_int{i},c2_int{i},pv_int{i}{4},x1(:),x2(:),method);
    a5{i}    = griddata(c1_int{i},c2_int{i},pv_int{i}{5},x1(:),x2(:),method);
    a6{i}    = griddata(c1_int{i},c2_int{i},pv_int{i}{6},x1(:),x2(:),method);
end

%Energy
g            = [];
Npc          = [];
phs          = [];
for i = 1:length(c1_int)
    g        = [g ; a1{i} + a2{i}.*x1 + a3{i}.*x2 + a4{i}.*x1.^2 + a5{i}.*x1.*x2 + a6{i}.*x2.^2 ];
    Npc      = [Npc , [x1';x2';x3'] ];
    phs      = [phs , i*ones(1,length(x1))];
end
Nsys         = [C1 C2 C3];

%Remove zeros6
id           = g<10;
g            = g(id);
Npc          = Npc(:,id);
LB           = zeros(1,length(g));
phs          = phs(id);
sol          = linprog(g,[],[],Npc,Nsys,LB,[]);
a            = find(sol>1e-5);


plot(sol)
%Display
for i = 1:length(c1_int)
    disp(['Phase: ',Phase{i}])
    disp(sum(sol(a(phs(a)==i)))/sum(sol))
    disp(sol(a(phs(a)==i))' * Npc(:,a(phs(a)==i))' / sum( sol(a(phs(a)==i))' * Npc(:,a(phs(a)==i))' ) );
end

% disp('Cpx')
% disp(sum(sol(a(phs(a)==1)))/sum(sol))
% disp(sol(a(phs(a)==1))' * Npc(:,a(phs(a)==1))' / sum( sol(a(phs(a)==1))' * Npc(:,a(phs(a)==1))' ) );
% 
% disp('Melt')
% disp(sum(sol(a(phs(a)==2)))/sum(sol))
% disp(sol(a(phs(a)==2))' * Npc(:,a(phs(a)==2))' / sum( sol(a(phs(a)==2))' * Npc(:,a(phs(a)==2))' ) );
% Npc(:,a)