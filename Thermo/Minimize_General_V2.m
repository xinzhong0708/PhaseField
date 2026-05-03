clear;figure(2);clf;addpath([cd,'\glpkmex'])
%Load data
Phase = {'Cpx','Grt','Olv','Opx','Melt'};
for i = 1:length(Phase)
    load(['Data_',Phase{i}])
end

%Bulk composition
C1           = 0.1;   % Fe
C2           = 0.3;   % Mg
C3           = 0.2;   % Al
C4           = 0.1;   % Ca
C5           = 1-C1-C2-C3-C4;

%Pseudocompound
x            = 1e-9:0.020:0.66;
xx           = 1e-9:0.020:0.66;
[x1,x2,x3,x4]= ndgrid(x,x,xx,x);
x1           = x1(:);
x2           = x2(:);
x3           = x3(:);
x3           = x3(:);
x4           = x4(:);
x5           = 1-x1-x2-x3-x4; 

%Function
F{1}         = f_Cpx;
F{2}         = f_Grt;
F{3}         = f_Olv;
F{4}         = f_Opx;
F{5}         = f_Melt;

%Energy
g            = [];
Npc          = [];
phs          = [];
for i = 1:length(Phase)
    %Composition
    g_temp   =  F{i}({x1,x2,x3,x4});
    comp     = [x1';x2';x3';x4';x5'];
    id       =  min(comp)>0;
    g        = [g ; g_temp(id)];
    Npc      = [Npc , comp(:,id)];
    phs      = [phs , i*ones(1,sum(id))];
end
Nsys         = [C1 C2 C3 C4 C5];

%Shorten
id           = find( (imag(g)>1e-12));
g(id)        = [];
Npc(:,id)    = [];
phs(:,id)    = [];
g            = real(g);
id           = find(g>1e9);
g(id)        = [];
Npc(:,id)    = [];
phs(:,id)    = [];

%Remove zeros
id           = find( (imag(g)<1e-10).*isfinite(g)  ==1);
LB           = zeros(1,length(g));
[sol,gmin]   = glpk(g,Npc,Nsys,LB,[],repmat('S',1,size(Npc,1)));     
a            = find(sol>1e-5);

%Display
for i = 1:length(Phase)
    disp(['Phase: ',Phase{i}])
    disp(sum(sol(a(phs(a)==i)))/sum(sol))
    disp(sol(a(phs(a)==i))' * Npc(:,a(phs(a)==i))' / sum( sol(a(phs(a)==i))' * Npc(:,a(phs(a)==i))' ) );
end

disp(gmin)
