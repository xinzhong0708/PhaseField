clear;figure(1);clf
%Melt
x_m  = linspace(0.001,0.999,1000);
f_m  = -3*(1-x_m) + x_m.*log(x_m) + (1-x_m).*log(1-x_m);
%Solid
x_s  = 0.6;
f_s  = interp1(x_m,f_m,x_s)-0.2;

Nsys = [0.3 0.7];

plot(x_m,f_m,x_s,f_s,'+')

g            = [f_m(:) ; f_s];
Npc          =[[x_m,x_s];[1-x_m,1-x_s]];

LB           = zeros(1,length(g));
UB           = 1*ones(1,length(g));
sol          = linprog(g,[],[],Npc,Nsys,LB,[]);
a            = find(sol>1e-5);
sol(a) 
plot(x_m,f_m,x_s,f_s,'+',Npc(1,a(1)),g(a(1)),'o',Npc(1,a(2)),g(a(2)),'o')
