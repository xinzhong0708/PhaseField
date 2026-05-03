clear,addpath ../ ./Solutions/ ./EOS
T     = 1500 + 273.15;
P     = 1e9;
phase = {'Orthopyroxene'};
Cname         = {'Si'  ,'Al' , 'Cr',    'Ti'    ,'Fe'   ,'Mn',    'Mg',   'Ca',   'Na',    'K',    'H'     ,'O' };
td    = init_thermo(phase,Cname,'solution_models_H18');
p     = props_generate(td);
[g,Npc,pc_id] = tl_gibbs_energy(T,P,phase,td,p);

%Find only Fe, Mg, Al
id            = (Npc(3,:)+Npc(4,:)+Npc(6,:)+Npc(8,:)+Npc(9,:)+Npc(10,:)+Npc(11,:))>0;
Npc(:,id==1)  =  []; 
g(id)         =  [];

%Normalize to cation
g             =  g./sum(Npc(1:end-1,:))';

Npc           =  Npc(1:end-1,:);

X_Fe          =  Npc(3,:)./sum(Npc);
X_Mg          =  Npc(4,:)./sum(Npc);
X_Al          =  Npc(2,:)./sum(Npc);
X_Si          =  Npc(1,:)./sum(Npc);

load Data_Opx
A = real(f_Opx({X_Fe,X_Mg,X_Al}));
clf
plot(g-A'*2,'.')
% hold on
% plot(A*2,'o')

