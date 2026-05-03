clear,figure(1);clf;addpath ./bin ./Utilities/ ./Solutions/ ./EOS
%Energy scaling
E_sc          =  1e5;
%Choose P-T
T             = [1500  1500  ] + 273.15;
P             = [1     1+1e-5] * 1e9;
solmod        = 'solution_models_H18';
Cname         = {'Si'  ,'Al' , 'Cr',    'Ti'    ,'Fe'   ,'Mn',    'Mg',   'Ca',   'Na',    'K',    'H'     ,'O' };
Oxname        = {'SiO2','Al2O3','Cr2O3','TiO2'  ,'FeO'  ,'MnO',   'MgO',  'CaO',  'Na2O',  'K2O',  'H2O'        };
mass          = [28.085, 26.9815, 51.9961, 47.867, 55.845, 54.9380 , 24.305, 40.078,22.9898, 39.0983,1.0079, 15.999];
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name      = {'Orthopyroxene'};

% Calculate free energy
td            =  init_thermo(phs_name,Cname,solmod);
[g0,v0]       =  tl_g0(T(1),P(1),td);
p             =  props_generate(td);     % generate endmember proportions
[g,Npc,pc_id] =  tl_gibbs_energy(T(1),P(1),phs_name,td,p,g0,v0);    

%Find only Fe, Mg, Al
id            = (Npc(3,:)+Npc(4,:)+Npc(6,:)+Npc(9,:)+Npc(10,:)+Npc(11,:))>0;
Npc(:,id==1)  =  [];
g(id)         =  [];

%Normalize to 1 mole of cation
g             =  g ./sum(Npc(1:end-1,:))';

%Volume per mole of cation
v             =  2e-05;

%Energy per volume (J/m3)
g             =  g/v;

%Energy scaling
g             =  g/E_sc;

%Density
Npc           =  Npc./repmat(sum(Npc(1:end-1,:)),12,1);
rho           =  sum(Npc.*repmat(mass',1,length(Npc)))'./v/1e3;

%Concentration
X_Fe          =  Npc(5,:);%./sum(Npc(1:end-1,:),1);
X_Mg          =  Npc(7,:);%./sum(Npc(1:end-1,:),1);
X_Al          =  Npc(2,:);%./sum(Npc(1:end-1,:),1);
X_Si          =  Npc(1,:);%./sum(Npc(1:end-1,:),1);

%Find minimal and sorting the raw data
M             = [X_Fe(:), X_Mg(:), X_Al(:)];
[M,g]         =  Unique_Min(M,g);
X_Fe          =  M(:,1)';
X_Mg          =  M(:,2)';
X_Al          =  M(:,3)';
X_Si          =  1 - sum(M,2)';
g             =  g';

%Add penalty
w             =  1e5;
pp            =  polyfitn([X_Fe(:),X_Mg(:)],X_Al(:),1);
normdir       =  1./[-pp.Coefficients(1:end-1) 1];
r             =  0.1*(rand(size(X_Fe))-0.5);
X_Fe_fit      = [X_Fe , X_Fe + r.*normdir(1)];
X_Mg_fit      = [X_Mg , X_Mg + r.*normdir(1)];
X_Al_fit      = [X_Al , X_Al + r.*normdir(1)];
X_Si_fit      =  1-X_Fe_fit-X_Mg_fit-X_Al_fit;
g_fit         = [g , g+r.^2*w];
scatter3(X_Fe_fit,X_Mg_fit,X_Al_fit,40,g_fit,'.')

%Use fitting function to fit
[f,mu,S]      =  Fit_Energy(g_fit,T,{X_Fe_fit,X_Mg_fit,X_Al_fit,X_Si_fit});

%Final result
f_Opx         =  f;  % energy
S_Opx         =  S;  % suscetibility
mu_Opx        =  mu; % mu
c_Opx         = {'Fe','Mg','Al','Si'};

save('Data_Opx','f_Opx','S_Opx','mu_Opx','c_Opx')

%Checking
A = f({X_Fe,X_Mg,X_Al});
clf;hold on
plot(real(A)-g,'o')
% plot(g,'.')
