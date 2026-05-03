clear,figure(1);clf;addpath ./bin ./Utilities/ ./Solutions/ ./EOS
%Energy scaling
E_sc          =  5e4;
%Choose P-T
T             = [1400  1400  ] + 273.15;
P             = [1     1+1e-5] * 1e9;
solmod        = 'solution_models_H18';
Cname         = {'Fe' 'Mg' 'Al' 'Ca' 'Si'  'O'};
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name      = {'Orthopyroxene'};

% Calculate free energy
td            =  init_thermo(phs_name,Cname,solmod);
[g0,v0]       =  tl_g0(T(1),P(1),td);
p             =  props_generate(td);     % generate endmember proportions
[g,Npc,pc_id] =  tl_gibbs_energy(T(1),P(1),phs_name,td,p,g0,v0);    

%Normalize to 1 mole of cation
g             =  g ./sum(Npc(1:end-1,:))';

%Volume per mole of cation
v             =  2e-05;

%Energy per volume (J/m3)
g             =  g/v;

%Energy scaling
g             =  g/E_sc;

%Remove
Npc           =  Npc(1:end-1,:);
Npc           =  Npc./repmat(sum(Npc),size(Npc,1),1);

%Find minimal and sorting the raw data
[M,g]         =  Unique_Min(Npc',g);

%Add penalty
w             =  2e8;
A             =  td.n_em(:,1:end-1);
A             =  A./repmat(sum(A,2),1,size(A,2));
W             =  null(A, 'r');
for icm = 1:size(W,2)
    rr        =  0.1*(rand(size(M,1),1)-0.5);
    dM        =  repmat(W(:,icm)',length(rr),1).*repmat(rr,1,size(W,1));
    M         =  M+dM;
    g         =  g + w*rr.^2;
end

%Reasign
X_Fe          =  M(:,1);
X_Mg          =  M(:,2);
X_Al          =  M(:,3);
X_Ca          =  M(:,4);
X_Si          =  M(:,5);

%Use fitting function to fit
[f,mu,S,H]    =  Fit_Energy_Poly(g(:),T,{X_Fe(:),X_Mg(:),X_Al(:),X_Ca(:),X_Si(:)});

%Checking
c             = {M(:,1) M(:,2) M(:,3) M(:,4)};
plot(real(f(c))-g,'.')

%Final result
f_Opx         =  f;  % energy
S_Opx         =  S;  % suscetibility
mu_Opx        =  mu; % mu
H_Opx         =  H;  % H
n_Opx         =  td.n_em(:,1:end-1);
c_Opx         = {'Fe','Mg','Al','Ca','Si'};

save('Data_Opx','f_Opx','S_Opx','mu_Opx','c_Opx','H_Opx','n_Opx')

