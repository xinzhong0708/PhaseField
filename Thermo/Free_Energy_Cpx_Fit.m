clear,figure(1);clf;addpath ./bin ./Utilities/ ./Solutions/ ./EOS
%Energy scaling
E_sc          =  5e4;
%Choose P-T
T             = [800   800   ] + 273.15;
P             = [1     1+1e-5] * 1e9;
solmod        = 'solution_models_H18';
Cname         = {'Fe' 'Mg'   'Ca'  'Si'   'O'};
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name      = {'Clinopyroxene'};

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
w             =  1e9;
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
X_Ca          =  M(:,3);
X_Si          =  M(:,4);

%Use fitting function to fit
[f,mu,S,H]    =  Fit_Energy_Poly(g(:),T,{X_Fe(:),X_Mg(:),X_Ca(:),X_Si(:)});

%Checking
c             = {M(:,1) M(:,2) M(:,3) };
plot(real(f(c))-g,'.')

%Final result
f_Cpx        =  f;  % energy
S_Cpx        =  S;  % suscetibility
mu_Cpx       =  mu; % mu
H_Cpx        =  H;  % H
n_Cpx        =  td.n_em(:,1:end-1);
c_Cpx        = {'Fe','Mg','Al','Ca','Si'};

save('Data_Cpx','f_Cpx','S_Cpx','mu_Cpx','c_Cpx','H_Cpx','n_Cpx')

% f({0.1,0.1,0.31})
A=cell2mat(S_Cpx({0.225+0.05,0.13,0.145-0.05}))
det(A)