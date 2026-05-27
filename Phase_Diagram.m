clear
addpath ./ ./Thermo/Utilities/ ./Thermo/Solutions/ ./EOS ./glpkmex

T             = linspace(800, 800,1) + 273.15;
P             = linspace(0.5 ,0.5,1) * 1e9;
solmod        = 'solution_models_PFM';

Cname         = {'Fe'      'Mg'      'Ca'      'Al'      'Na'    'Si'    'O'};
Nsys          = [0.168833333333333   0.156416666666666   0.056183333333333   0.176600000000001   0.022400000000000  ];
Nsys          = [Nsys 1-sum(Nsys)];

%Choose phases considered in Gibbs minimization
phs_name      = {'Garnet','Olivine','Clinopyroxene','Feldspar'};
phs_name      = {'Garnet','Clinopyroxene','Feldspar','Olivine'};
td            =  init_thermo(phs_name,Cname,solmod);
p             =  props_generate(td);     % generate endmember proportions

% Minimization refinement
[g0,v0]       =  tl_g0(T,P,td);
[g,Npc,pc_id] =  tl_gibbs_energy(T,P,phs_name,td,p,g0,v0);

%Normalize
g             =  g./sum(Npc(1:end-1,:))';
Npc           =  Npc(1:end-1,:)./sum(Npc(1:end-1,:));

%Normalize based on Npc
LB            =  zeros(1,length(g));
UB            =   ones(1,length(g));
[alph,gmin]   =  linprog(g,[],[],Npc,Nsys,LB,UB);     




alph        = alph/sum(alph);
id          = find(alph > 1e-6);
disp(Npc(:,id))

disp('Active phase proportion')
disp(alph(id))

disp('Active phase ID')
disp(pc_id(id))

disp(Npc(:,id))

disp('Resulting phase proportion and endmember proportion')

for iph = 1:length(phs_name)
    idp      = find(pc_id == iph);
    ph_prop  = sum(alph(idp));
    if ph_prop > 1e-5
        ph_comp = alph(idp).' * p{iph} / ph_prop;
        fprintf('\n%s\n',phs_name{iph})
        fprintf('Phase proportion:\n')
        disp(ph_prop)
        fprintf('Endmember proportion:\n')
        disp(ph_comp)
    end
end