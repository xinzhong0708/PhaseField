clear
addpath ./ ./Thermo/Utilities/ ./Thermo/Solutions/ ./EOS ./glpkmex

T        = linspace(1000,1000,1) + 273.15;
P        = linspace(2.5,2.5,1) * 1e9;
solmod   = 'solution_models_PFM';

Cname    = {'Si' 'Fe' 'Mg' 'Ca' 'Al' 'Na' 'O'};
Nsys     = [0    0.117000000000000   0.173250000000000   0.053450000000000   0.166200000000001   0.028800000000000];
Nsys(1)  =  1-sum(Nsys);
Nsys     = [Nsys Nsys(1)*2 + Nsys(2) + Nsys(3) + Nsys(4) + Nsys(5)/2*3 + Nsys(6)/2];

%Choose phases considered in Gibbs minimization
phs_name = {'Clinopyroxene','Garnet','Feldspar'};

td       = init_thermo(phs_name,Cname,solmod); % initialize thermodynamic data
for ip = 1:length(phs_name),td(ip).nc(:) = 3;end
p        = props_generate(td);                 % generate pseudocompounds
options.eps_dg     = 1e-12;
options.fsolve     = 1;
options.use_pgrid  = 1;
options.show_react = 1;
% One minimization
[alph,Npc,pc_id,p,g_min] = tl_minimizer(T,P,Nsys,phs_name,p,td);



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