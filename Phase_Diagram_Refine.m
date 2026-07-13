clear
addpath ./ ./bin ./Thermo ./Thermo/Utilities ./Thermo/Solutions ./Thermo/EOS ./Thermo/glpkmex

T        = 573;
P        = 3e8;
phs_name = {'Orthopyroxene','Garnet','Kyanite','Quartz'};

Nsys     = [ 0.091610625000002   0.087311249999999   0.005008125000000   0.210873333333337];

options.nref       = 50;
options.dz         = 1/4;
options.z_window   = 0.1;
options.dz_fact    = 1.5;
options.ref_fact   = 1.25;
options.active_tol = 1e-9;
options.verbose    = 1;
options.solmod     = 'solution_models_PFM';
options.Cname      = {'Fe','Mg','Ca','Al','Si','O'};


OUT = Equilibrium_Linprog_Refine(T,P,phs_name,Nsys,options);