clear
addpath ./ ./bin ./Thermo ./Thermo/Utilities ./Thermo/Solutions ./Thermo/EOS ./Thermo/glpkmex

T        = 673;
P        = 3e8;
phs_name = {'Clinopyroxene','Orthopyroxene','Garnet','Kyanite','Quartz'};

Nsys     = [ 0.1142    0.1086    0.0527    0.1556];

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