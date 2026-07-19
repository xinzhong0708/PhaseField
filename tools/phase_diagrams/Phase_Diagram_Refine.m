clear
addpath ./ ./bin ./Thermo ./Thermo/Utilities ./Thermo/Solutions ./Thermo/EOS ./Thermo/glpkmex

T        = 773;
P        = 4e8;
phs_name      = {'Garnet','Orthopyroxene','Kyanite','Quartz','Staurolite'};

Cname         = {'Fe'         'Mg'       'Ca'    'Al'     'Mn'      'H'     'Si'    'O'};
Nsys          = [ 0.0896    0.0753   -0.0000    0.2030    0.0238    0.0183 ];

options.nref       = 50;
options.dz         = 1/4;
options.z_window   = 0.1;
options.dz_fact    = 1.5;
options.ref_fact   = 1.25;
options.active_tol = 1e-9;
options.verbose    = 1;
options.solmod     = 'metapelite_PFM';
options.Cname      = {'Fe','Mg','Ca','Al','Mn','H','Si','O'};


OUT = Equilibrium_Linprog_Refine(T,P,phs_name,Nsys,options.Cname,options);