function [PHYS,NUM,PARAM] = PFM_Defaults(PHYS,NUM,PARAM)

% Scaling and physics
PHYS.E_sc = 1e9;
PHYS.L_sc = 1e-6;
PHYS.t_sc = 1;

% LE active set
PARAM.LE.p_tail = 1e-6;
PARAM.LE.p_full = 2e-4;
PARAM.LE.p_on   = 1e-5;
PARAM.LE.p_off  = 1e-7;
PARAM.LE.Pmax   = 3;

% LE nonlinear solver
PARAM.LE.alpha_LE = [0.6 0.4 0.3 0.2];
PARAM.LE.iter_LE  = [100 100 100 100];
PARAM.LE.c_tol    = 1e-6;
PARAM.LE.max_ls   = 10;
PARAM.LE.lam_c    = 1e-7;
PARAM.LE.dc_cap   = 1e-2;

% Other groups...
end