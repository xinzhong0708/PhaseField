%% Clear and restart
clear; figure(3); clf
addpath('..\bin')
addpath('..\ThermoData')

%% ------------------------------------------------------------------------
%  Scaling / basic physical constants
% -------------------------------------------------------------------------
PHYS             = struct();
PHYS.E_sc        = 1e9;
PHYS.L_sc        = 1;

E_sc             = PHYS.E_sc;
L_sc             = PHYS.L_sc;
eta0             = 1000e10/E_sc;      %Must be identical in Run_2D
PHYS.eta         = eta0;

%% ------------------------------------------------------------------------
%  Thermodynamic phases and Gibbs-minimization input
%
%  Use full-precision phase proportions and compositions from
%  Phase_Diagram.m. This is the initial guess and defines E_target.
% -------------------------------------------------------------------------
phs_name         = {'Grt','Olv','Cpx'};
pars_phase       = Load_Data(phs_name);
Nphase           = length(pars_phase);

phase_prop_GM    = [0.30 0.40 0.30];
phase_prop_GM    = phase_prop_GM/sum(phase_prop_GM);

c_value          = cell(1,Nphase);
c_value{1}       = [0.38 0.51];                  %Grt
c_value{2}       = [0.02 0.48 0.34];             %Olv
c_value{3}       = [0.03 0.34 0.12 0.39];       %Cpx

c_guess          = cell(1,Nphase);

for ip = 1:Nphase
    Nc           = Count_Independent_Endmembers(pars_phase{ip});
    c_guess{ip}  = num2cell(c_value{ip}(:).');
end

%% ------------------------------------------------------------------------
%  Target conserved bulk composition
%
%  E_target is the PF composition represented by the Gibbs-minimization
%  phase amounts and compositions. It remains fixed during the refinement.
% -------------------------------------------------------------------------
[~,E_target]     = Evaluate_Mixture(pars_phase,c_guess,phase_prop_GM);
Ne               = length(E_target);

%% ------------------------------------------------------------------------
%  PF penalty-equilibrium refinement
%
%  IMPORTANT:
%  LE_Run does not solve exact Gibbs minimization. It minimizes a finite
%  penalty model with E-Emix = mu/eta. Therefore the initial state must be
%  refined using the SAME penalty model if initial mu_e and delta omega are
%  required to be flat while preserving E_target.
% -------------------------------------------------------------------------
OPT_INIT              = struct();
OPT_INIT.Display      = 'iter';
OPT_INIT.p_min        = 1e-10;
OPT_INIT.c_lower      = -0.10;
OPT_INIT.c_upper      =  1.10;
OPT_INIT.MaxIter      = 1000;
OPT_INIT.MaxFunEval   = 5e4;
OPT_INIT.OptimalityTol= 1e-10;
OPT_INIT.StepTol      = 1e-12;
OPT_INIT.FunctionTol  = 1e-12;

[phase_prop,c_ref,mu_ref,DIAG_INIT] = Refine_Penalty_Equilibrium_PF(pars_phase,phase_prop_GM,c_guess,E_target,eta0,OPT_INIT);

fprintf('\nGibbs-minimization input and PF penalty-refined phase proportions:\n')
fprintf('Phase           Gibbs input       PF refined       difference\n')

for ip = 1:Nphase
    fprintf('%-8s        %.10f      %.10f     %+.3e\n', ...
            phs_name{ip},phase_prop_GM(ip),phase_prop(ip), ...
            phase_prop(ip)-phase_prop_GM(ip))
end

fprintf('\nPF penalty-refined independent endmember proportions:\n')

for ip = 1:Nphase
    fprintf('%s\n',phs_name{ip})
    disp(cell2mat(c_ref{ip}(:)).')
end

%% ------------------------------------------------------------------------
%  Grid
% -------------------------------------------------------------------------
Lx               = 5e-6;
Ly               = 5e-6;

nx               = 300;
ny               = 4;

x                = linspace(0,Lx,nx)/L_sc;
y                = linspace(0,Ly,ny)/L_sc;

dx               = x(2)-x(1);
dy               = y(2)-y(1);

GRID             = struct();
GRID.x           = x;
GRID.y           = y;
GRID.dx          = dx;
GRID.dy          = dy;
GRID.nx          = nx;
GRID.ny          = ny;
GRID.Lx          = Lx/L_sc;
GRID.Ly          = Ly/L_sc;

%% ------------------------------------------------------------------------
%  Spatially uniform refined compositions for every phase
%
%  Absent phases are initialized with the same thermodynamically refined
%  c as present phases. This avoids a poor absent-phase c producing a
%  large omega spike when that phase becomes active.
% -------------------------------------------------------------------------
c_phase          = cell(1,Nphase);

for ip = 1:Nphase

    c_phase{ip} = cell(size(c_ref{ip}));

    for ic = 1:length(c_ref{ip})
        c_phase{ip}{ic} = c_ref{ip}{ic}*ones(ny,nx);
    end

end

%% ------------------------------------------------------------------------
%  Pseudo-2D repeating 1D phase bands
%
%  Nrepeat = 1: Grt | Olv | Cpx
%  Nrepeat = 2: Grt | Olv | Cpx | Grt | Olv | Cpx
% -------------------------------------------------------------------------
Nrepeat          = 1;
band_order       = repmat(1:Nphase,1,Nrepeat);

Noccur           = accumarray(band_order(:),1,[Nphase,1]).';
band_fraction    = phase_prop(band_order)./Noccur(band_order);

grain_phase      = band_order;
Ngrain           = length(grain_phase);
Np               = Ngrain;

x_edges          = 1 + round(nx*[0 cumsum(band_fraction)]);
x_edges(end)     = nx + 1;

phi              = zeros(ny,nx,Ngrain);
phase_ID         = zeros(ny,nx);
grain_ID         = zeros(ny,nx);
seed_xy          = zeros(Ngrain,2);

for ig = 1:Ngrain

    iph                      = grain_phase(ig);
    ix1                      = x_edges(ig);
    ix2                      = x_edges(ig+1)-1;

    phi(:,ix1:ix2,ig)       = 1;
    phase_ID(:,ix1:ix2)     = iph;
    grain_ID(:,ix1:ix2)     = ig;

    seed_xy(ig,1)           = mean(x(ix1:ix2));
    seed_xy(ig,2)           = mean(y);

end

grain_phase      = grain_phase(:);

%% ------------------------------------------------------------------------
%  Grain-resolved MODEL
% -------------------------------------------------------------------------
MODEL                    = struct();
MODEL.phs_name           = phs_name;
MODEL.phase_index        = grain_phase(:).';
MODEL.pars               = pars_phase(MODEL.phase_index);

eps_phi                  = 1e-14;
MODEL.p_fun              = @(a,phi)    phi(:,:,a).^2 ./ (sum(phi.^2,3) + eps_phi);
MODEL.dpdphi             = @(a,b,phi) (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) - 2*phi(:,:,a).*phi(:,:,b).^2 ./ (sum(phi.^2,3) + eps_phi).^2;

pars                     = MODEL.pars;
F                        = MODEL;

figure(3); clf
pcolor(x*1e6,y*1e6,phase_ID)
shading flat
axis equal tight
colorbar
title('Initial repeating 1D thermodynamic phase bands')
xlabel('x \mum')
ylabel('y \mum')
drawnow

%% ------------------------------------------------------------------------
%  Construct penalty-consistent local fields
%
%  At the refined penalty equilibrium:
%       E_local = Emix_local + mu_ref/eta0
%
%  The spatial average equals E_target when the grid-realized phase
%  fractions equal phase_prop. The diagnostic below reports rounding error.
% -------------------------------------------------------------------------
c                = Expand_c_By_Phase(c_phase,MODEL.phase_index);
p                = Calc_p(MODEL,phi);
e                = Calc_e(pars,c);
eta              = eta0*ones(ny,nx);

E                = Calc_E_Tot(e,p);

for ie = 1:Ne
    E{ie} = E{ie} + mu_ref{ie}/eta0;
end

mu_e             = cell(1,Ne);

for ie = 1:Ne
    mu_e{ie} = mu_ref{ie}*ones(ny,nx);
end

chi              = repmat({zeros(ny,nx)},Ne,Ne);
omg              = Calc_Omega_Grain(pars,c,e,mu_e,ny,nx);

%% ------------------------------------------------------------------------
%  PARAM and STATE before LE_Run
% -------------------------------------------------------------------------
PARAM            = struct();
PARAM.Np         = Np;
PARAM.Ne         = Ne;
PARAM.eta        = eta;
PARAM.use_WScale = 0;

STATE            = struct();
STATE.c          = c;
STATE.e          = e;
STATE.E          = E;
STATE.mu_e       = mu_e;
STATE.chi        = chi;
STATE.omg        = omg;
STATE.phi        = phi;
STATE.p          = p;
STATE.mask       = ones(ny,nx,Np);
STATE.LE_state   = [];

STATE_INI        = STATE;
omega_ini        = Extract_Omega_By_Phase(STATE_INI.omg,MODEL.phase_index,Nphase);

%% ------------------------------------------------------------------------
%  First LE_Run check
% -------------------------------------------------------------------------
STATE            = LE_Run(STATE,PARAM,MODEL);
omega_le         = Extract_Omega_By_Phase(STATE.omg,MODEL.phase_index,Nphase);

p_phase          = Collapse_p_By_Phase(STATE.p,MODEL.phase_index,Nphase);
phase_prop_map   = zeros(1,Nphase);

for ip = 1:Nphase
    tmp                = p_phase(:,:,ip);
    phase_prop_map(ip) = mean(tmp(:));
end

E_mean           = zeros(1,Ne);
E_target_vec     = cell2mat(E_target(:)).';

for ie = 1:Ne
    E_mean(ie) = mean(STATE.E{ie}(:));
end

dc_max = 0;

for ig = 1:Ngrain
    for ic = 1:length(STATE.c{ig})
        dc_max = max(dc_max,max(abs(STATE.c{ig}{ic}(:)-STATE_INI.c{ig}{ic}(:))));
    end
end

fprintf('\nRefined and grid-realized phase proportions:\n')
fprintf('Phase           refined          map             difference\n')

for ip = 1:Nphase
    fprintf('%-8s        %.10f      %.10f     %+.3e\n', ...
            phs_name{ip},phase_prop(ip),phase_prop_map(ip), ...
            phase_prop_map(ip)-phase_prop(ip))
end

fprintf('\nGrand-potential mismatch before LE_Run:\n')
Print_Omega_Difference(omega_ini,phs_name)

fprintf('\nGrand-potential mismatch after first LE_Run:\n')
Print_Omega_Difference(omega_le,phs_name)

fprintf('max|c after first LE - c before LE| = %.8e\n',dc_max)
fprintf('max|mean(E map) - E_target|         = %.8e\n',max(abs(E_mean-E_target_vec)))

%% ------------------------------------------------------------------------
%  Plot initial grand-potential differences
% -------------------------------------------------------------------------
figure(4); clf

subplot(2,1,1)
Plot_Omega_Difference(x,omega_ini,phs_name)
title('Initial \Delta\omega before LE\_Run')

subplot(2,1,2)
Plot_Omega_Difference(x,omega_le,phs_name)
title('Initial \Delta\omega after first LE\_Run')

%% ------------------------------------------------------------------------
%  Save structured + legacy variables
% -------------------------------------------------------------------------
c                = STATE.c;
e                = STATE.e;
E                = STATE.E;
p                = STATE.p;
mu_e             = STATE.mu_e;
chi              = STATE.chi;

save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'E_sc','L_sc','eta','pars','Np','Ne','Nphase','Ngrain', ...
    'phs_name','phase_prop_GM','phase_prop','phase_prop_map', ...
    'c_value','c_ref','mu_ref','E_target','DIAG_INIT', ...
    'phi','p','c','e','E','mu_e','chi', ...
    'phase_ID','grain_ID','grain_phase','seed_xy','Nrepeat','band_order')

%% ========================================================================
%  Local helper functions
% ========================================================================

function Nc = Count_Independent_Endmembers(par)

Nall = numel(par.g0);

if isfield(par,'gN') && ~isempty(par.gN)
    Nall = Nall + numel(par.gN);
end

Nc = Nall-1;

end


function [phase_prop,c_ref,mu_ref,DIAG] = Refine_Penalty_Equilibrium_PF(pars,phase_prop,c_guess,E_target,eta,OPT)
%REFINE_PENALTY_EQUILIBRIUM_PF Refine the initial state for current LE_Run.
%
% Minimize the SAME finite-eta thermodynamic objective represented by
% LE_Calculator, but allow both compositions and phase proportions to
% adjust:
%
%   F = sum_ip p_ip*g_ip(c_ip)
%       + eta/2 * |E_target - sum_ip p_ip*e_ip(c_ip)|^2
%
% subject to:
%   sum(p_ip) = 1,   p_ip >= p_min.
%
% The optimum satisfies:
%   mu_ref     = eta*(E_target - Emix)
%   omega_ip   = g_ip - e_ip'*mu_ref = common value for present phases
%
% This avoids the ill-conditioned single-phase mu_e initialization that
% caused the previous lsqnonlin refinement to stall.

Np       = length(pars);
Nc       = cellfun(@length,c_guess);
Nct      = sum(Nc);

z0       = zeros(Nct+Np,1);
i0       = 0;

for ip = 1:Np
    ids     = i0+(1:Nc(ip));
    z0(ids) = cell2mat(c_guess{ip}(:));
    i0      = i0+Nc(ip);
end

z0(Nct+(1:Np)) = phase_prop(:)/sum(phase_prop);

lb       = [OPT.c_lower*ones(Nct,1); OPT.p_min*ones(Np,1)];
ub       = [OPT.c_upper*ones(Nct,1); ones(Np,1)];

Aeq      = [zeros(1,Nct),ones(1,Np)];
beq      = 1;

options  = optimoptions('fmincon', ...
           'Algorithm','sqp', ...
           'Display',OPT.Display, ...
           'MaxIterations',OPT.MaxIter, ...
           'MaxFunctionEvaluations',OPT.MaxFunEval, ...
           'OptimalityTolerance',OPT.OptimalityTol, ...
           'StepTolerance',OPT.StepTol, ...
           'FunctionTolerance',OPT.FunctionTol);

fun      = @(z) Penalty_Objective(z,pars,Nc,E_target,eta);

[z,Fmin,exitflag,output] = fmincon(fun,z0,[],[],Aeq,beq,lb,ub,[],options);

if exitflag <= 0
    error('PF penalty-equilibrium refinement did not converge. fmincon exitflag = %d.',exitflag)
end

[c_ref,phase_prop] = Unpack_Penalty_State(z,Nc,Np);
phase_prop         = phase_prop(:).';

[Fchk,PART]        = Penalty_Objective(z,pars,Nc,E_target,eta);
mu_ref             = num2cell(PART.mu(:).');

active             = phase_prop > max(OPT.p_min*10,1e-8);
omega_active       = PART.omega(active);
omega_res          = omega_active-omega_active(1);

DIAG.Fmin              = Fmin;
DIAG.Fcheck            = Fchk;
DIAG.exitflag          = exitflag;
DIAG.output            = output;
DIAG.E_mix             = PART.Emix;
DIAG.E_residual        = PART.Etarget-PART.Emix;
DIAG.mu                = PART.mu;
DIAG.mu_stationarity   = PART.mu_stationarity;
DIAG.omega             = PART.omega;
DIAG.omega_difference  = omega_res;
DIAG.phase_sum_error   = sum(phase_prop)-1;

fprintf('\nPF finite-eta equilibrium refinement:\n')
fprintf('exitflag                      = %d\n',exitflag)
fprintf('Fmin                          = %.12e\n',Fmin)
fprintf('max |Etarget-Emix|            = %.8e\n',max(abs(DIAG.E_residual)))
fprintf('max chemical stationarity res = %.8e\n',max(abs(DIAG.mu_stationarity)))
fprintf('max |omega-omega(ref)|        = %.8e\n',max(abs(omega_res)))
fprintf('sum(phase_prop)-1             = %.8e\n',DIAG.phase_sum_error)

end


function [F,PART] = Penalty_Objective(z,pars,Nc,E_target,eta)

Np                   = length(pars);
Ne                   = length(E_target);
[c,phase_prop]       = Unpack_Penalty_State(z,Nc,Np);

Etarget              = cell2mat(E_target(:));
Emix                 = zeros(Ne,1);
Gmix                 = 0;
omega                = zeros(Np,1);
R                    = cell(1,Np);

for ip = 1:Np

    R{ip}              = PhaseThermo(pars{ip},c{ip});
    e_ip               = cell2mat(R{ip}.e(:));

    Gmix                = Gmix + phase_prop(ip)*R{ip}.g(1);
    Emix                = Emix + phase_prop(ip)*e_ip;

end

res                  = Etarget-Emix;
mu                   = eta*res;
F                    = Gmix + 0.5*eta*sum(res.^2);
mu_stationarity      = [];

for ip = 1:Np

    e_ip               = cell2mat(R{ip}.e(:));
    omega(ip)          = R{ip}.g(1)-e_ip.'*mu;

    if phase_prop(ip) > 1e-8 && ~isempty(R{ip}.mu_c)
        mu_c           = cell2mat(R{ip}.mu_c(:));
        J              = R{ip}.Jac(:,:,1);
        mu_stationarity= [mu_stationarity; mu_c-J.'*mu]; %#ok<AGROW>
    end

end

if ~isfinite(F)
    F = realmax/1e100;
end

if nargout > 1
    PART.Etarget         = Etarget;
    PART.Emix            = Emix;
    PART.mu              = mu;
    PART.omega           = omega;
    PART.mu_stationarity = mu_stationarity;
end

end


function [c,phase_prop] = Unpack_Penalty_State(z,Nc,Np)

c        = cell(1,Np);
i0       = 0;

for ip = 1:Np
    ids   = i0+(1:Nc(ip));
    c{ip} = num2cell(z(ids).');
    i0    = i0+Nc(ip);
end

phase_prop = z(i0+(1:Np));

end


function [Gmix,Emix] = Evaluate_Mixture(pars,c,phase_prop)

Np   = length(pars);
Ne   = length(PhaseThermo(pars{1},c{1}).e);
Gmix = 0;
Emix = cell(1,Ne);

for ie = 1:Ne
    Emix{ie} = 0;
end

for ip = 1:Np

    R     = PhaseThermo(pars{ip},c{ip});
    Gmix  = Gmix + phase_prop(ip)*R.g(1);

    for ie = 1:Ne
        Emix{ie} = Emix{ie} + phase_prop(ip)*R.e{ie};
    end

end

end


function c = Expand_c_By_Phase(c_phase,phase_index)

Ngrain = length(phase_index);
c      = cell(1,Ngrain);

for ig = 1:Ngrain
    c{ig} = c_phase{phase_index(ig)};
end

end


function p_phase = Collapse_p_By_Phase(p_grain,phase_index,Nphase)

[ny,nx,~] = size(p_grain);
p_phase   = zeros(ny,nx,Nphase);

for ip = 1:Nphase

    grains = find(phase_index == ip);

    if ~isempty(grains)
        p_phase(:,:,ip) = sum(p_grain(:,:,grains),3);
    end

end

end


function omg = Calc_Omega_Grain(pars,c,e,mu_e,ny,nx)
%CALC_OMEGA_GRAIN Calculate omega from spatial composition fields.
%
% PhaseG requires composition in the form Nc x Npoint. The map stores
% c{ip}{ic} as ny x nx fields, so flatten them before calling PhaseG.

Np      = length(pars);
Ne      = length(mu_e);
omg     = zeros(ny,nx,Np);
c_flat  = Unpack_c(c);

for ip = 1:Np

    g_ip          = PhaseG(pars{ip},c_flat{ip});
    omg(:,:,ip)   = reshape(g_ip,ny,nx);

    for ie = 1:Ne
        omg(:,:,ip) = omg(:,:,ip)-e{ip}{ie}.*mu_e{ie};
    end

end

end


function omega_phase = Extract_Omega_By_Phase(omg,phase_index,Nphase)

[ny,nx,~]  = size(omg);
omega_phase = zeros(ny,nx,Nphase);

for ip = 1:Nphase
    ig                  = find(phase_index == ip,1,'first');
    omega_phase(:,:,ip) = omg(:,:,ig);
end

end


function Print_Omega_Difference(omega_phase,phs_name)

Nphase = length(phs_name);

for ip = 1:Nphase-1
    for jp = ip+1:Nphase

        domg = omega_phase(:,:,ip)-omega_phase(:,:,jp);

        fprintf('%s - %s: max|domg| = %.8e, mean|domg| = %.8e\n', ...
                phs_name{ip},phs_name{jp}, ...
                max(abs(domg(:))),mean(abs(domg(:))))

    end
end

end


function Plot_Omega_Difference(x,omega_phase,phs_name)

Nphase = length(phs_name);
label  = cell(1,Nphase-1);

hold on

for ip = 2:Nphase

    domg = mean(omega_phase(:,:,ip)-omega_phase(:,:,1),1);
    plot(x*1e6,domg,'LineWidth',1.2)

    label{ip-1} = [phs_name{ip},' - ',phs_name{1}];

end

hold off
grid on
xlim([min(x) max(x)]*1e6)
ylabel('\Delta\omega')
xlabel('x (\mum)')
legend(label,'Location','best')

end
