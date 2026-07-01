%% Create_Map2D_Melt_Garnet.m
% Simple Map2d.mat with two phases and two grains:
%   grain 1 = Melt(H18), background
%   grain 2 = Garnet, circular seed in the middle
%
% Full thermodynamic phase names are used.
% Thermodynamic pars are built on the fly from T and P.
% No Data_*.mat files are loaded.

clear; clf

addpath('..\bin')
addpath('..\')
addpath('..\Thermo')
addpath('..\Thermo\Solutions')

%% ========================================================================
%  User controls
% ========================================================================

map_mode = 'melt_garnet_seed';

% Pressure-temperature initial condition
T              = 1210 + 273.15;       % K
P              = 0.5e9;              % Pa

% Elements and solution model
Cname          = {'Fe' 'Mg' 'Ca' 'Al' 'Si' 'O'};
solmod         = 'solution_models_PFM';

% Only two thermodynamic phases
phs_name       = {'Melt(H18)','Garnet'};

% Domain
Lx = 25e-6;
Ly = 25e-6;
nx = 120;
ny = 120;

% Central garnet seed
seed_phase_name  = 'Garnet';
seed_radius      = 2.00e-6;      % meter
seed_center_real = [];           % [] means domain center, or [x y] in meter

% Scaling / penalty
PHYS        = struct();
PHYS.E_sc   = 1e9;
PHYS.L_sc   = 1;
PHYS.vref   = 2e-5;

E_sc        = PHYS.E_sc;
L_sc        = PHYS.L_sc;
vref        = PHYS.vref;
eta0        = 4000e10/E_sc;

%% ========================================================================
%  Grid
% ========================================================================

x  = linspace(0,Lx,nx)/L_sc;
y  = linspace(0,Ly,ny)/L_sc;
dx = x(2)-x(1);
dy = y(2)-y(1);

GRID       = struct();
GRID.x     = x;
GRID.y     = y;
GRID.dx    = dx;
GRID.dy    = dy;
GRID.nx    = nx;
GRID.ny    = ny;
GRID.Lx    = Lx/L_sc;
GRID.Ly    = Ly/L_sc;

%% ========================================================================
%  Two-grain map
% ========================================================================

Nphase = 2;
Ngrain = 2;
Np     = 2;

grain_phase = [1 2];     % grain 1 = Melt(H18), grain 2 = Garnet

[X,Y] = meshgrid(GRID.x,GRID.y);

if isempty(seed_center_real)
    seed_center = [mean(GRID.x),mean(GRID.y)];
else
    seed_center = seed_center_real/L_sc;
end

seed_radius_sc = seed_radius/L_sc;
garnet_mask    = (X-seed_center(1)).^2 + ...
                 (Y-seed_center(2)).^2 <= seed_radius_sc^2;

if ~any(garnet_mask(:))
    error('Central garnet seed contains no grid points. Increase seed_radius.')
end

grain_ID = ones(ny,nx);
grain_ID(garnet_mask) = 2;

phase_ID = grain_ID;

phi = zeros(ny,nx,Ngrain);
phi(:,:,1) = grain_ID == 1;
phi(:,:,2) = grain_ID == 2;

seed_xy = zeros(Ngrain,2);
seed_xy(1,:) = [mean(X(grain_ID==1)),mean(Y(grain_ID==1))];
seed_xy(2,:) = seed_center;

phase_prop      = Compute_Phase_Fraction_From_Map(phase_ID,Nphase);
phase_prop_geom = phase_prop;

theta_grain      = [0 0];
theta_grain_user = theta_grain;
theta_rng_seed   = [];

% Legacy placeholders kept for save-file compatibility
grain_size          = NaN;
grain_size_real     = NaN;
Ngrain_user         = Ngrain;
rng_seed            = NaN;
periodic_map        = 0;
use_center_seed     = 1;
single_seed_only    = 1;
separate_same_phase = 0;
phase_sep_weight    = NaN;

fprintf('Central Garnet seed: radius %.4e m, area fraction %.8f\n', ...
    seed_radius,phase_prop(2))

%% ========================================================================
%  Build thermodynamics and reference LE
% ========================================================================

pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref);

% Initial independent endmember compositions.
% Melt(H18): use one interior composition generated from Thermolab props.
% Garnet: use the same simple Fe-Mg starting point as before.
td_tmp       = init_thermo(phs_name,Cname,solmod);
p_tmp        = props_generate(td_tmp);
c_value      = cell(1,Nphase);
c_value{1}   = mean(p_tmp{1}(:,1:end-1),1);     % Melt(H18)
c_value{2}   = [0.45 0.45];                     % Garnet

c_guess = cell(1,Nphase);

for ip = 1:Nphase
    c_guess{ip} = num2cell(c_value{ip});
end

% Reference LE uses the actual geometric phase proportions.
p_ref    = reshape(phase_prop,1,1,Nphase);
e_guess  = Calc_e(pars_phase,c_guess);
E_target = Calc_E_Tot(e_guess,p_ref);
Ne       = numel(E_target);

[c_ref,mu_ref] = LE_Calculator(pars_phase,p_ref,c_guess,E_target,eta0,[0.1,2000]);

e_ref     = Calc_e(pars_phase,c_ref);
E_mix_ref = Calc_E_Tot(e_ref,p_ref);

E_offset  = cell(1,Ne);
omega_ref = zeros(1,Nphase);

for ie = 1:Ne
    E_offset{ie} = E_target{ie} - E_mix_ref{ie};
end

for ip = 1:Nphase

    omega_ref(ip) = PhaseG(pars_phase{ip},c_ref{ip});

    for ie = 1:Ne
        omega_ref(ip) = omega_ref(ip) - e_ref{ip}{ie}.*mu_ref{ie};
    end
end

%% ========================================================================
%  Uniform phase compositions from reference LE
% ========================================================================

c_phase = cell(1,Nphase);

for ip = 1:Nphase
    c_phase{ip} = cell(size(c_ref{ip}));

    for ic = 1:numel(c_ref{ip})
        c_phase{ip}{ic} = c_ref{ip}{ic}*ones(ny,nx);
    end
end

%% ========================================================================
%  Grain-resolved MODEL
% ========================================================================

MODEL             = struct();
MODEL.phs_name    = phs_name;
MODEL.phase_index = grain_phase(:).';
MODEL.pars        = pars_phase(MODEL.phase_index);
MODEL.eta         = eta0;

% Full thermo information for on-the-fly P-T update
MODEL.Cname       = Cname;
MODEL.solmod      = solmod;
MODEL.T           = T;
MODEL.P           = P;
MODEL.E_sc        = E_sc;
MODEL.vref        = vref;

eps_phi = 1e-14;

MODEL.p_fun = @(a,phi) ...
    phi(:,:,a).^2 ./ (sum(phi.^2,3) + eps_phi);

MODEL.dpdphi = @(a,b,phi) ...
    (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) ...
    - 2*phi(:,:,a).*phi(:,:,b).^2 ./ (sum(phi.^2,3) + eps_phi).^2;

pars = MODEL.pars;
F    = MODEL;

%% ========================================================================
%  Plot map
% ========================================================================

figure(3); clf
imagesc(x*1e6,y*1e6,phase_ID)
set(gca,'YDir','normal')
axis image tight
colorbar
title('Initial Melt(H18) matrix with central Garnet seed')
xlabel('x (\mum)')
ylabel('y (\mum)')
drawnow

fprintf('\nInitial map mode: %s\n',map_mode)
fprintf('Number of phases = %d\n',Nphase)
fprintf('Number of grains = %d\n',Ngrain)

fprintf('\nPhase       realised area fraction\n')
for ip = 1:Nphase
    fprintf('%-18s %.8f\n',phs_name{ip},phase_prop_geom(ip))
end

fprintf('\nGrain orientations theta_grain:\n')
disp(theta_grain)

%% ========================================================================
%  Construct locally penalty-consistent conserved composition field
% ========================================================================

c   = Expand_c_By_Phase(c_phase,MODEL.phase_index);
p   = Calc_p(MODEL,phi);
e   = Calc_e(pars,c);
E   = Calc_E_Tot(e,p);
eta = eta0*ones(ny,nx);

% Local penalty consistency: E - Emix = mu_ref / eta.
for ie = 1:Ne
    E{ie} = E{ie} + E_offset{ie};
end

E_bulk_shift = zeros(1,Ne);

mu_e = repmat({zeros(ny,nx)},1,Ne);
chi  = repmat({zeros(ny,nx)},Ne,Ne);

%% ========================================================================
%  PARAM and STATE
% ========================================================================

PARAM             = struct();
PARAM.Np          = Np;
PARAM.Ne          = Ne;
PARAM.eta         = eta;
PARAM.use_WScale  = 0;
PARAM.T           = T;
PARAM.P           = P;
PARAM.theta_grain = theta_grain;

STATE          = struct();
STATE.c        = c;
STATE.e        = e;
STATE.E        = E;
STATE.mu_e     = mu_e;
STATE.chi      = chi;
STATE.omg      = zeros(ny,nx,Np);
STATE.phi      = phi;
STATE.p        = p;
STATE.mask     = ones(ny,nx,Np);
STATE.LE_state = [];

%% ========================================================================
%  Completeness check before LE
% ========================================================================

if numel(MODEL.phase_index) ~= Ngrain
    error('MODEL.phase_index length does not match Ngrain.')
end

if numel(MODEL.pars) ~= Ngrain
    error('MODEL.pars length does not match Ngrain.')
end

if size(STATE.phi,3) ~= Ngrain
    error('STATE.phi third dimension does not match Ngrain.')
end

if size(STATE.p,3) ~= Ngrain
    error('STATE.p third dimension does not match Ngrain.')
end

if numel(PARAM.theta_grain) ~= Ngrain
    error('PARAM.theta_grain length does not match Ngrain.')
end

if PARAM.Np ~= Ngrain
    error('PARAM.Np does not match Ngrain.')
end

if PARAM.Ne ~= Ne
    error('PARAM.Ne does not match Ne.')
end

%% ========================================================================
%  Initial LE check
% ========================================================================

STATE_INI = STATE;

PARAM.LE_mode = 'LE';
STATE         = LE_Run_Mode_New(STATE,PARAM,MODEL);

p_phase        = Collapse_p_By_Phase(STATE.p,MODEL.phase_index,Nphase);
phase_prop_map = zeros(1,Nphase);

for ip = 1:Nphase
    tmp = p_phase(:,:,ip);
    phase_prop_map(ip) = mean(tmp(:));
end

omega_phase = zeros(ny,nx,Nphase);

for ip = 1:Nphase
    grains = find(MODEL.phase_index == ip);

    if ~isempty(grains)
        omega_phase(:,:,ip) = mean(STATE.omg(:,:,grains),3);
    end
end

fprintf('\nInitial LE omega differences:\n')

domg = omega_phase(:,:,1) - omega_phase(:,:,2);
fprintf('%s - %s: max|domg| = %.8e, mean|domg| = %.8e\n', ...
    phs_name{1},phs_name{2},max(abs(domg(:))),mean(abs(domg(:))))

dc_max = 0;

for ig = 1:Ngrain
    for ic = 1:numel(STATE.c{ig})
        dc_max = max(dc_max, ...
            max(abs(STATE.c{ig}{ic}(:)-STATE_INI.c{ig}{ic}(:))));
    end
end

E_mean        = zeros(Ne,1);
E_target_vec  = zeros(Ne,1);
mu_jump       = zeros(1,Ne);

for ie = 1:Ne
    E_mean(ie)       = mean(STATE.E{ie}(:));
    E_target_vec(ie) = E_target{ie};
    mu_jump(ie)      = max(STATE.mu_e{ie}(:)) - min(STATE.mu_e{ie}(:));
end

fprintf('\nmax|c after first LE - c before LE| = %.8e\n',dc_max)
fprintf('max|mean(E map) - E target|         = %.8e\n', ...
    max(abs(E_mean-E_target_vec)))
fprintf('initial mu_e jump after LE:\n')
disp(mu_jump)

% Keep legacy variables consistent
c    = STATE.c;
e    = STATE.e;
E    = STATE.E;
p    = STATE.p;
mu_e = STATE.mu_e;
chi  = STATE.chi;
phi  = STATE.phi;

%% ========================================================================
%  Save
% ========================================================================

save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'E_sc','L_sc','vref','eta','pars','Np','Ne','Nphase','Ngrain', ...
    'T','P','Cname','solmod', ...
    'phs_name','phase_prop','phase_prop_geom','phase_prop_map', ...
    'map_mode','grain_ID','phase_ID','grain_phase','seed_xy','theta_grain', ...
    'theta_grain_user','theta_rng_seed', ...
    'grain_size','grain_size_real','Ngrain_user','rng_seed','periodic_map', ...
    'use_center_seed','seed_phase_name','seed_radius','seed_center_real','single_seed_only', ...
    'separate_same_phase','phase_sep_weight', ...
    'E_target','E_offset','E_bulk_shift','c_ref','mu_ref','omega_ref', ...
    'c_value','c','e','E','p','mu_e','chi','phi')

fprintf('\nSaved Map2d.mat\n')
fprintf('Saved full-name phases:\n')
disp(MODEL.phs_name)
fprintf('Saved MODEL.phase_index:\n')
disp(MODEL.phase_index)
fprintf('Saved PARAM.theta_grain:\n')
disp(PARAM.theta_grain)

%% ========================================================================
%  Local helper functions
% ========================================================================

function pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref)

%BUILD_PARS_PHASES Build pars for present phases at given T and P.
%
% Full-name-only version. No Data_*.mat files are loaded.
%
% Recomputes:
%   td = init_thermo(...)
%   g0 = tl_g0(T,P,td)
%   n  = td.n_em(:,1:end-1)

Nphase     = numel(phs_name);
pars_phase = cell(1,Nphase);

for ip = 1:Nphase

    phase_name = phs_name(ip);

    td = init_thermo(phase_name,Cname,solmod);
    g0 = cell2mat(tl_g0(T,P,td));
    n  = td.n_em(:,1:end-1);

    td.n_em(:,1:end-1) = n;

    pars            = td;
    pars.n          = n;
    pars.P          = P;
    pars.T          = T;
    pars.g0         = g0;
    pars.E_sc       = E_sc;
    pars.vref       = vref;
    pars.phase_name = phase_name;

    pars_phase{ip}  = pars;

end

end


function phase_prop_geom = Compute_Phase_Fraction_From_Map(phase_ID,Nphase)

phase_prop_geom = zeros(1,Nphase);
N = numel(phase_ID);

for ip = 1:Nphase
    phase_prop_geom(ip) = nnz(phase_ID == ip)/N;
end

end


function c = Expand_c_By_Phase(c_phase,phase_index)

Ngrain = numel(phase_index);
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
