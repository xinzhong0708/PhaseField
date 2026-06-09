%% Create_Map2D_NoLE_Sphere.m
% Simple 2D circular solid core map without any initial LE.
% The initial phase compositions are exactly the user-specified c_value.
%
% Phase 1 = circular solid core
% Phase 2 = surrounding matrix/liquid

clear; clf

addpath('..\bin')
addpath('..\')
addpath('..\ThermoData')

%% ========================================================================
%  User controls
% ========================================================================

% Thermodynamic phases
phs_name = {'A','B'};

% Initial independent endmember compositions
c_value    = cell(1,numel(phs_name));
c_value{1} = [1.1e-3];        % A, solid core
c_value{2} = [0.30];          % B, matrix/liquid

% Domain
Lx = 1;
Ly = 1;
nx = 100;
ny = 100;

% Circular core geometry, in the same physical unit as Lx and Ly
core_radius = 0.15;
core_center = [0.50 0.50];     % [x0 y0]

% Scaling / penalty
PHYS        = struct();
PHYS.E_sc   = 1;
PHYS.L_sc   = 1;

E_sc        = PHYS.E_sc;
L_sc        = PHYS.L_sc;
eta0        = 5000;

%% ========================================================================
%  Load thermodynamics
% ========================================================================

pars_phase = Load_Data(phs_name);
Nphase     = numel(pars_phase);

if Nphase ~= 2
    error('This simple sphere/circle initializer assumes exactly two phases.')
end

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
%  Circular core map
% ========================================================================

[X,Y] = meshgrid(x,y);

x0 = core_center(1)/L_sc;
y0 = core_center(2)/L_sc;
r0 = core_radius/L_sc;

core_mask = (X-x0).^2 + (Y-y0).^2 <= r0^2;

Ngrain      = 2;
Np          = Ngrain;
grain_phase = [1;2];

phi      = zeros(ny,nx,Ngrain);
phase_ID = zeros(ny,nx);
grain_ID = zeros(ny,nx);
seed_xy  = zeros(Ngrain,2);

% Phase/grain 1: circular core
phi(:,:,1)        = core_mask;
phase_ID(core_mask) = 1;
grain_ID(core_mask) = 1;
seed_xy(1,:)      = [x0 y0];

% Phase/grain 2: matrix
phi(:,:,2)           = ~core_mask;
phase_ID(~core_mask) = 2;
grain_ID(~core_mask) = 2;
seed_xy(2,:)         = [mean(x(~any(core_mask,1))) mean(y)];

if any(~isfinite(seed_xy(2,:)))
    seed_xy(2,:) = [mean(x) mean(y)];
end

phase_prop_geom = zeros(1,Nphase);
for ip = 1:Nphase
    phase_prop_geom(ip) = mean(phase_ID(:) == ip);
end

% For compatibility with old scripts, set requested phase_prop to realised one.
phase_prop     = phase_prop_geom;
phase_prop_map = phase_prop_geom;

%% ========================================================================
%  Grain-resolved MODEL
% ========================================================================

MODEL             = struct();
MODEL.phs_name    = phs_name;
MODEL.phase_index = grain_phase(:).';
MODEL.pars        = pars_phase(MODEL.phase_index);

eps_phi = 1e-14;

MODEL.p_fun = @(a,phi) ...
    phi(:,:,a).^2 ./ (sum(phi.^2,3) + eps_phi);

MODEL.dpdphi = @(a,b,phi) ...
    (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) ...
    - 2*phi(:,:,a).*phi(:,:,b).^2 ./ (sum(phi.^2,3) + eps_phi).^2;

pars = MODEL.pars;
F    = MODEL;

%% ========================================================================
%  Initial c, e, E directly from user-set c_value
% ========================================================================

c_phase = cell(1,Nphase);

for ip = 1:Nphase
    c_phase{ip} = cell(1,numel(c_value{ip}));
    for ic = 1:numel(c_value{ip})
        c_phase{ip}{ic} = c_value{ip}(ic)*ones(ny,nx);
    end
end

c   = Expand_c_By_Phase(c_phase,MODEL.phase_index);
p   = Calc_p(MODEL,phi);
e   = Calc_e(pars,c);
E   = Calc_E_Tot(e,p);
Ne  = numel(E);
eta = eta0*ones(ny,nx);

mu_e = repmat({zeros(ny,nx)},1,Ne);
chi  = repmat({zeros(ny,nx)},Ne,Ne);

E_target     = cell(1,Ne);
E_offset     = cell(1,Ne);
E_bulk_shift = zeros(1,Ne);

for ie = 1:Ne
    E_target{ie} = mean(E{ie}(:));
    E_offset{ie} = zeros(ny,nx);
end

%% ========================================================================
%  PARAM and STATE
% ========================================================================

PARAM            = struct();
PARAM.Np         = Np;
PARAM.Ne         = Ne;
PARAM.eta        = eta;
PARAM.use_WScale = 0;

STATE          = struct();
STATE.c        = c;
STATE.e        = e;
STATE.E        = E;
STATE.mu_e     = mu_e;
STATE.chi      = chi;
STATE.phi      = phi;
STATE.p        = p;
STATE.mask     = ones(ny,nx,Np);
STATE.LE_state = [];

%% ========================================================================
%  Plot and print
% ========================================================================

figure(3); clf
imagesc(x*1e6,y*1e6,phase_ID)
set(gca,'YDir','normal')
axis image tight
colorbar
title('Initial circular solid core')
xlabel('x (\mum)')
ylabel('y (\mum)')
drawnow

fprintf('\nInitial circular core map without LE\n')
fprintf('Number of grains = %d\n',Ngrain)
fprintf('Core radius = %.8e\n',core_radius)
fprintf('Core center = [%.8e %.8e]\n',core_center(1),core_center(2))

fprintf('\nPhase       realised area fraction\n')
for ip = 1:Nphase
    fprintf('%-8s    %.8f\n',phs_name{ip},phase_prop_geom(ip))
end

fprintf('\nInitial c values used directly:\n')
for ip = 1:Nphase
    fprintf('%s: ',phs_name{ip})
    fprintf('%.8e ',c_value{ip})
    fprintf('\n')
end

%% ========================================================================
%  Save
% ========================================================================

map_mode        = 'circle';
grain_size      = NaN;
grain_size_real = 2*core_radius;
Ngrain_user     = [];
rng_seed        = [];
periodic_map    = 0;
c_ref           = c_phase;

save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'E_sc','L_sc','eta','pars','Np','Ne','Nphase','Ngrain', ...
    'phs_name','phase_prop','phase_prop_geom','phase_prop_map', ...
    'map_mode','grain_ID','phase_ID','grain_phase','seed_xy', ...
    'grain_size','grain_size_real','Ngrain_user','rng_seed','periodic_map', ...
    'E_target','E_offset','E_bulk_shift','c_ref', ...
    'c','e','E','p','mu_e','chi','phi', ...
    'core_radius','core_center')

fprintf('\nSaved Map2d.mat\n')

%% ========================================================================
%  Local helper functions
% ========================================================================

function c = Expand_c_By_Phase(c_phase,phase_index)

Ngrain = numel(phase_index);
c      = cell(1,Ngrain);

for ig = 1:Ngrain
    c{ig} = c_phase{phase_index(ig)};
end

end
