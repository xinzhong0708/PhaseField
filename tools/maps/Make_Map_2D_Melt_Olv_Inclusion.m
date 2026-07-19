%% Create_Map2D_Inclusion_Clean.m
% Generate Map2d.mat with one circular inclusion in the domain center.
%
% Grain 1 = inclusion phase
% Grain 2 = matrix phase
%
% The initial phase compositions are obtained from the reference LE state
% using the requested phase proportions.
%
% Grain orientation:
%   PARAM.theta_grain(1) is the crystal orientation of grain 1.
%   PARAM.theta_grain(2) is the matrix/melt orientation and is not used.

clear; clf

addpath('..\bin')
addpath('..\')
addpath('..\ThermoData')

%% ========================================================================
%  User controls
% ========================================================================

% Thermodynamic phases
phs_name = {'Olv','Melt'};

% Requested phase proportions
% phase_prop(1) is used to determine the inclusion radius when
% incl_radius = [].
phase_prop = [0.01 0.995];
phase_prop = phase_prop/sum(phase_prop);

% Inclusion controls
incl_phase_name   = phs_name{1};     % inclusion phase
matrix_phase_name = phs_name{2};     % matrix phase
incl_radius       = [];              % [] means radius from phase_prop of inclusion phase
incl_center_real  = [];              % [] means domain center, or [x y] in meter
incl_smooth_width = 0;               % 0 = sharp; e.g. 4*Lx/nx = smooth

% Grain orientation controls
% Grain 1 = inclusion/crystal.
% Grain 2 = matrix/melt, orientation is not used.
theta0 = 0.0;                         % rad
% theta0 = pi/8;                      % example
% theta0 = pi*rand;                   % random orientation

% Initial independent endmember compositions
c_value    = cell(1,numel(phs_name));
c_value{1} = [0.45 0.45];            % Olv
c_value{2} = [0.35 ];                % Melt

% Domain
Lx = 20e-6;
Ly = 20e-6;
nx = 100;
ny = 100;

% Scaling / penalty
PHYS      = struct();
PHYS.E_sc = 1e9;
PHYS.L_sc = 1;

E_sc = PHYS.E_sc;
L_sc = PHYS.L_sc;
eta0 = 4000e10/E_sc;

%% ========================================================================
%  Load thermodynamics and reference LE
% ========================================================================

pars_phase = Load_Data(phs_name);
Nphase     = numel(pars_phase);

incl_phase   = find(strcmp(phs_name,incl_phase_name),1,'first');
matrix_phase = find(strcmp(phs_name,matrix_phase_name),1,'first');

if isempty(incl_phase)
    error('incl_phase_name %s is not in phs_name.',incl_phase_name)
end
if isempty(matrix_phase)
    error('matrix_phase_name %s is not in phs_name.',matrix_phase_name)
end
if incl_phase == matrix_phase
    error('incl_phase_name and matrix_phase_name must be different.')
end

c_guess = cell(1,Nphase);

for ip = 1:Nphase
    c_guess{ip} = num2cell(c_value{ip});
end

% Reference LE uses requested phase proportions.
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
%  Central inclusion map
% ========================================================================

% Grain 1 is the inclusion, grain 2 is the matrix.
grain_phase = [incl_phase; matrix_phase];
Ngrain      = 2;
Np          = Ngrain;

[X,Y] = meshgrid(GRID.x,GRID.y);

if isempty(incl_center_real)
    incl_center = [mean(GRID.x),mean(GRID.y)];
else
    incl_center = incl_center_real/L_sc;
end

if isempty(incl_radius)
    incl_radius_sc = sqrt(phase_prop(incl_phase)*GRID.Lx*GRID.Ly/pi);
else
    incl_radius_sc = incl_radius/L_sc;
end

r = sqrt((X-incl_center(1)).^2 + (Y-incl_center(2)).^2);

phi      = zeros(ny,nx,Ngrain);
phase_ID = matrix_phase*ones(ny,nx);
grain_ID = 2*ones(ny,nx);
seed_xy  = zeros(Ngrain,2);

if incl_smooth_width > 0

    w  = incl_smooth_width/L_sc;
    p1 = 0.5*(1 - tanh((r - incl_radius_sc)/w));
    p2 = 1 - p1;

    phi(:,:,1) = sqrt(p1);
    phi(:,:,2) = sqrt(p2);

    mask_incl = p1 >= p2;

else

    mask_incl = r <= incl_radius_sc;

    phi(:,:,1) = double(mask_incl);
    phi(:,:,2) = double(~mask_incl);

end

phase_ID(mask_incl) = incl_phase;
grain_ID(mask_incl) = 1;

seed_xy(1,:) = incl_center;
seed_xy(2,:) = [NaN NaN];

% Grain orientations used by Calc_AC_Anisotropy_FacetedStiffness.
% Only crystal grain 1 matters here.
theta_grain    = zeros(1,Ngrain);
theta_grain(1) = theta0;
theta_grain(2) = 0.0;

phase_prop_geom = Compute_Phase_Fraction_From_Map(phase_ID,Nphase);
phase_prop_map  = phase_prop_geom;
map_mode        = 'incl';

fprintf('Central %s inclusion in %s matrix\n',incl_phase_name,matrix_phase_name)
fprintf('Inclusion radius = %.4e m\n',incl_radius_sc*L_sc)
fprintf('Inclusion area fraction = %.8f\n',nnz(mask_incl)/numel(mask_incl))
fprintf('Seed orientation theta0 = %.6f rad = %.3f degree\n', ...
    theta0,theta0*180/pi)

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
%  Plot map
% ========================================================================

figure(3); clf
imagesc(x*1e6,y*1e6,phase_ID)
set(gca,'YDir','normal')
axis image tight
colorbar
title('Initial central inclusion map')
xlabel('x (\mum)')
ylabel('y (\mum)')
drawnow

fprintf('\nInitial map mode: %s\n',map_mode)
fprintf('Number of grains = %d\n',Ngrain)
fprintf('\nPhase       requested        realised area fraction      difference\n')
for ip = 1:Nphase
    fprintf('%-8s    %.8f         %.8f                 %+.3e\n', ...
        phs_name{ip},phase_prop(ip),phase_prop_geom(ip), ...
        phase_prop_geom(ip)-phase_prop(ip))
end

%% ========================================================================
%  Construct locally penalty-consistent conserved composition field
% ========================================================================

c   = Expand_c_By_Phase(c_phase,MODEL.phase_index);
p   = Calc_p(MODEL,phi);
e   = Calc_e(pars,c);
E   = Calc_E_Tot(e,p);
eta = eta0*ones(ny,nx);

% Local penalty consistency:
%   E - Emix = mu_ref / eta
for ie = 1:Ne
    E{ie} = E{ie} + E_offset{ie};
end

E_bulk_shift = zeros(1,Ne);

mu_e = repmat({zeros(ny,nx)},1,Ne);
chi  = repmat({zeros(ny,nx)},Ne,Ne);

%% ========================================================================
%  PARAM and STATE
% ========================================================================

PARAM              = struct();
PARAM.Np           = Np;
PARAM.Ne           = Ne;
PARAM.eta          = eta;
PARAM.use_WScale   = 0;
PARAM.theta_grain  = theta_grain;

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

STATE_INI     = STATE;
PARAM.LE_mode = 'LE';
STATE         = LE_Run_Mode_New(STATE,PARAM,MODEL);

p_phase        = Collapse_p_By_Phase(STATE.p,MODEL.phase_index,Nphase);
phase_prop_map = zeros(1,Nphase);

for ip = 1:Nphase
    phase_prop_map(ip) = mean(p_phase(:,:,ip),'all');
end

omega_phase = zeros(ny,nx,Nphase);

for ip = 1:Nphase
    ig = find(MODEL.phase_index == ip,1,'first');

    if ~isempty(ig)
        omega_phase(:,:,ip) = STATE.omg(:,:,ig);
    end
end

fprintf('\nInitial LE omega differences:\n')

for ip = 1:Nphase-1
    for jp = ip+1:Nphase
        domg = omega_phase(:,:,ip) - omega_phase(:,:,jp);

        fprintf('%s - %s: max|domg| = %.8e, mean|domg| = %.8e\n', ...
            phs_name{ip},phs_name{jp}, ...
            max(abs(domg(:))),mean(abs(domg(:))))
    end
end

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
%  Completeness check before save
% ========================================================================

if ~isfield(PARAM,'theta_grain')
    error('PARAM.theta_grain is missing.')
end

if numel(PARAM.theta_grain) ~= Ngrain
    error('PARAM.theta_grain length does not match Ngrain before save.')
end

if ~isfield(STATE,'omg') || size(STATE.omg,3) ~= Ngrain
    error('STATE.omg is missing or has wrong number of grains.')
end

if ~isfield(STATE,'mu_e') || numel(STATE.mu_e) ~= Ne
    error('STATE.mu_e is missing or has wrong number of elements.')
end

if ~isfield(STATE,'chi') || size(STATE.chi,1) ~= Ne || size(STATE.chi,2) ~= Ne
    error('STATE.chi is missing or has wrong size.')
end

%% ========================================================================
%  Save
% ========================================================================

grain_size      = NaN;
grain_size_real = NaN;
Ngrain_user     = [];
rng_seed        = [];
periodic_map    = 0;

save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'E_sc','L_sc','eta','pars','Np','Ne','Nphase','Ngrain', ...
    'phs_name','phase_prop','phase_prop_geom','phase_prop_map', ...
    'map_mode','grain_ID','phase_ID','grain_phase','seed_xy','theta_grain','theta0', ...
    'incl_phase_name','matrix_phase_name','incl_radius','incl_center_real','incl_smooth_width', ...
    'grain_size','grain_size_real','Ngrain_user','rng_seed','periodic_map', ...
    'E_target','E_offset','E_bulk_shift','c_ref','mu_ref','omega_ref', ...
    'c','e','E','p','mu_e','chi','phi')

fprintf('\nSaved Map2d.mat\n')
fprintf('Saved PARAM.theta_grain:\n')
disp(PARAM.theta_grain)

%% ========================================================================
%  Local helper functions
% ========================================================================

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