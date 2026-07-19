%% Make_Map_2D_Garnet_Seed.m
% Generate Map2d.mat for either:
%   map_mode = 'bands'
%   map_mode = 'polygon'
%
% Full thermodynamic phase names are used.
% No Data_*.mat files are loaded.
%
% Thermodynamic pars are built on the fly from:
%
%   init_thermo(...)
%   tl_g0(T,P,...)
%
% Each grain has an orientation:
%
%   PARAM.theta_grain(ig)
%
% This is used later by:
%
%   Calc_AC_Anisotropy_FacetedStiffness

clear; clf

addpath('..\bin')
addpath('..\')
addpath('..\Thermo')
addpath('..\Thermo\Solutions')

%% ========================================================================
%  User controls
% ========================================================================

map_mode = 'bands';     % 'bands' or 'polygon'

% Pressure-temperature initial condition.  Keep this synchronized with the
% PF metadata by default so map construction, Phase_Diagram.m, and
% Run_2D_Scaled.m use the same thermodynamic state.
metadata_file          = '..\Metadata.xlsx';
sync_PT_from_metadata  = 1;
T              = 500 + 273.15;       % K
P              = 0.3e9;              % Pa

if sync_PT_from_metadata == 1
    [T,P] = Read_First_PT_From_Metadata(metadata_file,T,P);
end

% Elements and solution model
Cname          = {'Fe' 'Mg' 'Ca' 'Al' 'Si' 'O'};
solmod         = 'solution_models_PFM';

% Full thermodynamic phase names used in this PFM run
phs_name = {'Orthopyroxene','Garnet','Kyanite','Quartz'};

% Requested phase proportions
phase_prop = [0.3 0.04 0.2 0.1];
phase_prop = phase_prop/sum(phase_prop);

% Initial independent endmember compositions
c_value = cell(1,numel(phs_name));
c_value{1} = [0.05 0.15 0  0.4];          % Clinopyroxene
c_value{2} = [0.45 0.45];                 % Garnet
c_value{3} = 1.0;                         % Kyanite
c_value{4} = 1.0;                         % Quartz

% Domain
Lx = 100e-6;
Ly = 100e-6;

switch lower(map_mode)
    case 'bands'
        nx = 300;
        ny = 4;
    case 'polygon'
        nx = 80;
        ny = 80;
    otherwise
        error('Unknown map_mode: %s',map_mode)
end

% Band controls
Nrepeat    = 1;
band_order = repmat(1:numel(phs_name),1,Nrepeat);

% Polygon controls
grain_size   = 55.0e-6;   % mean equivalent-circle diameter, meter
Ngrain_user  = [];        % [] means estimate from grain_size
rng_seed     = 6;
periodic_map = 0;         % keep 0 unless the PF solver is periodic

% Clean polygon controls
seed_min_dist_fac  = 0.45;
lloyd_iter         = 12;
min_grain_pixels   = 8;
island_clean_iter  = 6;
quad_clean_iter    = 6;
phase_assign_tries = 80;

% Central Garnet seed controls, polygon mode only.
% [] for seed_radius means choose the radius from the requested Garnet area.
use_center_seed  = 1;
seed_phase_name  = 'Garnet';
seed_radius      = [];
seed_center_real = [];    % [] means domain center, or [x y] in meter
single_seed_only = 1;     % do not assign Garnet to any non-seed grain

% Phase separation controls, polygon mode only
separate_same_phase = 1;
phase_sep_weight    = 1000; % kept for save-file compatibility

% Smooth polygonal phi to create a diffuse initial interface.
smooth_phi            = 1;
phi_smooth_sigma_grid = 2.0;
phi_tail_cut          = 1e-10;

% Grain orientation controls
% [] means zero orientation for bands and random orientation for polygon.
% To set manually, use a row vector with length Ngrain.
theta_grain_user = [];
theta_rng_seed   = 1002;

% Scaling / penalty
PHYS        = struct();
PHYS.E_sc   = 1e9;
PHYS.L_sc   = 1;
PHYS.vref   = 2e-5;

E_sc        = PHYS.E_sc;
L_sc        = PHYS.L_sc;
vref        = PHYS.vref;
eta0        = 50000e10/E_sc;

% Initial eta used to build E.  Match this to the first eta map in
% Run_2D_Scaled; otherwise the first LE sees an artificial eta jump.
sync_eta_controls_from_metadata = 1;
init_eta_mode          = 'interface_damped'; % 'uniform' or 'interface_damped'
init_eta_damp_factor   = 0.05;
init_eta_q2            = 4;
init_eta_p02           = 3e-3;
init_eta_q3            = 4;
init_eta_p03           = 3e-3;
init_eta_nsmooth       = 1;
init_eta_halo_cut      = 1e-3;

% Reference LE controls.  Use the requested modal proportion and bulk eta so
% bands and polygon maps start from the same thermodynamic c/mu reference.
init_ref_phase_prop_mode = 'requested'; % 'requested' or 'realised_p'
init_ref_eta_mode        = 'bulk';      % 'bulk' or 'harmonic_init'

if sync_eta_controls_from_metadata == 1
    [eta0,init_eta_damp_factor] = Read_Init_Eta_From_Metadata( ...
        metadata_file,eta0,init_eta_damp_factor,E_sc);
end

%% ========================================================================
%  Build thermodynamics on the fly
% ========================================================================

pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref);
Nphase     = numel(pars_phase);

c_guess = cell(1,Nphase);

for ip = 1:Nphase
    c_guess{ip} = num2cell(c_value{ip});
end

% Reference LE is computed after the actual map p and initial eta are known.
% This avoids building an E field that is compatible with one eta and then
% immediately running the PF code with another.
e_guess = Calc_e(pars_phase,c_guess);
E_guess = Calc_E_Tot(e_guess,reshape(phase_prop,1,1,Nphase));
Ne      = numel(E_guess);

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
%  Generate phase / grain map
% ========================================================================

switch lower(map_mode)

    case 'bands'

        Noccur        = accumarray(band_order(:),1,[Nphase,1]).';
        band_fraction = phase_prop(band_order)./Noccur(band_order);
        band_fraction = band_fraction/sum(band_fraction);

        grain_phase = band_order(:);
        Ngrain      = numel(grain_phase);
        Np          = Ngrain;

        band_width = Local_Band_Widths(band_fraction,nx);

        phi      = zeros(ny,nx,Ngrain);
        phase_ID = zeros(ny,nx);
        grain_ID = zeros(ny,nx);
        seed_xy  = zeros(Ngrain,2);

        ix1 = 1;

        for ig = 1:Ngrain

            iph = grain_phase(ig);
            ix2 = ix1 + band_width(ig) - 1;

            phi(:,ix1:ix2,ig)   = 1;
            phase_ID(:,ix1:ix2) = iph;
            grain_ID(:,ix1:ix2) = ig;

            seed_xy(ig,1) = mean(x(ix1:ix2));
            seed_xy(ig,2) = mean(y);

            ix1 = ix2 + 1;

        end

        Check_Initial_Phi_Map(phi,grain_ID)

    case 'polygon'

        grain_size_sc  = grain_size/L_sc;
        grain_area_tar = pi/4*grain_size_sc^2;

        if isempty(Ngrain_user)
            Ngrain = max(nnz(phase_prop > 0), ...
                round(GRID.Lx*GRID.Ly/grain_area_tar));
        else
            Ngrain = Ngrain_user;
        end

        rng(rng_seed,'twister')
        min_seed_dist = seed_min_dist_fac*grain_size_sc;

        [grain_ID,seed_xy] = Generate_Polygonal_Grains_2D_Clean( ...
            GRID,Ngrain,periodic_map,min_seed_dist,lloyd_iter, ...
            min_grain_pixels,island_clean_iter,quad_clean_iter);

        fixed_grain = [];
        fixed_phase = [];

        if use_center_seed == 1

            seed_phase = find(strcmp(phs_name,seed_phase_name),1,'first');

            if isempty(seed_phase)
                error('seed_phase_name %s is not in phs_name.',seed_phase_name)
            end

            [X,Y] = meshgrid(GRID.x,GRID.y);

            if isempty(seed_center_real)
                seed_center = [mean(GRID.x),mean(GRID.y)];
            else
                seed_center = seed_center_real/L_sc;
            end

            if isempty(seed_radius)
                seed_radius_use = sqrt(phase_prop(seed_phase)*Lx*Ly/pi);
            else
                seed_radius_use = seed_radius;
            end

            seed_radius_sc = seed_radius_use/L_sc;
            seed_mask      = (X-seed_center(1)).^2 + ...
                             (Y-seed_center(2)).^2 <= seed_radius_sc^2;

            if ~any(seed_mask(:))
                error('Central seed contains no grid points. Increase seed_radius.')
            end

            fixed_grain = size(seed_xy,1) + 1;
            fixed_phase = seed_phase;

            grain_ID(seed_mask) = fixed_grain;
            seed_xy(fixed_grain,:) = seed_center;

            fprintf('Central %s seed: radius %.4e m, raw area fraction %.8f\n', ...
                seed_phase_name,seed_radius_use,nnz(seed_mask)/numel(seed_mask))
        end

        [grain_ID,seed_xy,fixed_grain] = Clean_Grain_Map_Topology( ...
            grain_ID,seed_xy,GRID,min_grain_pixels,island_clean_iter, ...
            quad_clean_iter,fixed_grain);

        Ngrain = size(seed_xy,1);
        Np     = Ngrain;

        phase_prop_assign = phase_prop;
        forbidden_phase   = [];

        if use_center_seed == 1 && single_seed_only == 1 && ~isempty(fixed_phase)

            seed_frac  = nnz(grain_ID == fixed_grain)/numel(grain_ID);
            rest_phase = setdiff(1:Nphase,fixed_phase);
            rest_sum   = sum(phase_prop(rest_phase));

            phase_prop_assign(:) = 0;
            phase_prop_assign(fixed_phase) = seed_frac;

            if rest_sum > 0
                phase_prop_assign(rest_phase) = phase_prop(rest_phase)/rest_sum*(1-seed_frac);
            end

            forbidden_phase = fixed_phase;

            fprintf('Only one %s seed is allowed. Non-seed grains cannot be %s.\n', ...
                seed_phase_name,seed_phase_name)
            fprintf('Central seed clean area fraction = %.8f\n',seed_frac)
        end

        if separate_same_phase == 1

            [grain_phase,~] = Assign_Grain_Phases_By_Area_Separated_Strict( ...
                grain_ID,phase_prop_assign,fixed_grain,fixed_phase, ...
                forbidden_phase,phase_assign_tries,rng_seed+1000);

        elseif use_center_seed == 1

            [grain_phase,~] = Assign_Grain_Phases_By_Area_FixedSeed( ...
                grain_ID,phase_prop_assign,fixed_grain,fixed_phase,forbidden_phase);

        else

            [grain_phase,~] = Assign_Grain_Phases_By_Area(grain_ID,phase_prop);

        end

        same_contact = Count_Same_Phase_Contacts(grain_ID,grain_phase);
        quad_count   = Count_Quad_Junctions(grain_ID);

        fprintf('Same-phase grain contacts after assignment = %d\n',same_contact)
        fprintf('Remaining 2x2 quadruple junctions          = %d\n',quad_count)

        phi_raw  = zeros(ny,nx,Ngrain);
        phase_ID = zeros(ny,nx);

        for ig = 1:Ngrain
            mask              = grain_ID == ig;
            phi_raw(:,:,ig)   = mask;
            phase_ID(mask)    = grain_phase(ig);
        end

        Check_Initial_Phi_Map(phi_raw,grain_ID)

        if smooth_phi == 1
            phi = Smooth_Initial_Phi(phi_raw,phi_smooth_sigma_grid,phi_tail_cut);
            Check_Smoothed_Phi_Map(phi)
        else
            phi = phi_raw;
        end

    otherwise

        error('Unknown map_mode: %s',map_mode)

end

phase_prop_geom = Compute_Phase_Fraction_From_Map(phase_ID,Nphase);

% Grain orientations
theta_grain = zeros(1,Ngrain);

if isempty(theta_grain_user)

    if strcmpi(map_mode,'bands')
        theta_grain(:) = 0;
    else
        rng_state = rng;

        if ~isempty(theta_rng_seed)
            rng(theta_rng_seed,'twister')
        end

        theta_grain = pi*rand(1,Ngrain);
        rng(rng_state)
    end

else

    if numel(theta_grain_user) ~= Ngrain
        error('theta_grain_user must be [] or length Ngrain.')
    end

    theta_grain = theta_grain_user(:).';

end

fprintf('Initial %s map\n',map_mode)
fprintf('Number of grains = %d\n',Ngrain)

if strcmpi(map_mode,'bands')
    fprintf('\nBand positions and orientations:\n')
    for ig = 1:Ngrain
        fprintf('grain %4d: phase = %-16s x = %.4e m, theta = %.6f rad\n', ...
            ig,phs_name{grain_phase(ig)},seed_xy(ig,1)*L_sc,theta_grain(ig))
    end
else
    grain_size_real = sqrt(4*(GRID.Lx*GRID.Ly/Ngrain)/pi)*L_sc;
    fprintf('Requested grain diameter = %.4e m\n',grain_size)
    fprintf('Realised mean equivalent diameter = %.4e m\n',grain_size_real)
    fprintf('phi smoothing sigma = %.3f grids\n',phi_smooth_sigma_grid)
end

%% ========================================================================
%  Grain-resolved MODEL
% ========================================================================

MODEL             = struct();
MODEL.phs_name    = phs_name;
MODEL.phase_index = grain_phase(:).';
MODEL.pars        = pars_phase(MODEL.phase_index);

% Thermodynamic information for on-the-fly P-T update
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
if strcmpi(map_mode,'bands')
    title('Initial pseudo-1D phase bands')
else
    title('Initial polygonal phase map')
    hold on
    plot(seed_xy(:,1)*L_sc*1e6,seed_xy(:,2)*L_sc*1e6, ...
        'k.','MarkerSize',5)
    hold off
end
xlabel('x (\mum)')
ylabel('y (\mum)')
drawnow

fprintf('\nInitial map mode: %s\n',map_mode)
fprintf('Number of grains = %d\n',Ngrain)

if strcmpi(map_mode,'bands')
    grain_size_real = NaN;
end

fprintf('\nPhase             requested        realised area fraction      difference\n')
for ip = 1:Nphase
    fprintf('%-16s    %.8f         %.8f                 %+.3e\n', ...
        phs_name{ip},phase_prop(ip),phase_prop_geom(ip), ...
        phase_prop_geom(ip)-phase_prop(ip))
end

%% ========================================================================
%  Construct locally penalty-consistent conserved composition field
% ========================================================================

p   = Calc_p(MODEL,phi);

p_phase_ref   = Collapse_p_By_Phase(p,MODEL.phase_index,Nphase);
phase_prop_ref = zeros(1,Nphase);

for ip = 1:Nphase
    tmp = p_phase_ref(:,:,ip);
    phase_prop_ref(ip) = mean(tmp(:));
end

phase_prop_ref = phase_prop_ref/sum(phase_prop_ref);

eta_bulk = eta0*ones(ny,nx);

switch lower(init_eta_mode)
    case 'uniform'
        eta_init = eta_bulk;
    case 'interface_damped'
        eta_init = Eta_Damping_SmoothHalo( ...
            p_phase_ref,eta_bulk,init_eta_damp_factor*eta_bulk, ...
            init_eta_q2,init_eta_p02,init_eta_q3,init_eta_p03, ...
            init_eta_damp_factor*eta_bulk,init_eta_nsmooth, ...
            init_eta_damp_factor*eta_bulk,init_eta_halo_cut);
    otherwise
        error('Unknown init_eta_mode: %s',init_eta_mode)
end

switch lower(init_ref_phase_prop_mode)
    case 'requested'
        phase_prop_le_ref = phase_prop;
    case 'realised_p'
        phase_prop_le_ref = phase_prop_ref;
    otherwise
        error('Unknown init_ref_phase_prop_mode: %s',init_ref_phase_prop_mode)
end

phase_prop_le_ref = phase_prop_le_ref/sum(phase_prop_le_ref);

switch lower(init_ref_eta_mode)
    case 'bulk'
        eta_ref = eta0;
    case 'harmonic_init'
        eta_ref = 1/mean(1./eta_init(:));
    otherwise
        error('Unknown init_ref_eta_mode: %s',init_ref_eta_mode)
end

eta = eta_bulk;

p_ref    = reshape(phase_prop_le_ref,1,1,Nphase);
E_target = Calc_E_Tot(e_guess,p_ref);

[c_ref,mu_ref] = LE_Calculator(pars_phase,p_ref,c_guess,E_target,eta_ref,[0.1,2000]);

e_ref     = Calc_e(pars_phase,c_ref);
E_mix_ref = Calc_E_Tot(e_ref,p_ref);

E_offset_ref = cell(1,Ne);
E_offset     = cell(1,Ne);
omega_ref    = zeros(1,Nphase);

for ie = 1:Ne
    E_offset_ref{ie} = E_target{ie} - E_mix_ref{ie};
    E_offset{ie}     = mu_ref{ie}./eta_init;
end

for ip = 1:Nphase

    omega_ref(ip) = PhaseG(pars_phase{ip},c_ref{ip});

    for ie = 1:Ne
        omega_ref(ip) = omega_ref(ip) - e_ref{ip}{ie}.*mu_ref{ie};
    end

end

c_phase = cell(1,Nphase);

for ip = 1:Nphase
    c_phase{ip} = cell(size(c_ref{ip}));

    for ic = 1:numel(c_ref{ip})
        c_phase{ip}{ic} = c_ref{ip}{ic}*ones(ny,nx);
    end
end

c = Expand_c_By_Phase(c_phase,MODEL.phase_index);
e = Calc_e(pars,c);
E = Calc_E_Tot(e,p);

% Local penalty consistency for the actual initial eta map:
% E - Emix = mu_ref / eta(x,y).
for ie = 1:Ne
    E{ie} = E{ie} + E_offset{ie};
end

E_bulk_shift = zeros(1,Ne);

mu_e = cell(1,Ne);

for ie = 1:Ne
    mu_e{ie} = mu_ref{ie}*ones(ny,nx);
end

chi  = repmat({zeros(ny,nx)},Ne,Ne);

fprintf('Initial eta mode = %s, eta range = %.4e to %.4e, eta_ref = %.4e\n', ...
    init_eta_mode,min(eta_init(:)),max(eta_init(:)),eta_ref)
fprintf('Reference LE phase proportions (%s):',init_ref_phase_prop_mode)
fprintf(' %.8f',phase_prop_le_ref)
fprintf('\n')

%% ========================================================================
%  PARAM and STATE
% ========================================================================

PARAM            = struct();
PARAM.Np         = Np;
PARAM.Ne         = Ne;
PARAM.eta        = eta_bulk;
PARAM.eta_bulk   = eta_bulk;
PARAM.eta_init   = eta_init;
PARAM.use_WScale = 0;
PARAM.T          = T;
PARAM.P          = P;
PARAM.init_eta_mode        = init_eta_mode;
PARAM.init_eta_damp_factor = init_eta_damp_factor;
PARAM.init_eta_synced_E    = 1;
PARAM.init_ref_phase_prop_mode = init_ref_phase_prop_mode;
PARAM.init_ref_eta_mode        = init_ref_eta_mode;

% Grain orientations for anisotropic faceting
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
PARAM_LE      = PARAM;
PARAM_LE.eta  = eta_init;

if exist('LE_Run_Mode','file') == 2
    STATE = LE_Run_Mode(STATE,PARAM_LE,MODEL);
else
    STATE = LE_Run_Mode_New(STATE,PARAM_LE,MODEL);
end

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

mu_jump = zeros(1,Ne);
E_mu_consistency = zeros(1,Ne);
E_mix_check = Calc_E_Tot(STATE.e,STATE.p);

for ie = 1:Ne
    mu_jump(ie) = max(STATE.mu_e{ie}(:)) - min(STATE.mu_e{ie}(:));
    E_mu_consistency(ie) = max(abs( ...
        eta_init(:).*(STATE.E{ie}(:)-E_mix_check{ie}(:)) - STATE.mu_e{ie}(:)));
end

fprintf('\nmax|c after first LE - c before LE| = %.8e\n',dc_max)
fprintf('max|eta_init*(E-Emix)-mu_e|        = %.8e\n',max(E_mu_consistency))
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
    'E_sc','L_sc','vref','eta','eta_bulk','eta_init', ...
    'pars','Np','Ne','Nphase','Ngrain', ...
    'T','P','Cname','solmod', ...
    'phs_name','phase_prop','phase_prop_geom','phase_prop_ref', ...
    'phase_prop_le_ref','phase_prop_map', ...
    'map_mode','grain_ID','phase_ID','grain_phase','seed_xy','theta_grain', ...
    'theta_grain_user','theta_rng_seed', ...
    'grain_size','grain_size_real','Ngrain_user','rng_seed','periodic_map', ...
    'seed_min_dist_fac','lloyd_iter','min_grain_pixels','island_clean_iter', ...
    'quad_clean_iter','phase_assign_tries', ...
    'use_center_seed','seed_phase_name','seed_radius','seed_center_real', ...
    'single_seed_only','separate_same_phase','phase_sep_weight', ...
    'smooth_phi','phi_smooth_sigma_grid','phi_tail_cut', ...
    'init_ref_phase_prop_mode','init_ref_eta_mode', ...
    'init_eta_mode','init_eta_damp_factor','init_eta_q2','init_eta_p02', ...
    'init_eta_q3','init_eta_p03','init_eta_nsmooth','init_eta_halo_cut', ...
    'eta_ref','E_target','E_offset','E_offset_ref','E_bulk_shift', ...
    'c_ref','mu_ref','omega_ref', ...
    'c','e','E','p','mu_e','chi','phi')

fprintf('\nSaved Map2d.mat\n')
fprintf('Saved full-name phases:\n')
disp(MODEL.phs_name)
fprintf('Saved PARAM.theta_grain:\n')
disp(PARAM.theta_grain)

%% ========================================================================
%  Local helper functions
% ========================================================================

function [T,P] = Read_First_PT_From_Metadata(metadata_file,T_default,P_default)

T = T_default;
P = P_default;

if exist(metadata_file,'file') ~= 2
    if exist('Metadata.xlsx','file') == 2
        metadata_file = 'Metadata.xlsx';
    elseif exist('MetaData.xlsx','file') == 2
        metadata_file = 'MetaData.xlsx';
    end
end

if exist(metadata_file,'file') ~= 2
    warning('Metadata file "%s" was not found. Using map P-T controls.',metadata_file)
    return
end

try
    Cpt = readcell(metadata_file,'Sheet','PTt_Path');
catch ME
    warning('Could not read PTt_Path from metadata file "%s": %s. Using map P-T controls.', ...
        metadata_file,ME.message)
    return
end

header_row = [];

for i = 1:size(Cpt,1)

    row_txt = cell(1,size(Cpt,2));

    for j = 1:size(Cpt,2)
        row_txt{j} = strtrim(Cell_String_Local(Cpt{i,j}));
    end

    if any(strcmpi(row_txt,'P (Pa)')) && ...
       any(strcmpi(row_txt,'T (K)'))  && ...
       any(strcmpi(row_txt,'t (s)'))
        header_row = i;
        break
    end
end

if isempty(header_row)
    warning('PTt_Path header row not found in "%s". Using map P-T controls.',metadata_file)
    return
end

colP = [];
colT = [];
colt = [];

for j = 1:size(Cpt,2)

    txt = strtrim(Cell_String_Local(Cpt{header_row,j}));

    if strcmpi(txt,'P (Pa)')
        colP = j;
    elseif strcmpi(txt,'T (K)')
        colT = j;
    elseif strcmpi(txt,'t (s)')
        colt = j;
    end
end

pt_rows = [];

for i = header_row+1:size(Cpt,1)

    if Is_Empty_Cell_Local(Cpt{i,colP}) || ...
       Is_Empty_Cell_Local(Cpt{i,colT}) || ...
       Is_Empty_Cell_Local(Cpt{i,colt})
        continue
    end

    pt_rows(end+1,:) = [ ...
        Cell_Number_Local(Cpt{i,colt}), ...
        Cell_Number_Local(Cpt{i,colT}), ...
        Cell_Number_Local(Cpt{i,colP})]; %#ok<AGROW>
end

if isempty(pt_rows)
    warning('PTt_Path contains no valid rows in "%s". Using map P-T controls.',metadata_file)
    return
end

[~,ord] = sort(pt_rows(:,1));
pt_rows = pt_rows(ord,:);

T = pt_rows(1,2);
P = pt_rows(1,3);

fprintf('Map P-T from %s: T = %.12g K, P = %.12g Pa\n',metadata_file,T,P)

end


function [eta0,init_eta_damp_factor] = Read_Init_Eta_From_Metadata( ...
    metadata_file,eta0,init_eta_damp_factor,E_sc)

if exist(metadata_file,'file') ~= 2
    if exist('Metadata.xlsx','file') == 2
        metadata_file = 'Metadata.xlsx';
    elseif exist('MetaData.xlsx','file') == 2
        metadata_file = 'MetaData.xlsx';
    end
end

if exist(metadata_file,'file') ~= 2
    warning('Metadata file "%s" was not found. Using map eta controls.',metadata_file)
    return
end

try
    Cmain = readcell(metadata_file,'Sheet','Main');
catch ME
    warning('Could not read metadata file "%s": %s. Using map eta controls.', ...
        metadata_file,ME.message)
    return
end

eta_SI = Get_Metadata_Value_Any(Cmain,{'Penalty eta','eta'},[]);

if ~isempty(eta_SI) && isnumeric(eta_SI)
    eta0 = eta_SI/E_sc;
end

init_eta_damp_factor = Get_Metadata_Value_Any( ...
    Cmain,{'interface damping factor','Interface damping factor'}, ...
    init_eta_damp_factor);

if ~isnumeric(init_eta_damp_factor) || ~isfinite(init_eta_damp_factor)
    error('Metadata interface damping factor must be a finite number.')
end

init_eta_damp_factor = min(max(init_eta_damp_factor,eps),1);

end


function val = Get_Metadata_Value_Any(C,keys,default)

val = default;

for ik = 1:numel(keys)

    key = strtrim(keys{ik});

    for i = 1:size(C,1)

        key_i = strtrim(Cell_String_Local(C{i,1}));

        if strcmpi(key_i,key)

            if size(C,2) < 2 || Is_Empty_Cell_Local(C{i,2})
                return
            end

            raw = C{i,2};

            if isnumeric(raw)
                val = raw;
            elseif ischar(raw) || isstring(raw)

                tmp = str2double(raw);

                if isnan(tmp)
                    val = char(raw);
                else
                    val = tmp;
                end

            else
                val = raw;
            end

            return
        end
    end
end

end


function val = Cell_Number_Local(x)

if isnumeric(x)
    val = x;
elseif ischar(x) || isstring(x)
    val = str2double(x);
else
    val = NaN;
end

if ~isfinite(val)
    error('Metadata PTt_Path contains a nonnumeric P-T-time value.')
end

end


function s = Cell_String_Local(x)

if ismissing(x)
    s = '';
elseif ischar(x)
    s = x;
elseif isstring(x)
    s = char(x);
elseif isnumeric(x)
    s = num2str(x);
else
    s = '';
end

end


function tf = Is_Empty_Cell_Local(x)

if isempty(x)
    tf = true;
elseif ismissing(x)
    tf = true;
elseif ischar(x) || isstring(x)
    tf = strlength(string(x)) == 0;
else
    tf = false;
end

end


function pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref)

%BUILD_PARS_PHASES Build pars for present phases at given T and P.
%
% Full-name-only version.
%
% Recomputes:
%   td = init_thermo(...)
%   g0 = tl_g0(T,P,td)
%   n  = td.n_em(:,1:end-1)

Nphase     = numel(phs_name);
pars_phase = cell(1,Nphase);

for ip = 1:Nphase

    % ------------------------------------------------------------
    % Phase
    % ------------------------------------------------------------
    phase_name = phs_name(ip);

    td = init_thermo(phase_name,Cname,solmod);
    g0 = cell2mat(tl_g0(T,P,td));
    n  = td.n_em(:,1:end-1);

    % ------------------------------------------------------------
    % Thermodynamic values
    % ------------------------------------------------------------
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


function band_width = Local_Band_Widths(band_fraction,nx)

%LOCAL_BAND_WIDTHS Integer band widths with no missing positive band.

band_fraction = band_fraction(:).';
band_fraction = band_fraction/sum(band_fraction);

Nband = numel(band_fraction);
raw   = nx*band_fraction;

band_width = floor(raw);
band_width(raw > 0 & band_width < 1) = 1;

while sum(band_width) < nx

    deficit = nx - sum(band_width);
    res     = raw - band_width;
    [~,ord] = sort(res,'descend');

    for k = 1:min(deficit,Nband)
        band_width(ord(k)) = band_width(ord(k)) + 1;
    end

end

while sum(band_width) > nx

    excess = sum(band_width) - nx;
    res    = raw - band_width;
    [~,ord] = sort(res,'ascend');

    changed = 0;

    for k = 1:Nband

        ig = ord(k);

        if band_width(ig) > 1
            band_width(ig) = band_width(ig) - 1;
            changed = changed + 1;
        end

        if changed >= excess
            break
        end

    end

    if changed == 0
        error('Local_Band_Widths: cannot assign positive width to all bands.')
    end

end

end


function phi = Smooth_Initial_Phi(phi_raw,sigma_grid,tail_cut)

if nargin < 2 || isempty(sigma_grid) || sigma_grid <= 0
    phi = phi_raw;
    return
end

if nargin < 3 || isempty(tail_cut)
    tail_cut = 0;
end

[ny,nx,Ngrain] = size(phi_raw);
rad = max(1,ceil(3*sigma_grid));
[xk,yk] = meshgrid(-rad:rad,-rad:rad);
K = exp(-(xk.^2 + yk.^2)/(2*sigma_grid^2));
K = K/sum(K(:));

% Boundary-normalized convolution avoids artificial phi loss at non-periodic edges.
W = conv2(ones(ny,nx),K,'same');
W = max(W,eps);

phi = zeros(ny,nx,Ngrain);
for ig = 1:Ngrain
    phi(:,:,ig) = conv2(double(phi_raw(:,:,ig)),K,'same')./W;
end

if tail_cut > 0
    phi(phi < tail_cut) = 0;
end

S = sum(phi,3);
empty = S <= eps;

if any(empty(:))
    [~,gid] = max(phi_raw,[],3);
    for ig = 1:Ngrain
        tmp = phi(:,:,ig);
        tmp(empty & gid == ig) = 1;
        phi(:,:,ig) = tmp;
    end
    S = sum(phi,3);
end

phi = phi./max(S,eps);

end


function Check_Smoothed_Phi_Map(phi)

if any(~isfinite(phi(:)))
    error('Check_Smoothed_Phi_Map: phi contains non-finite values.')
end

if min(phi(:)) < -1e-14
    error('Check_Smoothed_Phi_Map: phi contains negative values.')
end

S = sum(phi,3);
err = max(abs(S(:)-1));

if err > 1e-10
    error('Check_Smoothed_Phi_Map: sum(phi,3) is not one. max error %.3e',err)
end

end


function [grain_ID,seed_xy] = Generate_Polygonal_Grains_2D_Clean( ...
    GRID,Ngrain,periodic_map,min_seed_dist,lloyd_iter, ...
    min_grain_pixels,island_clean_iter,quad_clean_iter)

seed_xy = Generate_Seeds_MinDist(GRID,Ngrain,min_seed_dist,periodic_map);

for it = 1:lloyd_iter

    grain_ID = Voronoi_From_Seeds(GRID,seed_xy,periodic_map);
    [grain_ID,seed_xy] = Remove_Empty_Grains(grain_ID,seed_xy);
    seed_xy = Grain_Centroids(grain_ID,seed_xy,GRID);

end

grain_ID = Voronoi_From_Seeds(GRID,seed_xy,periodic_map);
[grain_ID,seed_xy] = Remove_Empty_Grains(grain_ID,seed_xy);

[grain_ID,seed_xy] = Clean_Grain_Map_Topology( ...
    grain_ID,seed_xy,GRID,min_grain_pixels,island_clean_iter,quad_clean_iter,[]);

end


function seed_xy = Generate_Seeds_MinDist(GRID,Ngrain,min_seed_dist,periodic_map)

xmin = min(GRID.x);
ymin = min(GRID.y);
Lx   = max(GRID.x)-xmin;
Ly   = max(GRID.y)-ymin;

seed_xy = zeros(Ngrain,2);
nseed   = 0;
ntrial  = 0;
dmin    = min_seed_dist;

while nseed < Ngrain

    ntrial = ntrial + 1;
    cand   = [xmin+Lx*rand, ymin+Ly*rand];

    if nseed == 0

        ok = true;

    else

        dx = abs(seed_xy(1:nseed,1)-cand(1));
        dy = abs(seed_xy(1:nseed,2)-cand(2));

        if periodic_map == 1
            dx = min(dx,Lx-dx);
            dy = min(dy,Ly-dy);
        end

        ok = all(dx.^2 + dy.^2 >= dmin^2);

    end

    if ok
        nseed = nseed + 1;
        seed_xy(nseed,:) = cand;
    end

    if ntrial > 5000
        dmin   = 0.90*dmin;
        ntrial = 0;
    end
end

end


function grain_ID = Voronoi_From_Seeds(GRID,seed_xy,periodic_map)

[X,Y] = meshgrid(GRID.x,GRID.y);

xmin = min(GRID.x);
ymin = min(GRID.y);
Lx   = max(GRID.x)-xmin;
Ly   = max(GRID.y)-ymin;

Ngrain  = size(seed_xy,1);
dist2    = inf(GRID.ny,GRID.nx);
grain_ID = zeros(GRID.ny,GRID.nx);

for ig = 1:Ngrain

    ddx = abs(X-seed_xy(ig,1));
    ddy = abs(Y-seed_xy(ig,2));

    if periodic_map == 1
        ddx = min(ddx,Lx-ddx);
        ddy = min(ddy,Ly-ddy);
    end

    trial = ddx.^2 + ddy.^2;
    take  = trial < dist2;

    dist2(take)    = trial(take);
    grain_ID(take) = ig;
end

end


function seed_xy = Grain_Centroids(grain_ID,seed_xy,GRID)

[X,Y]  = meshgrid(GRID.x,GRID.y);
Ngrain = size(seed_xy,1);

for ig = 1:Ngrain

    mask = grain_ID == ig;

    if any(mask(:))
        seed_xy(ig,1) = mean(X(mask));
        seed_xy(ig,2) = mean(Y(mask));
    end
end

end


function [grain_ID,seed_xy,fixed_grain] = Clean_Grain_Map_Topology( ...
    grain_ID,seed_xy,GRID,min_grain_pixels,island_clean_iter,quad_clean_iter,fixed_grain)

if nargin < 7
    fixed_grain = [];
end

fixed_grain = fixed_grain(:).';

for it = 1:island_clean_iter
    grain_ID = Merge_Tiny_Grains(grain_ID,min_grain_pixels,fixed_grain);
    grain_ID = Remove_Isolated_Pixels(grain_ID,fixed_grain);
end

for it = 1:quad_clean_iter
    grain_ID = Break_Quad_Junctions(grain_ID,fixed_grain);
    grain_ID = Remove_Isolated_Pixels(grain_ID,fixed_grain);
end

[grain_ID,seed_xy,id_map] = Remove_Empty_Grains(grain_ID,seed_xy);

if ~isempty(fixed_grain)
    fixed_grain = id_map(fixed_grain);
    fixed_grain = fixed_grain(fixed_grain > 0);
end

seed_xy = Grain_Centroids(grain_ID,seed_xy,GRID);

Check_Grain_Map(grain_ID)

end


function grain_ID = Merge_Tiny_Grains(grain_ID,min_grain_pixels,fixed_grain)

Ngrain = max(grain_ID(:));
area   = accumarray(grain_ID(:),1,[Ngrain,1]);
small  = find(area > 0 & area < min_grain_pixels);
small  = setdiff(small,fixed_grain);

for ii = 1:numel(small)

    ig = small(ii);

    if ~any(grain_ID(:) == ig)
        continue
    end

    neigh = Neighbour_Grains(grain_ID,ig);

    if isempty(neigh)
        continue
    end

    area_now = accumarray(grain_ID(:),1,[max(grain_ID(:)),1]);
    [~,im]   = max(area_now(neigh));
    target   = neigh(im);

    grain_ID(grain_ID == ig) = target;

end

end


function neigh = Neighbour_Grains(grain_ID,ig)

neigh = [];

A = grain_ID(:,1:end-1);
B = grain_ID(:,2:end);

neigh = [neigh; B(A == ig); A(B == ig)];

A = grain_ID(1:end-1,:);
B = grain_ID(2:end,:);

neigh = [neigh; B(A == ig); A(B == ig)];
neigh = unique(neigh);
neigh = neigh(neigh > 0 & neigh ~= ig);

end


function grain_ID = Remove_Isolated_Pixels(grain_ID,fixed_grain)

[ny,nx] = size(grain_ID);

for iy = 2:ny-1
    for ix = 2:nx-1

        id = grain_ID(iy,ix);

        if ismember(id,fixed_grain)
            continue
        end

        nb = [grain_ID(iy,ix-1),grain_ID(iy,ix+1), ...
              grain_ID(iy-1,ix),grain_ID(iy+1,ix)];

        if ~any(nb == id)

            u   = unique(nb);
            cnt = zeros(size(u));

            for k = 1:numel(u)
                cnt(k) = nnz(nb == u(k));
            end

            [~,im] = max(cnt);
            grain_ID(iy,ix) = u(im);

        end
    end
end

end


function grain_ID = Break_Quad_Junctions(grain_ID,fixed_grain)

[ny,nx] = size(grain_ID);

for iy = 1:ny-1
    for ix = 1:nx-1

        block = grain_ID(iy:iy+1,ix:ix+1);
        u     = unique(block(:));

        if numel(u) < 4
            continue
        end

        if ~ismember(grain_ID(iy+1,ix+1),fixed_grain)
            grain_ID(iy+1,ix+1) = grain_ID(iy,ix+1);
        elseif ~ismember(grain_ID(iy,ix),fixed_grain)
            grain_ID(iy,ix) = grain_ID(iy+1,ix);
        end
    end
end

end


function n_quad = Count_Quad_Junctions(grain_ID)

[ny,nx] = size(grain_ID);
n_quad  = 0;

for iy = 1:ny-1
    for ix = 1:nx-1
        block = grain_ID(iy:iy+1,ix:ix+1);
        if numel(unique(block(:))) >= 4
            n_quad = n_quad + 1;
        end
    end
end

end


function Check_Grain_Map(grain_ID)

if any(grain_ID(:) < 1)
    error('Check_Grain_Map: some grid points have no grain.')
end

if any(~isfinite(grain_ID(:)))
    error('Check_Grain_Map: grain_ID contains non-finite values.')
end

end


function Check_Initial_Phi_Map(phi,grain_ID)

s = sum(phi,3);

if any(abs(s(:)-1) > 0)
    error('Check_Initial_Phi_Map: initial phi is not one-hot.')
end

[~,gid] = max(phi,[],3);

if any(gid(:) ~= grain_ID(:))
    error('Check_Initial_Phi_Map: phi and grain_ID are inconsistent.')
end

end


function [grain_ID,seed_xy,id_map] = Remove_Empty_Grains(grain_ID,seed_xy)

present = unique(grain_ID(:));
present = present(present > 0);

Nold   = size(seed_xy,1);
id_map = zeros(Nold,1);

if isempty(present)
    error('Remove_Empty_Grains: no grains remain.')
end

id_map(present) = 1:numel(present);

grain_ID = reshape(id_map(grain_ID(:)),size(grain_ID));
seed_xy  = seed_xy(present,:);

if numel(present) ~= Nold
    warning('%d empty grain seeds removed after rasterization.', ...
        Nold-numel(present))
end

end


function [grain_phase,phase_prop_real] = Assign_Grain_Phases_By_Area(grain_ID,phase_prop)

Ngrain = max(grain_ID(:));
Nphase = numel(phase_prop);

area   = accumarray(grain_ID(:),1,[Ngrain,1]);
Atot   = sum(area);
target = phase_prop(:).'*Atot;

grain_phase = zeros(Ngrain,1);
area_now    = zeros(1,Nphase);

[~,ord] = sort(area,'descend');

for kk = 1:Ngrain

    ig    = ord(kk);
    merit = zeros(1,Nphase);

    for ip = 1:Nphase
        trial     = area_now;
        trial(ip) = trial(ip) + area(ig);
        merit(ip) = sum((trial-target).^2);
    end

    [~,iph]         = min(merit);
    grain_phase(ig) = iph;
    area_now(iph)   = area_now(iph) + area(ig);
end

count = accumarray(grain_phase,1,[Nphase,1]).';

for ip = find(phase_prop > 0 & count == 0)

    donors = find(count > 1);

    if isempty(donors)
        error('Cannot allocate at least one grain to every requested phase.')
    end

    [~,jd] = max(area_now(donors)-target(donors));
    donor  = donors(jd);

    grains = find(grain_phase == donor);
    [~,jm] = min(area(grains));
    ig     = grains(jm);

    grain_phase(ig) = ip;

    count(donor)    = count(donor)-1;
    count(ip)       = count(ip)+1;

    area_now(donor) = area_now(donor)-area(ig);
    area_now(ip)    = area_now(ip)+area(ig);
end

phase_prop_real = area_now/Atot;

end


function [grain_phase,phase_prop_real] = Assign_Grain_Phases_By_Area_FixedSeed( ...
    grain_ID,phase_prop,fixed_grain,fixed_phase,forbidden_phase)

Ngrain = max(grain_ID(:));
Nphase = numel(phase_prop);

if nargin < 5 || isempty(forbidden_phase)
    forbidden_phase = [];
end

allowed_phase = setdiff(find(phase_prop > 0),forbidden_phase);

if isempty(allowed_phase)
    error('No allowed non-seed phases left for assignment.')
end

area   = accumarray(grain_ID(:),1,[Ngrain,1]);
Atot   = sum(area);
target = phase_prop(:).'*Atot;

grain_phase = zeros(Ngrain,1);
area_now    = zeros(1,Nphase);

fixed_grain = fixed_grain(:).';
fixed_phase = fixed_phase(:).';

for ii = 1:numel(fixed_grain)
    ig = fixed_grain(ii);
    ip = fixed_phase(ii);

    grain_phase(ig) = ip;
    area_now(ip)    = area_now(ip) + area(ig);
end

free = setdiff(1:Ngrain,fixed_grain);
[~,ord_local] = sort(area(free),'descend');
ord = free(ord_local);

for kk = 1:numel(ord)

    ig = ord(kk);

    merit = zeros(1,numel(allowed_phase));

    for ii = 1:numel(allowed_phase)
        ip = allowed_phase(ii);
        trial     = area_now;
        trial(ip) = trial(ip) + area(ig);
        merit(ii) = sum((trial-target).^2);
    end

    [~,ii_best]     = min(merit);
    iph             = allowed_phase(ii_best);
    grain_phase(ig) = iph;
    area_now(iph)   = area_now(iph) + area(ig);
end

phase_prop_real = area_now/Atot;

end


function [grain_phase,phase_prop_real] = Assign_Grain_Phases_By_Area_Separated_Strict( ...
    grain_ID,phase_prop,fixed_grain,fixed_phase,forbidden_phase,ntry,rng_seed)

if nargin < 3 || isempty(fixed_grain)
    fixed_grain = [];
end
if nargin < 4 || isempty(fixed_phase)
    fixed_phase = [];
end
if nargin < 5 || isempty(forbidden_phase)
    forbidden_phase = [];
end
if nargin < 6 || isempty(ntry)
    ntry = 50;
end
if nargin < 7 || isempty(rng_seed)
    rng_seed = 1;
end

Ngrain = max(grain_ID(:));
Nphase = numel(phase_prop);

allowed_phase = setdiff(find(phase_prop > 0),forbidden_phase);

if isempty(allowed_phase)
    error('No allowed non-seed phases left for assignment.')
end

area   = accumarray(grain_ID(:),1,[Ngrain,1]);
Atot   = sum(area);
target = phase_prop(:).'*Atot;
adj    = Build_Grain_Adjacency(grain_ID);

best_phase = [];
best_score = inf;
best_area  = [];

rng_state = rng;

for itry = 1:ntry

    rng(rng_seed+itry,'twister')

    grain_phase = zeros(Ngrain,1);
    area_now    = zeros(1,Nphase);
    count       = zeros(1,Nphase);

    for ii = 1:numel(fixed_grain)
        ig = fixed_grain(ii);
        ip = fixed_phase(ii);

        grain_phase(ig) = ip;
        area_now(ip)    = area_now(ip) + area(ig);
        count(ip)       = count(ip) + 1;
    end

    free = setdiff(1:Ngrain,fixed_grain);
    [~,ord_area] = sort(area(free),'descend');
    ord = free(ord_area);

    if mod(itry,2) == 0
        ord = ord(randperm(numel(ord)));
    end

    for kk = 1:numel(ord)

        ig    = ord(kk);
        neigh = find(adj(ig,:) & grain_phase(:).' > 0);

        missing = find(phase_prop > 0 & count == 0 & ...
                       ~ismember(1:Nphase,forbidden_phase));

        if numel(ord)-kk+1 <= numel(missing)
            candidates = missing;
        else
            candidates = allowed_phase;
        end

        bad_phase  = unique(grain_phase(neigh));
        candidates = setdiff(candidates,bad_phase);

        if isempty(candidates)
            candidates = allowed_phase;
        end

        merit = zeros(1,numel(candidates));

        for ii = 1:numel(candidates)

            ip = candidates(ii);

            trial     = area_now;
            trial(ip) = trial(ip) + area(ig);

            area_err   = sum(((trial-target)./max(Atot,1)).^2);
            same_neigh = nnz(grain_phase(neigh) == ip);

            merit(ii) = area_err + 1e6*same_neigh;

        end

        [~,im]          = min(merit);
        ip              = candidates(im);
        grain_phase(ig) = ip;
        area_now(ip)    = area_now(ip) + area(ig);
        count(ip)       = count(ip) + 1;

    end

    [grain_phase,area_now,count] = Repair_Same_Phase_Contacts( ...
        grain_phase,area_now,count,area,target,adj,allowed_phase,phase_prop,forbidden_phase);

    n_same   = Count_Same_Phase_Contacts(grain_ID,grain_phase);
    area_err = sum(((area_now-target)./max(Atot,1)).^2);
    score    = 1e6*n_same + area_err;

    if score < best_score
        best_score = score;
        best_phase = grain_phase;
        best_area  = area_now;
    end

    if n_same == 0
        break
    end

end

rng(rng_state)

grain_phase     = best_phase;
phase_prop_real = best_area/Atot;

n_same = Count_Same_Phase_Contacts(grain_ID,grain_phase);

fprintf('Same-phase grain contacts after strict assignment = %d\n',n_same)

if n_same > 0
    warning(['Strict zero same-phase contacts was not reached. ', ...
             'Increase Ngrain, reduce dominant phase fraction, or increase phase_assign_tries.'])
end

end


function [grain_phase,area_now,count] = Repair_Same_Phase_Contacts( ...
    grain_phase,area_now,count,area,target,adj,allowed_phase,phase_prop,forbidden_phase)

for iter = 1:40

    changed = 0;
    [i,j] = find(triu(adj,1));

    bad_pair = [];

    for k = 1:numel(i)
        if grain_phase(i(k)) == grain_phase(j(k))
            bad_pair(end+1,:) = [i(k),j(k)]; %#ok<AGROW>
        end
    end

    if isempty(bad_pair)
        break
    end

    for kk = 1:size(bad_pair,1)

        grains = bad_pair(kk,:);
        best_ig = 0;
        best_ip = 0;
        best_merit = inf;

        for jj = 1:2

            ig   = grains(jj);
            iph0 = grain_phase(ig);

            if ismember(iph0,forbidden_phase)
                continue
            end

            if count(iph0) <= 1 && phase_prop(iph0) > 0
                continue
            end

            neigh = find(adj(ig,:) & grain_phase(:).' > 0);

            for ip = allowed_phase

                if ip == iph0
                    continue
                end

                if any(grain_phase(neigh) == ip)
                    continue
                end

                trial       = area_now;
                trial(iph0) = trial(iph0) - area(ig);
                trial(ip)   = trial(ip)   + area(ig);

                merit = sum(((trial-target)./max(sum(area),1)).^2);

                if merit < best_merit
                    best_merit = merit;
                    best_ig    = ig;
                    best_ip    = ip;
                end
            end
        end

        if best_ig > 0

            iph0 = grain_phase(best_ig);

            grain_phase(best_ig) = best_ip;

            area_now(iph0)    = area_now(iph0) - area(best_ig);
            area_now(best_ip) = area_now(best_ip) + area(best_ig);

            count(iph0)    = count(iph0) - 1;
            count(best_ip) = count(best_ip) + 1;

            changed = 1;

        end
    end

    if changed == 0
        break
    end
end

end


function adj = Build_Grain_Adjacency(grain_ID)

Ngrain = max(grain_ID(:));
adj    = false(Ngrain,Ngrain);

A = grain_ID(:,1:end-1);
B = grain_ID(:,2:end);
mask = A ~= B & A > 0 & B > 0;
pair = [A(mask),B(mask)];

A = grain_ID(1:end-1,:);
B = grain_ID(2:end,:);
mask = A ~= B & A > 0 & B > 0;
pair = [pair; A(mask),B(mask)];

for ii = 1:size(pair,1)
    i = pair(ii,1);
    j = pair(ii,2);
    adj(i,j) = true;
    adj(j,i) = true;
end

end


function n_same = Count_Same_Phase_Contacts(grain_ID,grain_phase)

adj = Build_Grain_Adjacency(grain_ID);
[i,j] = find(triu(adj,1));

n_same = 0;

for k = 1:numel(i)
    if grain_phase(i(k)) == grain_phase(j(k))
        n_same = n_same + 1;
    end
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
