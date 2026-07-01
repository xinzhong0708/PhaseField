%% Create_Map2D_Clean_PT_Theta_Clean.m
% Generate Map2d.mat for either:
%   map_mode = 'bands'
%   map_mode = 'polygon'
%
% Full thermodynamic phase names are used.
% Thermodynamic pars are built on the fly from T and P.
% No Data_*.mat files are loaded.
%
% Polygon mode improvements:
%   1) seeds are generated with minimum spacing and Lloyd relaxation
%   2) tiny / isolated grains and 2x2 four-grain junctions are cleaned
%   3) neighbouring grains are assigned different thermodynamic phases when possible
%   4) each grid point belongs to exactly one grain in the initial map
%
% Grain orientation is saved as:
%   PARAM.theta_grain(ig)
%
% This is used later by:
%   Calc_AC_Anisotropy_FacetedStiffness

clear; clf

addpath('..\bin')
addpath('..\')
addpath('..\Thermo')
addpath('..\Thermo\Solutions')

%% ========================================================================
%  User controls
% ========================================================================

map_mode = 'polygon';     % 'bands' or 'polygon'

% Pressure-temperature initial condition
T              = 900 + 273.15;       % K
P              = 0.3e9;              % Pa

% Elements and solution model
% Use only elements that are present in this PFM run.
% Make sure c_value below is consistent with this choice.
Cname          = {'Fe' 'Mg' 'Ca' 'Al' 'Na' 'Si' 'O'};
% Cname        = {'Fe' 'Mg' 'Si' 'O'};
solmod         = 'solution_models_PFM';

% Thermodynamic phases used in this PFM run.
% Use full names accepted by init_thermo.
phs_name = {'Garnet','Clinopyroxene','Kyanite','Feldspar','Quartz'};

% Requested phase proportions
phase_prop = [0.01 0.3 0.1 0.25 0.15];
phase_prop = phase_prop/sum(phase_prop);

% Initial independent endmember compositions
c_value = cell(1,numel(phs_name));
c_value{1} = [0.45 0.45];                 % Garnet
c_value{2} = [0.05 0.15 0 0.04 0.4];      % Clinopyroxene
c_value{3} = [1.0];                       % Kyanite
c_value{4} = [0.5];                       % Feldspar
c_value{5} = [1.0];                       % Quartz

% Grain orientation controls
% [] means random orientation in [0,pi) for each grain.
% Otherwise use one value per grain. The length must be Ngrain.
theta_grain_user = [];
theta_rng_seed   = 1002;

% Domain
Lx = 25e-6;
Ly = 25e-6;

switch lower(map_mode)
    case 'bands'
        nx = 250;
        ny = 4;
    case 'polygon'
        nx = 90;
        ny = 90;
    otherwise
        error('Unknown map_mode: %s',map_mode)
end

% Polygon controls
grain_size   = 9.50e-6;   % mean equivalent-circle diameter, meter
Ngrain_user  = [];        % [] means estimate from grain_size
rng_seed     = 6;
periodic_map = 0;         % keep 0 unless the PF solver is periodic

% Clean polygon controls
seed_min_dist_fac  = 0.45;   % minimum seed spacing relative to grain_size
lloyd_iter         = 12;     % seed relaxation iterations
min_grain_pixels   = 8;      % merge grains smaller than this
island_clean_iter  = 6;      % remove isolated pixels / tiny defects
quad_clean_iter    = 6;      % remove 2x2 four-grain junctions
phase_assign_tries = 80;     % repeated phase-colouring attempts

% Initial diffuse-interface controls.
% This keeps polygonal grains, but avoids a one-hot interface at t = 0.
% The smoothing is grain-local and keeps only the two largest grain fields
% at each pixel, so LE does not start with artificial 4-5 phase overlap.
smooth_initial_phi = 1;
phi_smooth_iter    = 6;
phi_top_keep       = 2;

% Central seed controls, polygon mode only
use_center_seed  = 1;
seed_phase_name  = 'Garnet';
seed_radius      = 2.30e-6;   % meter
seed_center_real = [];        % [] means domain center, or [x y] in meter
single_seed_only = 1;         % do not assign this seed phase to any other grain

% Phase separation controls, polygon mode only
separate_same_phase = 1;       % avoid neighbouring grains with the same phase
phase_sep_weight    = 1000;    % kept for compatibility / save output

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
%  Build thermodynamics and reference LE
% ========================================================================

pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref);
Nphase     = numel(pars_phase);

c_guess = cell(1,Nphase);

for ip = 1:Nphase
    c_guess{ip} = num2cell(c_value{ip});
end

% Reference LE uses requested phase proportions.
% This same c_ref and E_offset are used for both bands and polygon.
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
%  Generate phase / grain map
% ========================================================================

switch lower(map_mode)

    case 'bands'

        Nrepeat    = 1;
        band_order = repmat(1:Nphase,1,Nrepeat);

        Noccur        = accumarray(band_order(:),1,[Nphase,1]).';
        band_fraction = phase_prop(band_order)./Noccur(band_order);

        grain_phase = band_order(:);
        Ngrain      = numel(grain_phase);
        Np          = Ngrain;

        x_edges      = 1 + round(nx*[0 cumsum(band_fraction)]);
        x_edges(1)   = 1;
        x_edges(end) = nx + 1;

        phi      = zeros(ny,nx,Ngrain);
        phase_ID = zeros(ny,nx);
        grain_ID = zeros(ny,nx);
        seed_xy  = zeros(Ngrain,2);

        for ig = 1:Ngrain

            iph = grain_phase(ig);
            ix1 = x_edges(ig);
            ix2 = x_edges(ig+1)-1;

            phi(:,ix1:ix2,ig)   = 1;
            phase_ID(:,ix1:ix2) = iph;
            grain_ID(:,ix1:ix2) = ig;

            seed_xy(ig,1) = mean(x(ix1:ix2));
            seed_xy(ig,2) = mean(y);
        end

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

        % Optional central seed grain. The rest of the domain remains a
        % normal polygonal grain map. The seed is treated as one additional
        % grain with fixed phase.
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

            seed_radius_sc = seed_radius/L_sc;
            seed_mask      = (X-seed_center(1)).^2 + ...
                             (Y-seed_center(2)).^2 <= seed_radius_sc^2;

            if ~any(seed_mask(:))
                error('Central seed contains no grid points. Increase seed_radius.')
            end

            fixed_grain = size(seed_xy,1) + 1;
            fixed_phase = seed_phase;

            grain_ID(seed_mask) = fixed_grain;
            seed_xy(fixed_grain,:) = seed_center;

            fprintf('Central %s seed: radius %.4e m, area fraction %.8f\n', ...
                seed_phase_name,seed_radius,nnz(seed_mask)/numel(seed_mask))
        end

        % Clean again after optional seed insertion and update fixed grain id.
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

        if same_contact > 0
            warning(['Some same-phase grain contacts remain. ', ...
                'Increase Ngrain_user, reduce dominant phase fraction, or increase phase_assign_tries.'])
        end

        phi      = zeros(ny,nx,Ngrain);
        phase_ID = zeros(ny,nx);

        for ig = 1:Ngrain
            mask           = grain_ID == ig;
            phi(:,:,ig)    = mask;
            phase_ID(mask) = grain_phase(ig);
        end
end

% ------------------------------------------------------------------------
% Smooth initial polygonal phi, if requested.
% Bands are left unchanged because the pseudo-1D case already relaxes well.
% ------------------------------------------------------------------------
if strcmpi(map_mode,'polygon') && smooth_initial_phi == 1

    phi = Smooth_Phi_TopN(phi,phi_smooth_iter,phi_top_keep);
    Check_Initial_Phi_Map_Soft(phi,grain_ID);

else

    Check_Initial_Phi_Map(phi,grain_ID)

end

% ------------------------------------------------------------------------
% Grain orientations
% ------------------------------------------------------------------------
if isempty(theta_grain_user)

    rng_state = rng;

    if ~isempty(theta_rng_seed)
        rng(theta_rng_seed,'twister')
    end

    % Crystallographic plane normals are unoriented, so [0,pi) is enough.
    theta_grain = pi*rand(1,Ngrain);

    rng(rng_state)

else

    if numel(theta_grain_user) ~= Ngrain
        error('theta_grain_user must be [] or length Ngrain.')
    end

    theta_grain = theta_grain_user(:).';

end

phase_prop_geom = Compute_Phase_Fraction_From_Map(phase_ID,Nphase);

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

if strcmpi(map_mode,'bands')
    title('Initial pseudo-1D phase bands')
else
    title('Initial cleaned polygonal phase map')
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

if strcmpi(map_mode,'polygon')
    grain_size_real = sqrt(4*(GRID.Lx*GRID.Ly/Ngrain)/pi)*L_sc;
    fprintf('Requested grain diameter = %.4e m\n',grain_size)
    fprintf('Realised mean equivalent diameter = %.4e m\n',grain_size_real)
else
    grain_size_real = NaN;
end

fprintf('\nPhase       requested        realised area fraction      difference\n')
for ip = 1:Nphase
    fprintf('%-18s %.8f         %.8f                 %+.3e\n', ...
        phs_name{ip},phase_prop(ip),phase_prop_geom(ip), ...
        phase_prop_geom(ip)-phase_prop(ip))
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

% Local penalty consistency:
% E - Emix = mu_ref / eta.
%
% Important:
% Do not add any polygon bulk correction here.
% A small domain-average E mismatch is acceptable.
% The important thing is local consistency for LE_Calculator.
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

% Mild LE defaults saved with the map. These only matter if the runtime
% metadata reader does not overwrite them. They help the first PF step when
% polygonal interfaces contain newly activated phases.
PARAM.LE_p_tail   = 1e-5;
PARAM.LE_p_full   = 1e-3;
PARAM.LE_p_on     = 5e-4;
PARAM.LE_p_off    = 2e-4;
PARAM.LE_Pmax     = 4;
PARAM.LE_alpha_LE = [0.3 0.2 0.1 0.03];
PARAM.LE_iter_LE  = [200 300 500 800];

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

if exist('LE_Run_Mode','file') == 2
    STATE = LE_Run_Mode(STATE,PARAM,MODEL);
else
    STATE = LE_Run_Mode_New(STATE,PARAM,MODEL);
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

E_mean        = zeros(Ne,1);
E_target_vec  = zeros(Ne,1);

for ie = 1:Ne
    E_mean(ie)       = mean(STATE.E{ie}(:));
    E_target_vec(ie) = E_target{ie};
end

mu_jump = zeros(1,Ne);

for ie = 1:Ne
    mu_jump(ie) = max(STATE.mu_e{ie}(:)) - min(STATE.mu_e{ie}(:));
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
    'seed_min_dist_fac','lloyd_iter','min_grain_pixels','island_clean_iter','quad_clean_iter','phase_assign_tries', ...
    'smooth_initial_phi','phi_smooth_iter','phi_top_keep', ...
    'use_center_seed','seed_phase_name','seed_radius','seed_center_real','single_seed_only', ...
    'separate_same_phase','phase_sep_weight', ...
    'E_target','E_offset','E_bulk_shift','c_ref','mu_ref','omega_ref', ...
    'c','e','E','p','mu_e','chi','phi')

fprintf('\nSaved Map2d.mat\n')
fprintf('Saved full-name phases:\n')
disp(MODEL.phs_name)
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

        % Break a 4-grain grid vertex into two 3-grain vertices.
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


function phi = Smooth_Phi_TopN(phi,niter,top_keep)
%SMOOTH_PHI_TOPN Create a diffuse but sparse initial grain field.
%
% The one-hot polygon map is smoothed by a compact 3x3 kernel. After each
% smoothing pass only the largest top_keep grain fields are kept at every
% pixel. This avoids artificial many-phase overlap at polygon junctions.

if nargin < 2 || isempty(niter)
    niter = 4;
end

if nargin < 3 || isempty(top_keep)
    top_keep = 2;
end

[ny,nx,Ngrain] = size(phi);

kernel = [1 2 1; 2 4 2; 1 2 1]/16;

for it = 1:niter

    phi_new = zeros(size(phi));

    for ig = 1:Ngrain
        phi_new(:,:,ig) = conv2(phi(:,:,ig),kernel,'same');
    end

    if top_keep < Ngrain
        ps = sort(phi_new,3,'descend');
        cut = ps(:,:,top_keep);
        keep = phi_new > 0 & phi_new >= repmat(cut,1,1,Ngrain);
        phi_new(~keep) = 0;
    end

    phi = Normalize_Phi(phi_new);

end

% Remove tiny roundoff tails once more.
phi(phi < 1e-12) = 0;
phi = Normalize_Phi(phi);

% Avoid exact zero/one only at interface is not necessary. Pure interiors can
% remain exactly one-hot.

end


function phi = Normalize_Phi(phi)

s = sum(phi,3);
bad = s <= 0;

if any(bad(:))
    % This should not happen for a valid one-hot starting map.
    error('Normalize_Phi: zero sum phi encountered.')
end

for ig = 1:size(phi,3)
    phi(:,:,ig) = phi(:,:,ig)./s;
end

end


function Check_Initial_Phi_Map_Soft(phi,grain_ID)

s = sum(phi,3);
err = max(abs(s(:)-1));

if err > 1e-10
    error('Check_Initial_Phi_Map_Soft: sum(phi) differs from 1 by %.3e.',err)
end

[~,gid] = max(phi,[],3);
mismatch = nnz(gid(:) ~= grain_ID(:))/numel(grain_ID);

if mismatch > 0.05
    warning('Check_Initial_Phi_Map_Soft: %.2f%% pixels changed dominant grain after smoothing.',100*mismatch)
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

        % Strictly remove phases already present in assigned neighbours.
        bad_phase  = unique(grain_phase(neigh));
        candidates = setdiff(candidates,bad_phase);

        % If impossible locally, fall back to allowed phases and penalize.
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
%BUILD_GRAIN_ADJACENCY Four-neighbour grain adjacency graph.

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
%COUNT_SAME_PHASE_CONTACTS Count unique neighbouring same-phase grain pairs.

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
