%% Create_Map2D_Clean_PT_Theta.m
% Generate Map2d.mat for either:
%   map_mode = 'bands'
%   map_mode = 'polygon'
%
% Full thermodynamic phase names are used.
% Thermodynamic pars are built on the fly from T and P.
% No Data_*.mat files are loaded.
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
T              = 400 + 273.15;       % K
P              = 0.4e9;              % Pa

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
phase_prop = [0.001 0.3 0.1 0.2 0.1];
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
Lx = 5e-6;
Ly = 5e-6;

switch lower(map_mode)
    case 'bands'
        nx = 250;
        ny = 4;
    case 'polygon'
        nx = 82;
        ny = 82;
    otherwise
        error('Unknown map_mode: %s',map_mode)
end

% Polygon controls
grain_size   = 1.5e-6;    % mean equivalent-circle diameter, meter
Ngrain_user  = [];        % [] means estimate from grain_size
rng_seed     = 2;
periodic_map = 0;         % keep 0 unless the PF solver is periodic

% Central seed controls, polygon mode only
use_center_seed  = 1;
seed_phase_name  = 'Garnet';
seed_radius      = 0.30e-6;   % meter
seed_center_real = [];        % [] means domain center, or [x y] in meter
single_seed_only = 1;         % do not assign this seed phase to any other grain

% Phase separation controls, polygon mode only
separate_same_phase = 1;       % avoid neighbouring grains with the same phase
phase_sep_weight    = 1000;    % large value makes separation dominate area fit

% Scaling / penalty
PHYS        = struct();
PHYS.E_sc   = 1e9;
PHYS.L_sc   = 1e-6;
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

        [grain_ID,seed_xy] = Generate_Polygonal_Grains_2D( ...
            GRID,Ngrain,periodic_map);

        [grain_ID,seed_xy] = Remove_Empty_Grains(grain_ID,seed_xy);

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

            [grain_phase,~] = Assign_Grain_Phases_By_Area_Separated( ...
                grain_ID,phase_prop_assign,fixed_grain,fixed_phase,phase_sep_weight,forbidden_phase);

            same_contact = Count_Same_Phase_Contacts(grain_ID,grain_phase);

            fprintf('Same-phase grain contacts after assignment = %d\n',same_contact)

            if same_contact > 0
                warning(['Some same-phase grain contacts remain. ', ...
                    'Use smaller grain_size or larger Ngrain_user if strict separation is required.'])
            end

        elseif use_center_seed == 1

            [grain_phase,~] = Assign_Grain_Phases_By_Area_FixedSeed( ...
                grain_ID,phase_prop_assign,fixed_grain,fixed_phase,forbidden_phase);

        else

            [grain_phase,~] = Assign_Grain_Phases_By_Area(grain_ID,phase_prop);

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
    title('Initial random polygonal phase map')
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

PARAM            = struct();
PARAM.Np         = Np;
PARAM.Ne         = Ne;
PARAM.eta        = eta;
PARAM.use_WScale = 0;
PARAM.T          = T;
PARAM.P          = P;
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

E_mean       = zeros(Ne,1);
E_target_vec = zeros(Ne,1);

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


function [grain_ID,seed_xy] = Generate_Polygonal_Grains_2D(GRID,Ngrain,periodic_map)

[X,Y] = meshgrid(GRID.x,GRID.y);

xmin = min(GRID.x);
ymin = min(GRID.y);
Lx   = max(GRID.x)-xmin;
Ly   = max(GRID.y)-ymin;

seed_xy = [xmin+Lx*rand(Ngrain,1), ymin+Ly*rand(Ngrain,1)];

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


function [grain_ID,seed_xy] = Remove_Empty_Grains(grain_ID,seed_xy)

present = unique(grain_ID(:));
present = present(present > 0);

if numel(present) == size(seed_xy,1)
    return
end

Nold         = size(seed_xy,1);
map          = zeros(Nold,1);
map(present) = 1:numel(present);

grain_ID = reshape(map(grain_ID(:)),size(grain_ID));
seed_xy  = seed_xy(present,:);

warning('%d empty grain seeds removed after rasterization.', ...
    Nold-numel(present))

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


function [grain_phase,phase_prop_real] = Assign_Grain_Phases_By_Area_FixedSeed(grain_ID,phase_prop,fixed_grain,fixed_phase,forbidden_phase)

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

    if area(ig) == 0
        [~,ii] = min(area_now(allowed_phase)-target(allowed_phase));
        iph = allowed_phase(ii);
        grain_phase(ig) = iph;
        continue
    end

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

count = accumarray(grain_phase,1,[Nphase,1]).';

for ip = find(phase_prop > 0 & count == 0 & ~ismember(1:Nphase,forbidden_phase))

    donors = find(count > 1 & ~ismember(1:Nphase,forbidden_phase));

    if isempty(donors)
        error('Cannot allocate at least one grain to every requested phase.')
    end

    movable = setdiff(find(ismember(grain_phase,donors)),fixed_grain);

    if isempty(movable)
        error('Cannot allocate phase %d without changing the fixed seed.',ip)
    end

    [~,jm] = min(area(movable));
    ig     = movable(jm);
    donor  = grain_phase(ig);

    grain_phase(ig) = ip;

    count(donor)    = count(donor)-1;
    count(ip)       = count(ip)+1;

    area_now(donor) = area_now(donor)-area(ig);
    area_now(ip)    = area_now(ip)+area(ig);
end

phase_prop_real = area_now/Atot;

end


function [grain_phase,phase_prop_real] = Assign_Grain_Phases_By_Area_Separated( ...
    grain_ID,phase_prop,fixed_grain,fixed_phase,phase_sep_weight,forbidden_phase)
%ASSIGN_GRAIN_PHASES_BY_AREA_SEPARATED Assign phases while avoiding neighbours.
%
% The area target is still used, but a large penalty is added if a grain is
% assigned the same phase as an already assigned neighbouring grain.

if nargin < 3 || isempty(fixed_grain)
    fixed_grain = [];
end
if nargin < 4 || isempty(fixed_phase)
    fixed_phase = [];
end
if nargin < 5 || isempty(phase_sep_weight)
    phase_sep_weight = 1000;
end
if nargin < 6 || isempty(forbidden_phase)
    forbidden_phase = [];
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

grain_phase = zeros(Ngrain,1);
area_now    = zeros(1,Nphase);
count       = zeros(1,Nphase);

fixed_grain = fixed_grain(:).';
fixed_phase = fixed_phase(:).';

for ii = 1:numel(fixed_grain)
    ig = fixed_grain(ii);
    ip = fixed_phase(ii);

    grain_phase(ig) = ip;
    area_now(ip)    = area_now(ip) + area(ig);
    count(ip)       = count(ip) + 1;
end

free = setdiff(1:Ngrain,fixed_grain);
[~,ord_local] = sort(area(free),'descend');
ord = free(ord_local);

for kk = 1:numel(ord)

    ig    = ord(kk);
    neigh = find(adj(ig,:) & grain_phase(:).' > 0);

    missing = find(phase_prop > 0 & count == 0 & ~ismember(1:Nphase,forbidden_phase));

    if numel(ord)-kk+1 <= numel(missing)
        candidates = missing;
    else
        candidates = allowed_phase;
    end

    candidates_sep = [];

    for ii = 1:numel(candidates)
        ip = candidates(ii);

        if ~any(grain_phase(neigh) == ip)
            candidates_sep(end+1) = ip; %#ok<AGROW>
        end
    end

    if ~isempty(candidates_sep)
        candidates = candidates_sep;
    end

    merit = zeros(1,numel(candidates));

    for ii = 1:numel(candidates)

        ip = candidates(ii);

        trial     = area_now;
        trial(ip) = trial(ip) + area(ig);

        area_err   = sum(((trial-target)./max(Atot,1)).^2);
        same_neigh = nnz(grain_phase(neigh) == ip);

        merit(ii) = area_err + phase_sep_weight*same_neigh;

    end

    [~,ii_best]     = min(merit);
    iph             = candidates(ii_best);
    grain_phase(ig) = iph;
    area_now(iph)   = area_now(iph) + area(ig);
    count(iph)      = count(iph) + 1;

end

for iter = 1:20

    changed = 0;

    for kk = 1:numel(ord)

        ig    = ord(kk);
        iph0  = grain_phase(ig);
        neigh = find(adj(ig,:) & grain_phase(:).' > 0);

        if ~any(grain_phase(neigh) == iph0)
            continue
        end

        best_phase = iph0;
        best_merit = inf;

        for ii_phase = 1:numel(allowed_phase)

            iph = allowed_phase(ii_phase);

            if iph == iph0
                continue
            end

            if count(iph0) <= 1 && phase_prop(iph0) > 0
                continue
            end

            if any(grain_phase(neigh) == iph)
                continue
            end

            trial       = area_now;
            trial(iph0) = trial(iph0) - area(ig);
            trial(iph)  = trial(iph)  + area(ig);

            merit = sum(((trial-target)./max(Atot,1)).^2);

            if merit < best_merit
                best_merit = merit;
                best_phase = iph;
            end
        end

        if best_phase ~= iph0
            grain_phase(ig) = best_phase;
            area_now(iph0)  = area_now(iph0) - area(ig);
            area_now(best_phase) = area_now(best_phase) + area(ig);
            count(iph0) = count(iph0) - 1;
            count(best_phase) = count(best_phase) + 1;
            changed = 1;
        end
    end

    if changed == 0
        break
    end
end

phase_prop_real = area_now/Atot;

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
