%% Create_Map2D_Clean.m
% Generate Map2d.mat for either:
%   map_mode = 'bands'
%   map_mode = 'polygon'

clear; clf

addpath('..\bin')
addpath('..\')
addpath('..\ThermoData')

%% ========================================================================
%  User controls
% ========================================================================

map_mode = 'polygon';     % 'bands' or 'polygon'

% Thermodynamic phases
phs_name = {'Grt','Cpx','Fel','Olv'};

% Requested phase proportions
phase_prop = [0.30 0.20 0.20 0.30];
phase_prop = phase_prop/sum(phase_prop);

% Initial independent endmember compositions
c_value = cell(1,numel(phs_name));
c_value{1} = [0.45 0.45];                 % Grt
c_value{2} = [0.04 0.38 0 0.04 0.4];      % Cpx
c_value{3} = [0.5];                       % Fel
c_value{4} = [0.20 0.25 0.35];            % Olv

% Domain
Lx = 5e-6;
Ly = 5e-6;

switch lower(map_mode)
    case 'bands'
        nx = 200;
        ny = 5;

    case 'polygon'
        nx = 70;
        ny = 70;

    otherwise
        error('Unknown map_mode: %s',map_mode)
end

% Polygon controls
grain_size   = 5.0e-6;    % mean equivalent-circle diameter, meter
Ngrain_user  = [];        % [] means estimate from grain_size
rng_seed     = 5;
periodic_map = 0;         % keep 0 unless the PF solver is periodic

% Scaling / penalty
PHYS        = struct();
PHYS.E_sc   = 1e9;
PHYS.L_sc   = 1;

E_sc        = PHYS.E_sc;
L_sc        = PHYS.L_sc;
eta0        = 5000e10/E_sc;

%% ========================================================================
%  Load thermodynamics and reference LE
% ========================================================================

pars_phase = Load_Data(phs_name);
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

[c_ref,mu_ref] = LE_Calculator( ...
    pars_phase,p_ref,c_guess,E_target,eta0,[0.02,2000]);

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

        Ngrain = size(seed_xy,1);
        Np     = Ngrain;

        [grain_phase,~] = Assign_Grain_Phases_By_Area(grain_ID,phase_prop);

        phi      = zeros(ny,nx,Ngrain);
        phase_ID = zeros(ny,nx);

        for ig = 1:Ngrain
            mask           = grain_ID == ig;
            phi(:,:,ig)    = mask;
            phase_ID(mask) = grain_phase(ig);
        end
end

phase_prop_geom = Compute_Phase_Fraction_From_Map(phase_ID,Nphase);

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
%  Initial LE check
% ========================================================================

STATE_INI = STATE;

PARAM.LE_mode = 'LE';
STATE = LE_Run_Mode(STATE,PARAM,MODEL);

p_phase        = Collapse_p_By_Phase(STATE.p,MODEL.phase_index,Nphase);
phase_prop_map = zeros(1,Nphase);

for ip = 1:Nphase
    tmp = p_phase(:,:,ip);
    phase_prop_map(ip) = mean(tmp(:));
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

%% ========================================================================
%  Save
% ========================================================================

save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'E_sc','L_sc','eta','pars','Np','Ne','Nphase','Ngrain', ...
    'phs_name','phase_prop','phase_prop_geom','phase_prop_map', ...
    'map_mode','grain_ID','phase_ID','grain_phase','seed_xy', ...
    'grain_size','grain_size_real','Ngrain_user','rng_seed','periodic_map', ...
    'E_target','E_offset','E_bulk_shift','c_ref','mu_ref','omega_ref', ...
    'c','e','E','p','mu_e','chi','phi')

fprintf('\nSaved Map2d.mat\n')

%% ========================================================================
%  Local helper functions
% ========================================================================

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

% Ensure every requested phase has at least one grain.
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