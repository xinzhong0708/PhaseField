%% Create_Map2D_Clean_PhaseProp.m
% Generate Map2d.mat with polygonal grains and phase proportions controlled
% by phase_prop. Phase separation is deliberately not used here because it
% can dominate the area target and over-assign small phases.

clear; clf

addpath('..\bin')
addpath('..\')
addpath('..\Thermo')
addpath('..\Thermo\Solutions')

%% User controls
T      = 500 + 273.15;
P      = 0.5e9;
Cname  = {'Fe' 'Mg' 'Ca' 'Al' 'Si' 'O'};
solmod = 'solution_models_PFM';

phs_name = {'Clinopyroxene','Orthopyroxene','Garnet','Kyanite','Quartz'};

phase_prop = [0.1 0.30 0.02 0.25 0.15];
phase_prop = phase_prop/sum(phase_prop);

c_value    = cell(1,numel(phs_name));
c_value{1} = [0.05 0.15 0 0.4];
c_value{2} = [0.2 0.2 0.1 0.1];
c_value{3} = [0.45 0.45];
c_value{4} = [1.0];
c_value{5} = [1.0];

Lx = 25e-6;
Ly = 25e-6;
nx = 120;
ny = 120;

grain_size           = 3.6e-6;
Ngrain_user          = [];
min_grain_per_phase  = 3;
rng_seed             = 5;
periodic_map         = 0;

use_center_seed  = 1;
seed_phase_name  = 'Garnet';
seed_radius      = [];        % []: derive radius from phase_prop
seed_center_real = [];        % []: domain center
single_seed_only = 1;

plot_map = 1;

PHYS        = struct();
PHYS.E_sc   = 1e9;
PHYS.L_sc   = 1;
PHYS.vref   = 2e-5;

E_sc  = PHYS.E_sc;
L_sc  = PHYS.L_sc;
vref  = PHYS.vref;
eta0  = 4000e10/E_sc;

%% Thermodynamics and reference LE
pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref);
Nphase     = numel(pars_phase);

c_guess = cell(1,Nphase);
for ip = 1:Nphase
    c_guess{ip} = num2cell(c_value{ip});
end

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

%% Grid
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

%% Uniform phase compositions from reference LE
c_phase = cell(1,Nphase);

for ip = 1:Nphase
    c_phase{ip} = cell(size(c_ref{ip}));

    for ic = 1:numel(c_ref{ip})
        c_phase{ip}{ic} = c_ref{ip}{ic}*ones(ny,nx);
    end
end

%% Polygonal grain map
phase_pos = phase_prop(phase_prop > 0);
Ngrain_min = ceil(min_grain_per_phase/min(phase_pos));
Ngrain_area = round((Lx/L_sc)*(Ly/L_sc)/(pi/4*(grain_size/L_sc)^2));

if isempty(Ngrain_user)
    Ngrain = max([Nphase,Ngrain_area,Ngrain_min]);
else
    Ngrain = Ngrain_user;
end

rng(rng_seed,'twister')
[grain_ID,seed_xy] = Generate_Polygonal_Grains_2D(GRID,Ngrain,periodic_map);
[grain_ID,seed_xy] = Remove_Empty_Grains(grain_ID,seed_xy);

fixed_grain    = [];
fixed_phase    = [];
forbidden_phase = [];

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
        seed_radius = sqrt(phase_prop(seed_phase)*Lx*Ly/pi);
    end

    seed_radius_sc = seed_radius/L_sc;
    seed_mask      = (X-seed_center(1)).^2 + ...
                     (Y-seed_center(2)).^2 <= seed_radius_sc^2;

    if ~any(seed_mask(:))
        error('Central seed contains no grid points. Increase seed_radius.')
    end

    fixed_grain_old = size(seed_xy,1) + 1;
    fixed_phase     = seed_phase;

    grain_ID(seed_mask) = fixed_grain_old;
    seed_xy(fixed_grain_old,:) = seed_center;

    [grain_ID,seed_xy,id_map] = Remove_Empty_Grains(grain_ID,seed_xy);
    fixed_grain = id_map(fixed_grain_old);

    if single_seed_only == 1
        forbidden_phase = fixed_phase;
    end
end

Ngrain = size(seed_xy,1);
Np     = Ngrain;

[grain_phase,phase_prop_geom] = Assign_Grain_Phases_By_Area_Target( ...
    grain_ID,phase_prop,fixed_grain,fixed_phase,forbidden_phase);

phi      = zeros(ny,nx,Ngrain);
phase_ID = zeros(ny,nx);

for ig = 1:Ngrain
    mask           = grain_ID == ig;
    phi(:,:,ig)    = mask;
    phase_ID(mask) = grain_phase(ig);
end

phase_prop_geom = Compute_Phase_Fraction_From_Map(phase_ID,Nphase);

theta_grain = pi*rand(1,Ngrain);

%% Grain-resolved MODEL
MODEL             = struct();
MODEL.phs_name    = phs_name;
MODEL.phase_index = grain_phase(:).';
MODEL.pars        = pars_phase(MODEL.phase_index);
MODEL.eta         = eta0;
MODEL.Cname       = Cname;
MODEL.solmod      = solmod;
MODEL.T           = T;
MODEL.P           = P;
MODEL.E_sc        = E_sc;
MODEL.vref        = vref;

EPS = 1e-14;
MODEL.p_fun = @(a,phi) ...
    phi(:,:,a).^2 ./ (sum(phi.^2,3) + EPS);

MODEL.dpdphi = @(a,b,phi) ...
    (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + EPS) ...
    - 2*phi(:,:,a).*phi(:,:,b).^2 ./ (sum(phi.^2,3) + EPS).^2;

pars = MODEL.pars;
F    = MODEL;

%% Initial STATE
c   = Expand_c_By_Phase(c_phase,MODEL.phase_index);
p   = Calc_p(MODEL,phi);
e   = Calc_e(pars,c);
E   = Calc_E_Tot(e,p);
eta = eta0*ones(ny,nx);

for ie = 1:Ne
    E{ie} = E{ie} + E_offset{ie};
end

E_bulk_shift = zeros(1,Ne);
mu_e         = repmat({zeros(ny,nx)},1,Ne);
chi          = repmat({zeros(ny,nx)},Ne,Ne);

PARAM            = struct();
PARAM.Np         = Np;
PARAM.Ne         = Ne;
PARAM.eta        = eta;
PARAM.use_WScale = 0;
PARAM.T          = T;
PARAM.P          = P;
PARAM.theta_grain = theta_grain;
PARAM.LE_mode    = 'LE';

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

STATE = LE_Run_Mode_New(STATE,PARAM,MODEL);

p_phase        = Collapse_p_By_Phase(STATE.p,MODEL.phase_index,Nphase);
phase_prop_map = zeros(1,Nphase);

for ip = 1:Nphase
    tmp = p_phase(:,:,ip);
    phase_prop_map(ip) = mean(tmp(:));
end

c    = STATE.c;
e    = STATE.e;
E    = STATE.E;
p    = STATE.p;
mu_e = STATE.mu_e;
chi  = STATE.chi;
phi  = STATE.phi;

%% Plot and save
if plot_map == 1
    figure(3); clf
    imagesc(x*1e6,y*1e6,phase_ID)
    set(gca,'YDir','normal')
    axis image tight
    colorbar
    title('Initial phase map')
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    drawnow
end

grain_size_real = sqrt(4*(GRID.Lx*GRID.Ly/Ngrain)/pi)*L_sc;
map_mode = 'polygon';
separate_same_phase = 0;
phase_sep_weight = 0;
theta_grain_user = [];
theta_rng_seed = [];
rng_seed = rng_seed;

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

%% Local helper functions
function pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref)

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

function [grain_ID,seed_xy,id_map] = Remove_Empty_Grains(grain_ID,seed_xy)

present = unique(grain_ID(:));
present = present(present > 0);

Nold = size(seed_xy,1);
id_map = zeros(Nold,1);
id_map(present) = 1:numel(present);

if numel(present) == Nold
    return
end

grain_ID = reshape(id_map(grain_ID(:)),size(grain_ID));
seed_xy  = seed_xy(present,:);

end

function [grain_phase,phase_prop_real] = Assign_Grain_Phases_By_Area_Target( ...
    grain_ID,phase_prop,fixed_grain,fixed_phase,forbidden_phase)

Ngrain = max(grain_ID(:));
Nphase = numel(phase_prop);

if nargin < 3 || isempty(fixed_grain)
    fixed_grain = [];
    fixed_phase = [];
end
if nargin < 5 || isempty(forbidden_phase)
    forbidden_phase = [];
end

area   = accumarray(grain_ID(:),1,[Ngrain,1]);
Atot   = sum(area);
target = phase_prop(:).'*Atot;

allowed_phase = setdiff(find(phase_prop > 0),forbidden_phase);

if isempty(allowed_phase)
    error('No allowed phase remains for free grains.')
end

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
[~,ord_local] = sort(area(free),'descend');
ord = free(ord_local);

for kk = 1:numel(ord)

    ig = ord(kk);

    missing = allowed_phase(count(allowed_phase) == 0);

    if numel(ord)-kk+1 <= numel(missing)
        candidates = missing;
    else
        candidates = allowed_phase;
    end

    merit = zeros(1,numel(candidates));

    for ii = 1:numel(candidates)
        ip = candidates(ii);
        trial     = area_now;
        trial(ip) = trial(ip) + area(ig);
        merit(ii) = sum((trial-target).^2);
    end

    [~,ibest] = min(merit);
    iph = candidates(ibest);

    grain_phase(ig) = iph;
    area_now(iph)   = area_now(iph) + area(ig);
    count(iph)      = count(iph) + 1;
end

for iter = 1:50

    changed = 0;

    for kk = 1:numel(free)

        ig   = free(kk);
        iph0 = grain_phase(ig);

        best_phase = iph0;
        best_merit = sum((area_now-target).^2);

        for ii = 1:numel(allowed_phase)

            iph = allowed_phase(ii);

            if iph == iph0
                continue
            end

            if count(iph0) <= 1 && phase_prop(iph0) > 0
                continue
            end

            trial       = area_now;
            trial(iph0) = trial(iph0) - area(ig);
            trial(iph)  = trial(iph)  + area(ig);
            merit       = sum((trial-target).^2);

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
