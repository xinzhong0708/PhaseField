%% Create_Map2D_RockTexture_InitFixed.m
% Rock-like initial map with phase proportions and stable LE initialization.
%
% Important change:
%   The phase map is generated first. The reference LE state is then built
%   from the realised phase fractions of that exact map. The initial c, E,
%   mu_e and omg are filled from this reference state, so the first solver
%   LE call should not need to change c when p is unchanged.

clear; clf

addpath('..\bin')
addpath('..\')
addpath('..\Thermo')
addpath('..\Thermo\Solutions')

%% User controls
T      = 450 + 273.15;
P      = 0.3e9;
Cname  = {'Fe' 'Mg' 'Ca' 'Al' 'Si' 'O'};
solmod = 'solution_models_PFM';

phs_name = {'Clinopyroxene','Garnet','Kyanite','Quartz'};

phase_prop = [0.45 0.007 0.25 0.15];
phase_prop = phase_prop/sum(phase_prop);

c_value    = cell(1,numel(phs_name));
c_value{1} = [0.05 0.15 0 0.4];
c_value{2} = [0.45 0.45];
c_value{3} = [1.0];
c_value{4} = [1.0];

Lx = 25e-6;
Ly = 25e-6;
nx = 80;
ny = 80;

% Larger grain_size gives fewer cleaner grains.
grain_size   = 15e-6;
Ngrain_user  = [];
rng_seed     = 6;
periodic_map = 0;
lloyd_iter   = 8;

% Phases below this fraction are made as isolated inclusion grains.
small_phase_th        = 0.04;
small_phase_ngrain    = ones(1,numel(phs_name));
small_phase_ngrain(1) = 1;

% Optional fixed central seed.
use_center_seed  = 1;
seed_phase_name  = 'Garnet';
seed_radius      = [];       % []: use phase_prop area
seed_center_real = [];       % []: domain center
single_seed_only = 1;

% Assignment weights.
area_weight       = 1;
overshoot_weight  = 4;
same_phase_weight = 25;
assign_iter       = 80;

plot_map = 1;

theta_grain_user = [];
theta_rng_seed   = 1002;

PHYS      = struct();
PHYS.E_sc = 1e9;
PHYS.L_sc = 1;
PHYS.vref = 2e-5;

E_sc = PHYS.E_sc;
L_sc = PHYS.L_sc;
vref = PHYS.vref;
eta0 = 6000e10/E_sc;

%% Grid
x  = linspace(0,Lx,nx)/L_sc;
y  = linspace(0,Ly,ny)/L_sc;
dx = x(2)-x(1);
dy = y(2)-y(1);

GRID    = struct();
GRID.x  = x;
GRID.y  = y;
GRID.dx = dx;
GRID.dy = dy;
GRID.nx = nx;
GRID.ny = ny;
GRID.Lx = Lx/L_sc;
GRID.Ly = Ly/L_sc;

%% Rock-like grain and phase map first
rng(rng_seed,'twister')

if isempty(seed_radius)
    seed_radius_sc = [];
else
    seed_radius_sc = seed_radius/L_sc;
end

if isempty(seed_center_real)
    seed_center_sc = [];
else
    seed_center_sc = seed_center_real/L_sc;
end

[grain_ID,seed_xy,grain_phase] = Generate_Rock_Texture_Map( ...
    GRID,phs_name,phase_prop,grain_size/L_sc,Ngrain_user,periodic_map, ...
    lloyd_iter,small_phase_th,small_phase_ngrain, ...
    use_center_seed,seed_phase_name,seed_radius_sc,seed_center_sc, ...
    single_seed_only,area_weight,overshoot_weight,same_phase_weight,assign_iter);

Nphase = numel(phs_name);
Ngrain = numel(grain_phase);
Np     = Ngrain;

phi      = zeros(ny,nx,Ngrain);
phase_ID = zeros(ny,nx);

for ig = 1:Ngrain
    mask         = grain_ID == ig;
    phi(:,:,ig)  = mask;
    phase_ID(mask) = grain_phase(ig);
end

phase_prop_geom = Compute_Phase_Fraction_From_Map(phase_ID,Nphase);
phase_prop_ref  = phase_prop_geom;

%% Thermodynamics and reference LE for the realised map
pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref);

c_guess = cell(1,Nphase);
for ip = 1:Nphase
    c_guess{ip} = num2cell(c_value{ip});
end

p_ref    = reshape(phase_prop_ref,1,1,Nphase);
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

%% Uniform phase compositions from reference LE
c_phase = cell(1,Nphase);

for ip = 1:Nphase
    c_phase{ip} = cell(size(c_ref{ip}));
    for ic = 1:numel(c_ref{ip})
        c_phase{ip}{ic} = c_ref{ip}{ic}*ones(ny,nx);
    end
end

%% Grain orientations
if isempty(theta_grain_user)
    rng_state = rng;
    if ~isempty(theta_rng_seed)
        rng(theta_rng_seed,'twister')
    end
    theta_grain = pi*rand(1,Ngrain);
    rng(rng_state)
else
    if numel(theta_grain_user) ~= Ngrain
        error('theta_grain_user must be [] or length Ngrain.')
    end
    theta_grain = theta_grain_user(:).';
end

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

eps_phi = 1e-14;

MODEL.p_fun = @(a,phi) ...
    phi(:,:,a).^2 ./ (sum(phi.^2,3) + eps_phi);

MODEL.dpdphi = @(a,b,phi) ...
    (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) ...
    - 2*phi(:,:,a).*phi(:,:,b).^2 ./ (sum(phi.^2,3) + eps_phi).^2;

pars = MODEL.pars;
F    = MODEL;

%% Initial conserved field and thermodynamic fields
c   = Expand_c_By_Phase(c_phase,MODEL.phase_index);
p   = Calc_p(MODEL,phi);
e   = Calc_e(pars,c);
E   = Calc_E_Tot(e,p);
eta = eta0*ones(ny,nx);

for ie = 1:Ne
    E{ie} = E{ie} + E_offset{ie};
end

E_bulk_shift = zeros(1,Ne);

mu_e = repmat({zeros(ny,nx)},1,Ne);
for ie = 1:Ne
    mu_e{ie}(:) = mu_ref{ie};
end

chi = repmat({zeros(ny,nx)},Ne,Ne);
chi = Fill_Initial_Chi(chi,pars_phase,c_ref,grain_ID,grain_phase,Ne);

omg = zeros(ny,nx,Np);
for ig = 1:Ngrain
    iph = MODEL.phase_index(ig);
    omg(:,:,ig) = omega_ref(iph);
end

%% PARAM and STATE
PARAM             = struct();
PARAM.Np          = Np;
PARAM.Ne          = Ne;
PARAM.eta         = eta;
PARAM.use_WScale  = 0;
PARAM.T           = T;
PARAM.P           = P;
PARAM.theta_grain = theta_grain;
PARAM.LE_mode     = 'LE';

STATE          = struct();
STATE.c        = c;
STATE.e        = e;
STATE.E        = E;
STATE.mu_e     = mu_e;
STATE.chi      = chi;
STATE.omg      = omg;
STATE.phi      = phi;
STATE.p        = p;
STATE.mask     = ones(ny,nx,Np);
STATE.LE_state = [];

phase_prop_map = phase_prop_geom;

%% Plot
if plot_map == 1
    Plot_Rock_Map(x,y,L_sc,phase_ID,grain_ID)
end

grain_size_real = sqrt(4*(GRID.Lx*GRID.Ly/Ngrain)/pi)*L_sc;
map_mode = 'rock_texture';

%% Save
save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'E_sc','L_sc','vref','eta','pars','Np','Ne','Nphase','Ngrain', ...
    'T','P','Cname','solmod', ...
    'phs_name','phase_prop','phase_prop_geom','phase_prop_map','phase_prop_ref', ...
    'map_mode','grain_ID','phase_ID','grain_phase','seed_xy','theta_grain', ...
    'theta_grain_user','theta_rng_seed', ...
    'grain_size','grain_size_real','Ngrain_user','rng_seed','periodic_map','lloyd_iter', ...
    'small_phase_th','small_phase_ngrain', ...
    'use_center_seed','seed_phase_name','seed_radius','seed_center_real','single_seed_only', ...
    'area_weight','overshoot_weight','same_phase_weight','assign_iter', ...
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

function [grain_ID,seed_xy,grain_phase] = Generate_Rock_Texture_Map( ...
    GRID,phs_name,phase_prop,grain_size,Ngrain_user,periodic_map, ...
    lloyd_iter,small_phase_th,small_phase_ngrain, ...
    use_center_seed,seed_phase_name,seed_radius,seed_center_real, ...
    single_seed_only,area_weight,overshoot_weight,same_phase_weight,assign_iter)

Nphase = numel(phs_name);
Npix   = GRID.nx*GRID.ny;
[X,Y]  = meshgrid(GRID.x,GRID.y);

phase_target_pix = Phase_Target_Pixels(phase_prop,Npix);
phase_now_pix    = zeros(1,Nphase);
grain_ID         = zeros(GRID.ny,GRID.nx);
seed_xy          = zeros(0,2);
grain_phase      = zeros(0,1);

fixed_phase = [];

if use_center_seed == 1

    fixed_phase = find(strcmp(phs_name,seed_phase_name),1,'first');

    if isempty(fixed_phase)
        error('seed_phase_name %s is not in phs_name.',seed_phase_name)
    end

    if isempty(seed_center_real)
        seed_center = [mean(GRID.x),mean(GRID.y)];
    else
        seed_center = seed_center_real;
    end

    if isempty(seed_radius)
        npix_seed = phase_target_pix(fixed_phase);
        seed_mask = Closest_Pixels_Mask(X,Y,grain_ID == 0,seed_center,npix_seed);
    else
        seed_mask = (X-seed_center(1)).^2 + ...
                    (Y-seed_center(2)).^2 <= seed_radius^2;
        seed_mask = seed_mask & grain_ID == 0;
    end

    if any(seed_mask(:))
        [grain_ID,seed_xy,grain_phase,phase_now_pix] = Add_Grain( ...
            grain_ID,seed_xy,grain_phase,phase_now_pix,seed_mask,seed_center,fixed_phase);
    end
end

small_phase = find(phase_prop < small_phase_th & phase_prop > 0);

if single_seed_only == 1 && ~isempty(fixed_phase)
    small_phase = setdiff(small_phase,fixed_phase);
end

for ii = 1:numel(small_phase)

    ip = small_phase(ii);
    ng = max(1,small_phase_ngrain(ip));

    for jj = 1:ng

        npix_left = phase_target_pix(ip)-phase_now_pix(ip);
        npix_inc  = max(1,round(npix_left/(ng-jj+1)));
        center    = Random_Free_Center(GRID,grain_ID == 0);
        inc_mask  = Closest_Pixels_Mask(X,Y,grain_ID == 0,center,npix_inc);

        if any(inc_mask(:))
            [grain_ID,seed_xy,grain_phase,phase_now_pix] = Add_Grain( ...
                grain_ID,seed_xy,grain_phase,phase_now_pix,inc_mask,center,ip);
        end
    end
end

forbidden_phase = [];

if single_seed_only == 1 && ~isempty(fixed_phase)
    forbidden_phase = [forbidden_phase fixed_phase];
end

forbidden_phase = unique([forbidden_phase small_phase]);

free_mask      = grain_ID == 0;
free_area      = nnz(free_mask)/Npix * GRID.Lx * GRID.Ly;
grain_area_tar = pi/4*grain_size^2;
allowed_phase  = setdiff(find(phase_target_pix > phase_now_pix),forbidden_phase);

if isempty(allowed_phase)
    allowed_phase = setdiff(1:Nphase,forbidden_phase);
end

if isempty(Ngrain_user)
    Nmatrix = max(numel(allowed_phase),round(free_area/grain_area_tar));
else
    Nmatrix = Ngrain_user;
end

Nmatrix = max(Nmatrix,numel(allowed_phase));

[grain_local,seed_local] = Generate_Lloyd_Grains_Masked( ...
    GRID,free_mask,Nmatrix,periodic_map,lloyd_iter);

offset = max(grain_ID(:));
for ig = 1:Nmatrix
    mask = grain_local == ig;
    if any(mask(:))
        grain_ID(mask) = offset + ig;
        seed_xy(offset+ig,:) = seed_local(ig,:);
        grain_phase(offset+ig,1) = 0;
    end
end

[grain_ID,seed_xy,grain_phase] = Remove_Empty_Grains_With_Phase(grain_ID,seed_xy,grain_phase);

grain_phase = Assign_Grain_Phases_Rock( ...
    grain_ID,grain_phase,phase_target_pix,forbidden_phase, ...
    area_weight,overshoot_weight,same_phase_weight,assign_iter);

end

function phase_target_pix = Phase_Target_Pixels(phase_prop,Npix)

Nphase = numel(phase_prop);
raw    = phase_prop(:).'*Npix;
phase_target_pix = floor(raw);

present = phase_prop > 0;
phase_target_pix(present & phase_target_pix == 0) = 1;

rem = Npix - sum(phase_target_pix);
frac = raw - floor(raw);

if rem > 0
    [~,ord] = sort(frac,'descend');
    for k = 1:rem
        phase_target_pix(ord(mod(k-1,Nphase)+1)) = ...
            phase_target_pix(ord(mod(k-1,Nphase)+1)) + 1;
    end
elseif rem < 0
    [~,ord] = sort(frac,'ascend');
    k = 1;
    while rem < 0
        ip = ord(mod(k-1,Nphase)+1);
        if phase_target_pix(ip) > present(ip)
            phase_target_pix(ip) = phase_target_pix(ip) - 1;
            rem = rem + 1;
        end
        k = k + 1;
    end
end

end

function [grain_ID,seed_xy,grain_phase,phase_now_pix] = Add_Grain( ...
    grain_ID,seed_xy,grain_phase,phase_now_pix,mask,center,iph)

ig = max(grain_ID(:)) + 1;
grain_ID(mask) = ig;
seed_xy(ig,:) = center;
grain_phase(ig,1) = iph;
phase_now_pix(iph) = phase_now_pix(iph) + nnz(mask);

end

function mask = Closest_Pixels_Mask(X,Y,free_mask,center,npix)

idx  = find(free_mask);
npix = min(npix,numel(idx));

if npix <= 0
    mask = false(size(free_mask));
    return
end

D2 = (X(idx)-center(1)).^2 + (Y(idx)-center(2)).^2;
[~,ord] = sort(D2,'ascend');

mask = false(size(free_mask));
mask(idx(ord(1:npix))) = true;

end

function center = Random_Free_Center(GRID,free_mask)

idx = find(free_mask);
idx = idx(randi(numel(idx)));
[iy,ix] = ind2sub(size(free_mask),idx); %#ok<ASGLU>
center = [GRID.x(ix),GRID.y(iy)];

end

function [grain_ID,seed_xy] = Generate_Lloyd_Grains_Masked(GRID,free_mask,Ngrain,periodic_map,lloyd_iter)

[X,Y] = meshgrid(GRID.x,GRID.y);
idx   = find(free_mask);

if numel(idx) < Ngrain
    Ngrain = numel(idx);
end

ord = randperm(numel(idx),Ngrain);
[iy,ix] = ind2sub(size(free_mask),idx(ord)); %#ok<ASGLU>

sx = GRID.x(ix);
sy = GRID.y(iy);
seed_xy = [sx(:),sy(:)];

grain_ID = zeros(GRID.ny,GRID.nx);

for iter = 1:max(1,lloyd_iter)

    grain_ID = Assign_To_Seeds_Masked(GRID,free_mask,seed_xy,periodic_map);

    for ig = 1:Ngrain
        mask = grain_ID == ig;
        if any(mask(:))
            seed_xy(ig,1) = mean(X(mask));
            seed_xy(ig,2) = mean(Y(mask));
        else
            ii = idx(randi(numel(idx)));
            [iy0,ix0] = ind2sub(size(free_mask),ii);
            seed_xy(ig,1) = GRID.x(ix0);
            seed_xy(ig,2) = GRID.y(iy0);
        end
    end
end

grain_ID = Assign_To_Seeds_Masked(GRID,free_mask,seed_xy,periodic_map);

end

function grain_ID = Assign_To_Seeds_Masked(GRID,free_mask,seed_xy,periodic_map)

[X,Y] = meshgrid(GRID.x,GRID.y);
Ngrain = size(seed_xy,1);

dist2    = inf(GRID.ny,GRID.nx);
grain_ID = zeros(GRID.ny,GRID.nx);

xmin = min(GRID.x);
ymin = min(GRID.y); %#ok<NASGU>
Lx   = max(GRID.x)-xmin;
Ly   = max(GRID.y)-min(GRID.y);

for ig = 1:Ngrain

    ddx = abs(X-seed_xy(ig,1));
    ddy = abs(Y-seed_xy(ig,2));

    if periodic_map == 1
        ddx = min(ddx,Lx-ddx);
        ddy = min(ddy,Ly-ddy);
    end

    trial = ddx.^2 + ddy.^2;
    take  = trial < dist2 & free_mask;

    dist2(take)    = trial(take);
    grain_ID(take) = ig;
end

end

function grain_phase = Assign_Grain_Phases_Rock( ...
    grain_ID,grain_phase,phase_target_pix,forbidden_phase, ...
    area_weight,overshoot_weight,same_phase_weight,assign_iter)

Ngrain = max(grain_ID(:));
Nphase = numel(phase_target_pix);
Npix   = numel(grain_ID);

area = accumarray(grain_ID(:),1,[Ngrain,1]);
area_now = zeros(1,Nphase);
count    = zeros(1,Nphase);

for ig = 1:Ngrain
    iph = grain_phase(ig);
    if iph > 0
        area_now(iph) = area_now(iph) + area(ig);
        count(iph)    = count(iph) + 1;
    end
end

allowed_phase = setdiff(find(phase_target_pix > area_now),forbidden_phase);
if isempty(allowed_phase)
    allowed_phase = setdiff(1:Nphase,forbidden_phase);
end

adjw = Build_Grain_Adjacency_Weighted(grain_ID);
free = find(grain_phase(:).' == 0);
[~,ord0] = sort(area(free),'descend');
free = free(ord0);

for kk = 1:numel(free)

    ig = free(kk);
    candidates = allowed_phase;

    missing = candidates(count(candidates) == 0 & phase_target_pix(candidates) > 0);
    if numel(free)-kk+1 <= numel(missing)
        candidates = missing;
    end

    best_phase = candidates(1);
    best_score = inf;

    for ii = 1:numel(candidates)
        iph = candidates(ii);

        trial = area_now;
        trial(iph) = trial(iph) + area(ig);

        area_err = sum(((trial-phase_target_pix)/Npix).^2);
        overshot = sum((max(trial-phase_target_pix,0)/Npix).^2);
        contact  = Contact_Length_To_Phase(adjw,grain_phase,ig,iph)/max(nnz(adjw),1);

        score = area_weight*area_err + overshoot_weight*overshot + same_phase_weight*contact;

        if score < best_score
            best_score = score;
            best_phase = iph;
        end
    end

    grain_phase(ig) = best_phase;
    area_now(best_phase) = area_now(best_phase) + area(ig);
    count(best_phase) = count(best_phase) + 1;
end

best_obj = Phase_Assign_Objective(adjw,grain_phase,area,phase_target_pix, ...
    area_weight,same_phase_weight);

for iter = 1:assign_iter

    changed = 0;
    free = free(randperm(numel(free)));

    for kk = 1:numel(free)

        ig   = free(kk);
        iph0 = grain_phase(ig);
        best_phase = iph0;
        obj0 = best_obj;

        for ii = 1:numel(allowed_phase)
            iph = allowed_phase(ii);

            if iph == iph0
                continue
            end

            if count(iph0) <= 1 && phase_target_pix(iph0) > 0
                continue
            end

            trial_phase = grain_phase;
            trial_phase(ig) = iph;

            obj = Phase_Assign_Objective(adjw,trial_phase,area,phase_target_pix, ...
                area_weight,same_phase_weight);

            if obj < obj0
                obj0 = obj;
                best_phase = iph;
            end
        end

        if best_phase ~= iph0
            grain_phase(ig) = best_phase;
            count(iph0) = count(iph0) - 1;
            count(best_phase) = count(best_phase) + 1;
            best_obj = obj0;
            changed = 1;
        end
    end

    if changed == 0
        break
    end
end

grain_phase = grain_phase(:).';

end

function obj = Phase_Assign_Objective(adjw,grain_phase,area,phase_target_pix,area_weight,same_phase_weight)

Npix = sum(area);
Nphase = numel(phase_target_pix);
area_now = zeros(1,Nphase);

for ig = 1:numel(grain_phase)
    iph = grain_phase(ig);
    if iph > 0
        area_now(iph) = area_now(iph) + area(ig);
    end
end

area_err = sum(((area_now-phase_target_pix)/Npix).^2);
same_len = Same_Phase_Contact_Length(adjw,grain_phase);
tot_len  = sum(adjw(:))/2;

obj = area_weight*area_err + same_phase_weight*same_len/max(tot_len,1);

end

function len = Contact_Length_To_Phase(adjw,grain_phase,ig,iph)

neigh = find(adjw(ig,:) > 0);
len = 0;

for ii = 1:numel(neigh)
    jg = neigh(ii);
    if grain_phase(jg) == iph
        len = len + adjw(ig,jg);
    end
end

end

function len = Same_Phase_Contact_Length(adjw,grain_phase)

[i,j,w] = find(triu(adjw,1));
len = 0;

for k = 1:numel(i)
    if grain_phase(i(k)) == grain_phase(j(k))
        len = len + w(k);
    end
end

end

function adjw = Build_Grain_Adjacency_Weighted(grain_ID)

Ngrain = max(grain_ID(:));
adjw   = sparse(Ngrain,Ngrain);

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
    adjw(i,j) = adjw(i,j) + 1;
    adjw(j,i) = adjw(j,i) + 1;
end

end

function [grain_ID,seed_xy,grain_phase] = Remove_Empty_Grains_With_Phase(grain_ID,seed_xy,grain_phase)

present = unique(grain_ID(:));
present = present(present > 0);

Nold = max(grain_ID(:));
map  = zeros(Nold,1);
map(present) = 1:numel(present);

grain_ID = reshape(map(grain_ID(:)),size(grain_ID));
seed_xy = seed_xy(present,:);
grain_phase = grain_phase(present);

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

function chi = Fill_Initial_Chi(chi,pars_phase,c_ref,grain_ID,grain_phase,Ne)

for ig = 1:numel(grain_phase)

    iph  = grain_phase(ig);
    mask = grain_ID == ig;

    try
        R = PhaseThermo(pars_phase{iph},c_ref{iph});

        if isfield(R,'chi')
            for ie = 1:Ne
                for je = 1:Ne
                    val = R.chi{ie,je};
                    if iscell(val)
                        val = val{1};
                    end
                    if numel(val) > 1
                        val = val(1);
                    end
                    chi{ie,je}(mask) = val;
                end
            end
        end
    catch
        % If PhaseThermo does not return chi for this phase, leave zero.
        % The solver LE update will fill it later.
    end
end

end

function Plot_Rock_Map(x,y,L_sc,phase_ID,grain_ID)

figure(3); clf
imagesc(x*L_sc*1e6,y*L_sc*1e6,phase_ID)
set(gca,'YDir','normal')
axis image tight
colorbar
hold on

B = false(size(grain_ID));
B(:,1:end-1) = B(:,1:end-1) | grain_ID(:,1:end-1) ~= grain_ID(:,2:end);
B(1:end-1,:) = B(1:end-1,:) | grain_ID(1:end-1,:) ~= grain_ID(2:end,:);
[iy,ix] = find(B);
plot(x(ix)*L_sc*1e6,y(iy)*L_sc*1e6,'k.','MarkerSize',1)

xlabel('x (\mum)')
ylabel('y (\mum)')
title('Initial rock-like phase map')
drawnow

end
