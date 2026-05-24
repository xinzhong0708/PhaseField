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
eta0             = 1000e10/E_sc;      %Use the same eta in Run_2D

%% ------------------------------------------------------------------------
%  Thermodynamic phases and Gibbs-equilibrium starting values
%  Edit this block when changing the phase assemblage.
% -------------------------------------------------------------------------
phs_name         = {'Grt','Cpx','Fel'};
pars_phase       = Load_Data(phs_name);
Nphase           = length(pars_phase);

%Target phase proportions from Gibbs minimization, same order as phs_name
phase_prop       = [0.20 0.40 0.40];
phase_prop       = phase_prop/sum(phase_prop);

%Independent endmember proportions from Gibbs minimization
c_value          = cell(1,Nphase);
c_value{1}       = [0.42 0.48];                    %Grt
c_value{2}       = [0.11 0.26 -0.00 0.00 0.52];  %Cpx
c_value{3}       = [0.48];                         %Fel

c_guess          = cell(1,Nphase);

if numel(phase_prop) ~= Nphase || numel(c_value) ~= Nphase
    error('phs_name, phase_prop and c_value must contain the same number of phases.')
end

for ip = 1:Nphase

    Nc = size(pars_phase{ip}.n,1)-1;

    if numel(c_value{ip}) ~= Nc
        error('%s requires %d independent endmember proportions; %d were supplied.', ...
              phs_name{ip},Nc,numel(c_value{ip}))
    end

    c_guess{ip} = num2cell(c_value{ip});

    if any([c_value{ip},1-sum(c_value{ip})] < 0)
        warning('%s initial endmember proportions include a negative value.',phs_name{ip})
    end

end

%% ------------------------------------------------------------------------
%  Grid
% -------------------------------------------------------------------------
Lx               = 5e-6;
Ly               = 5e-6;

nx               = 70;
ny               = 70;

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
%  Random polygonal grain map
%
%  MAP_OPT.grain_size controls the approximate polygon size.
%  MAP_OPT.phase_fraction controls the target areal proportion of phases.
%  'balanced' assigns phases to the generated polygons to approximately
%  match phase_fraction while preserving complete grains.
% -------------------------------------------------------------------------
MAP_OPT                  = struct();
MAP_OPT.rng_seed         = 15;
MAP_OPT.grain_size       = 3e-6/L_sc;
MAP_OPT.n_grain          = [];
MAP_OPT.phase_mode       = 'balanced';        %'balanced', 'random', 'cycle'
MAP_OPT.phase_fraction   = phase_prop;
MAP_OPT.periodic         = false;             %Set true only for periodic geometry

[phi,phase_ID,grain_ID,seed_xy,grain_phase,phase_prop_map] = ...
    Make_RandomPolygon_GrainMap(GRID,Nphase,MAP_OPT);

Ngrain           = size(phi,3);
Np               = Ngrain;

fprintf('\nTarget and realized map phase proportions:\n')
fprintf('Phase       target          map\n')

for ip = 1:Nphase
    fprintf('%-8s    %.8f      %.8f\n',phs_name{ip},phase_prop(ip),phase_prop_map(ip))
end

fprintf('Maximum phase-proportion mismatch = %.8e\n', ...
        max(abs(phase_prop-phase_prop_map)))

%% ------------------------------------------------------------------------
%  One reference LE solve using the realized map phase proportions
%
%  A polygonal map cannot generally match the requested phase proportions
%  exactly. Using phase_prop_map here ensures that the generated map and
%  the reference conserved composition represent the same initial system.
% -------------------------------------------------------------------------
p_ref            = reshape(phase_prop_map,1,1,Nphase);
e_guess          = Calc_e(pars_phase,c_guess);
E_target         = Calc_E_Tot(e_guess,p_ref);
Ne               = length(E_target);

[c_ref,mu_ref]   = LE_Calculator( ...
                    pars_phase,p_ref,c_guess,E_target,eta0,[0.1,2000]);

e_ref            = Calc_e(pars_phase,c_ref);
E_mix_ref        = Calc_E_Tot(e_ref,p_ref);
E_offset         = cell(1,Ne);
omega_ref        = zeros(1,Nphase);

for ie = 1:Ne
    E_offset{ie} = E_target{ie} - E_mix_ref{ie};
end

for ip = 1:Nphase

    omega_ref(ip) = PhaseG(pars_phase{ip},c_ref{ip});

    for ie = 1:Ne
        omega_ref(ip) = omega_ref(ip) - e_ref{ip}{ie}.*mu_ref{ie};
    end

end

fprintf('\nReference LE state based on realized map proportions:\n')

for ip = 1:Nphase
    fprintf('%-8s  phase = %.8f,  independent endmembers = ', ...
            phs_name{ip},phase_prop_map(ip))
    disp(cell2mat(c_ref{ip}(:)).')
end

fprintf('max|omega - omega(1)| = %.8e\n',max(abs(omega_ref-omega_ref(1))))

fprintf('\nFinite-eta offset check: E_offset and mu/eta should agree:\n')

for ie = 1:Ne
    fprintf('Element %d: E_offset = % .8e, mu/eta = % .8e, difference = %.3e\n', ...
            ie,E_offset{ie},mu_ref{ie}/eta0,abs(E_offset{ie}-mu_ref{ie}/eta0))
end

%% ------------------------------------------------------------------------
%  Spatially uniform LE-consistent compositions for every phase
% -------------------------------------------------------------------------
c_phase          = cell(1,Nphase);

for ip = 1:Nphase

    c_phase{ip} = cell(size(c_ref{ip}));

    for ic = 1:length(c_ref{ip})
        c_phase{ip}{ic} = c_ref{ip}{ic}*ones(ny,nx);
    end

end

%% ------------------------------------------------------------------------
%  Grain-resolved MODEL
% -------------------------------------------------------------------------
MODEL                    = struct();
MODEL.phs_name           = phs_name;
MODEL.phase_index        = grain_phase(:).';
MODEL.pars               = pars_phase(MODEL.phase_index);

eps_phi                  = 1e-14;
MODEL.p_fun              = @(a,phi) phi(:,:,a).^2 ./ ...
                            (sum(phi.^2,3) + eps_phi);
MODEL.dpdphi             = @(a,b,phi) ...
                            (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) ...
                            - 2*phi(:,:,a).*phi(:,:,b).^2 ./ ...
                            (sum(phi.^2,3) + eps_phi).^2;

pars                     = MODEL.pars;
F                        = MODEL;

figure(3); clf
pcolor(x*1e6,y*1e6,phase_ID)
shading flat
axis equal tight
colorbar
title('Initial random polygonal thermodynamic phase map')
xlabel('x \mum')
ylabel('y \mum')
drawnow

%% ------------------------------------------------------------------------
%  Construct the penalty-consistent conserved composition field
% -------------------------------------------------------------------------
c                = Expand_c_By_Phase(c_phase,MODEL.phase_index);
p                = Calc_p(MODEL,phi);
e                = Calc_e(pars,c);
E                = Calc_E_Tot(e,p);
eta              = eta0*ones(ny,nx);

%At finite eta, stationary LE requires E - Emix = mu/eta.
%The common offset keeps all grains consistent with the reference state.
for ie = 1:Ne
    E{ie} = E{ie} + E_offset{ie};
end

mu_e             = repmat({zeros(ny,nx)},1,Ne);
chi              = repmat({zeros(ny,nx)},Ne,Ne);

%% ------------------------------------------------------------------------
%  PARAM and STATE
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
STATE.omg        = zeros(ny,nx,Np);
STATE.phi        = phi;
STATE.p          = p;
STATE.mask       = ones(ny,nx,Np);
STATE.LE_state   = [];

%% ------------------------------------------------------------------------
%  Initial LE check before dynamic evolution
% -------------------------------------------------------------------------
STATE_INI        = STATE;
STATE            = LE_Run(STATE,PARAM,MODEL);

omega_phase      = zeros(ny,nx,Nphase);

for ip = 1:Nphase
    ig                  = find(MODEL.phase_index == ip,1,'first');
    omega_phase(:,:,ip) = STATE.omg(:,:,ig);
end

fprintf('\nInitial map check after LE_Run:\n')

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
    for ic = 1:length(STATE.c{ig})
        dc_max = max(dc_max,max(abs(STATE.c{ig}{ic}(:)-STATE_INI.c{ig}{ic}(:))));
    end
end

E_mean           = zeros(Ne,1);
E_target_vec     = zeros(Ne,1);

for ie = 1:Ne
    E_mean(ie)       = mean(STATE.E{ie}(:));
    E_target_vec(ie) = E_target{ie};
end

fprintf('max|c after first LE - c before LE| = %.8e\n',dc_max)
fprintf('max|mean(E map) - E target|         = %.8e\n',max(abs(E_mean-E_target_vec)))

%Keep legacy variables consistent with STATE after LE evaluation
c                = STATE.c;
e                = STATE.e;
E                = STATE.E;
p                = STATE.p;
mu_e             = STATE.mu_e;
chi              = STATE.chi;

%% ------------------------------------------------------------------------
%  Save structured + legacy variables
% -------------------------------------------------------------------------
save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'pars','pars_phase','c_value','c_guess','c_phase','c_ref','mu_ref', ...
    'E_target','E_offset','phase_prop','phase_prop_map','phs_name', ...
    'E_sc','Lx','Ly','c','E','e','p','phi','eta','mu_e','chi', ...
    'x','dx','nx','y','dy','ny','L_sc','F','Np','Ne','Nphase','Ngrain', ...
    'phase_ID','grain_ID','grain_phase','seed_xy','MAP_OPT')

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


function [phi,phase_ID,grain_ID,seed_xy,grain_phase,phase_prop_map] = ...
         Make_RandomPolygon_GrainMap(GRID,Nphase,opt)
%MAKE_RANDOMPOLYGON_GRAINMAP Generate a random polygonal grain map.
%
% phi(:,:,ig)       order parameter for grain ig.
% grain_phase(ig)   thermodynamic phase index assigned to grain ig.
% phase_prop_map    realized areal fraction of each thermodynamic phase.

if ~isfield(opt,'rng_seed') || isempty(opt.rng_seed)
    opt.rng_seed = 1;
end

if ~isfield(opt,'grain_size') || isempty(opt.grain_size)
    opt.grain_size = 0.5*min(GRID.Lx,GRID.Ly);
end

if ~isfield(opt,'n_grain')
    opt.n_grain = [];
end

if ~isfield(opt,'phase_mode') || isempty(opt.phase_mode)
    opt.phase_mode = 'balanced';
end

if ~isfield(opt,'phase_fraction') || isempty(opt.phase_fraction)
    opt.phase_fraction = ones(1,Nphase)/Nphase;
end

if ~isfield(opt,'periodic') || isempty(opt.periodic)
    opt.periodic = false;
end

opt.phase_fraction = opt.phase_fraction(:).';
opt.phase_fraction = opt.phase_fraction/sum(opt.phase_fraction);

rng(opt.rng_seed)

x       = GRID.x(:).';
y       = GRID.y(:);
[X,Y]   = meshgrid(x,y);
[ny,nx] = size(X);

xmin    = min(x);
xmax    = max(x);
ymin    = min(y);
ymax    = max(y);

Lx      = xmax-xmin;
Ly      = ymax-ymin;

if isempty(opt.n_grain)
    Nseed = max(Nphase,round(Lx*Ly/max(opt.grain_size^2,eps)));
else
    Nseed = max(Nphase,round(opt.n_grain));
end

sx      = xmin + Lx*rand(Nseed,1);
sy      = ymin + Ly*rand(Nseed,1);
seed_xy = [sx,sy];

grain_ID = zeros(ny,nx);
bestD2   = inf(ny,nx);

for ig = 1:Nseed

    if opt.periodic
        dX = abs(X-sx(ig));
        dY = abs(Y-sy(ig));
        dX = min(dX,Lx-dX);
        dY = min(dY,Ly-dY);
    else
        dX = X-sx(ig);
        dY = Y-sy(ig);
    end

    D2   = dX.^2 + dY.^2;
    mask = D2 < bestD2;

    bestD2(mask)   = D2(mask);
    grain_ID(mask) = ig;

end

%Remove seeds that receive no grid point and renumber grains consecutively
used     = unique(grain_ID(:)).';
seed_xy  = seed_xy(used,:);
grain_ID_old = grain_ID;
grain_ID = zeros(ny,nx);

for ig = 1:numel(used)
    grain_ID(grain_ID_old == used(ig)) = ig;
end

Ngrain = numel(used);

switch lower(opt.phase_mode)

    case 'balanced'
        grain_phase = Assign_Phases_By_Area( ...
                      grain_ID,Nphase,opt.phase_fraction);

    case 'random'
        edge        = [0 cumsum(opt.phase_fraction)];
        r           = rand(Ngrain,1);
        grain_phase = zeros(Ngrain,1);

        for ig = 1:Ngrain
            grain_phase(ig) = find(r(ig) <= edge(2:end),1,'first');
        end

    case 'cycle'
        grain_phase = mod((1:Ngrain).'-1,Nphase)+1;

    otherwise
        error('Unknown MAP_OPT.phase_mode: %s',opt.phase_mode)
end

if any(accumarray(grain_phase,1,[Nphase,1]) == 0)
    error('At least one thermodynamic phase is absent from the generated grain map.')
end

phase_ID = zeros(ny,nx);
phi      = zeros(ny,nx,Ngrain);

for ig = 1:Ngrain
    mask               = grain_ID == ig;
    phase_ID(mask)     = grain_phase(ig);
    phi(:,:,ig)        = double(mask);
end

phase_prop_map = zeros(1,Nphase);

for ip = 1:Nphase
    phase_prop_map(ip) = nnz(phase_ID == ip)/numel(phase_ID);
end

end


function grain_phase = Assign_Phases_By_Area(grain_ID,Nphase,phase_fraction)
%ASSIGN_PHASES_BY_AREA Assign whole grains close to target phase fractions.

Ngrain      = max(grain_ID(:));
grain_area  = accumarray(grain_ID(:),1,[Ngrain,1]).';
grain_area  = grain_area/sum(grain_area);

[~,order]   = sort(grain_area,'descend');

grain_phase = zeros(Ngrain,1);
area_now    = zeros(1,Nphase);
count_now   = zeros(1,Nphase);

for ii = 1:Ngrain

    ig      = order(ii);
    missing = find(count_now == 0);

    if Ngrain-ii+1 == numel(missing)
        candidate = missing;
    else
        candidate = 1:Nphase;
    end

    loss = inf(1,Nphase);

    for ip = candidate
        area_try     = area_now;
        area_try(ip) = area_try(ip) + grain_area(ig);
        loss(ip)     = sum((area_try-phase_fraction).^2);
    end

    [~,ip]          = min(loss);
    grain_phase(ig) = ip;
    area_now(ip)    = area_now(ip) + grain_area(ig);
    count_now(ip)   = count_now(ip) + 1;

end

%Improve the area match by moving whole grains when this reduces error
changed = true;

while changed

    changed = false;

    for ig = 1:Ngrain

        ip0 = grain_phase(ig);

        if count_now(ip0) <= 1
            continue
        end

        loss0 = sum((area_now-phase_fraction).^2);

        for ip = setdiff(1:Nphase,ip0)

            area_try      = area_now;
            area_try(ip0) = area_try(ip0)-grain_area(ig);
            area_try(ip)  = area_try(ip) +grain_area(ig);
            loss          = sum((area_try-phase_fraction).^2);

            if loss < loss0

                grain_phase(ig) = ip;
                area_now        = area_try;
                count_now(ip0)  = count_now(ip0)-1;
                count_now(ip)   = count_now(ip)+1;
                changed         = true;
                break

            end

        end

    end

end

end
