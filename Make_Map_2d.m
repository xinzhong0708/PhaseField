%% Clear and restart
clear; figure(3); clf
addpath([cd,'\bin'])
addpath([cd,'\Thermo'])
addpath([cd,'\Thermo\Solutions'])

%% ------------------------------------------------------------------------
%  Scaling / basic physical constants used by the map
% -------------------------------------------------------------------------
PHYS             = struct();
PHYS.E_sc        = 1e9;
PHYS.L_sc        = 1;

E_sc             = PHYS.E_sc;     % legacy save
L_sc             = PHYS.L_sc;     % legacy save

eta0             = 1200e10/E_sc;

%% ------------------------------------------------------------------------
%  Thermodynamic data
% -------------------------------------------------------------------------
phs_name         = {'Olv','Cpx','Grt','Qtz','Crd','Kya','Opx'};
pars             = Load_Data(phs_name);
% pars             = Load_Data({'Qtz','Crd'});

Np               = length(pars);

%% ------------------------------------------------------------------------
%  Model functions
% -------------------------------------------------------------------------
MODEL            = struct();
MODEL.pars       = pars;
MODEL.phs_name   = phs_name;

eps_phi          = 1e-14;

MODEL.p_fun      = @(a,phi) phi(:,:,a).^2 ./ (sum(phi.^2,3) + eps_phi);

MODEL.dpdphi     = @(a,b,phi) (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) - 2*phi(:,:,a).*phi(:,:,b).^2 ./ (sum(phi.^2,3) + eps_phi).^2;

% legacy name, for older functions
F                = MODEL;

%% ------------------------------------------------------------------------
%  Grid
% -------------------------------------------------------------------------
Lx               = 5e-6;
Ly               = 5e-6;

nx               = 60;
ny               = 60;

x                = linspace(0,Lx,nx);
y                = linspace(0,Ly,ny);

dx               = x(2)-x(1);
dy               = y(2)-y(1);

% scale grid
x                = x/L_sc;
y                = y/L_sc;
dx               = dx/L_sc;
dy               = dy/L_sc;

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
%  Initial endmember compositions
% -------------------------------------------------------------------------
c                = cell(1,Np);


% Olv
c{1}{1}          = 0.0063*ones(ny,nx);
c{1}{2}          = 0.2597*ones(ny,nx);
c{1}{3}          = 0.5133*ones(ny,nx);

% Cpx
c{2}{1}          = 0.0200*ones(ny,nx);
c{2}{2}          = 0.1200*ones(ny,nx);
c{2}{3}          = 0.0300*ones(ny,nx);
c{2}{4}          = 0.6000*ones(ny,nx);

% Grt
c{3}{1}          = 0.4818*ones(ny,nx);
c{3}{2}          = 0.3943*ones(ny,nx);

% Qtz
c{4}{1}          = 1.0000*ones(ny,nx);

% Crd
c{5}{1}          = 1.0000*ones(ny,nx);

% Kya
c{6}{1}          = 1.0000*ones(ny,nx);

% Opx
c{7}{1}          = 0.4000*ones(ny,nx);
c{7}{2}          = 0.3000*ones(ny,nx);
c{7}{3}          = 0.1000*ones(ny,nx);
c{7}{4}          = 0.0500*ones(ny,nx);

% Model
MODEL.phase_index= [1,2,3,4,5,6,7];
c_phase          = c;

% % Qtz
% c{1}{1}          = 1.0000*ones(ny,nx);
% 
% % Crd
% c{2}{1}          = 1.0000*ones(ny,nx);



%% ------------------------------------------------------------------------
%  Initial phase field
% -------------------------------------------------------------------------
phi              = zeros(ny,nx,Np);

% md               = nx/7;
% phi(:,1:md,4)              = 1;
% phi(:,1*md+1:2*md,2)       = 1;
% phi(:,2*md+1:3*md,3)       = 1;
% phi(:,3*md+1:4*md,1)       = 1;
% phi(:,4*md+1:5*md,5)       = 1;
% phi(:,5*md+1:6*md,7)       = 1;
% phi(:,6*md+1:7*md,6)       = 1;




%%  Initial phase field: random polygonal grains
% -------------------------------------------------------------------------
%% ------------------------------------------------------------------------
%  Initial phase field: random polygonal grains
% -------------------------------------------------------------------------
Nphase                   = length(pars);
MAP_OPT                  = struct();
MAP_OPT.rng_seed         = 7;
MAP_OPT.grain_size       = 1e-6/L_sc;
MAP_OPT.n_grain          = [];              % [] = estimate from grain_size

% Phase assignment to grains
% Options: 'cycle', 'random', 'fractions'
MAP_OPT.phase_mode       = 'random';
MAP_OPT.phase_fraction   =  ones(1,Nphase)/Nphase;

[phi,phase_ID,grain_ID,seed_xy,grain_phase] = ...
    Make_RandomPolygon_GrainMap(GRID,Nphase,MAP_OPT);

Ngrain                   = size(phi,3);

MODEL.phase_index        = grain_phase(:).';        % 1 x Ngrain
MODEL.pars               = pars(MODEL.phase_index);

Np                       = Ngrain;                  % legacy: number of order parameters/grains
pars                     = MODEL.pars;              % repeated pars, grain-sized

% Expand compositions from thermodynamic phase to grain
c = cell(1,Ngrain);

for ig = 1:Ngrain
    iph   = MODEL.phase_index(ig);
    c{ig} = c_phase{iph};
end

figure(3); clf
pcolor(GRID.x*1e6,GRID.y*1e6,phase_ID);
shading flat
axis equal tight
colorbar
title('Initial random polygonal thermodynamic phase map')
xlabel('x \mum')
ylabel('y \mum')
drawnow




%% ------------------------------------------------------------------------
%  Initial p, e, E
% -------------------------------------------------------------------------
p                = Calc_p(MODEL,phi);
e                = Calc_e(pars,c);
E                = Calc_E_Tot(e,p);

Ne               = length(E);

%% ------------------------------------------------------------------------
%  Initial thermodynamic state containers
% -------------------------------------------------------------------------
mu_e             = repmat({zeros(ny,nx)},1,Ne);
chi              = repmat({zeros(ny,nx)},Ne,Ne);
eta              = eta0*ones(ny,nx);

%% ------------------------------------------------------------------------
%  PARAM structure for local-equilibrium initialization
% -------------------------------------------------------------------------
PARAM            = struct();
PARAM.Np         = Np;
PARAM.Ne         = Ne;
PARAM.eta        = eta;

%% ------------------------------------------------------------------------
%  Pair-pair initialization in pure phase regions
% -------------------------------------------------------------------------
tol_pure         = 1e-12;
p_pair           = 0.05;

for ip = 1:Np

    % Pure region of phase ip
    mask = p(:,:,ip) > 1 - tol_pure;

    if ~any(mask(:))
        continue
    end

    id    = find(mask);
    Nmask = numel(id);

    pid_other     = 1:Np;
    pid_other(ip) = [];

    for io = pid_other

        fprintf('Checking pair initialization: %d %d\n',ip,io)

        % -------------------------------------------------------------
        % Artificial pair p only in pure ip region
        % -------------------------------------------------------------
        pp          = p;

        pp_ip       = pp(:,:,ip);
        pp_io       = pp(:,:,io);

        pp_ip(mask) = 1 - p_pair;
        pp_io(mask) = p_pair;

        pp(:,:,ip)  = pp_ip;
        pp(:,:,io)  = pp_io;

        % -------------------------------------------------------------
        % Local bulk composition from current c and artificial pp
        % -------------------------------------------------------------
        e_test      = Calc_e(pars,c);
        E_test      = Calc_E_Tot(e_test,pp);

        % -------------------------------------------------------------
        % Slice p into 1 x Nmask x Np
        % -------------------------------------------------------------
        p_slice     = zeros(1,Nmask,Np);

        for jp = 1:Np
            tmp             = pp(:,:,jp);
            p_slice(1,:,jp) = tmp(mask).';
        end

        % -------------------------------------------------------------
        % Slice c, E, eta
        % -------------------------------------------------------------
        c_slice     = Slice_c_Local(c,mask);
        E_slice     = Slice_E_Local(E_test,mask);

        eta_slice   = eta(mask).';
        eta_slice   = reshape(eta_slice,1,Nmask);

        mu_slice    = repmat({zeros(1,Nmask)},1,Ne);
        chi_slice   = repmat({zeros(1,Nmask)},Ne,Ne);

        % -------------------------------------------------------------
        % Build local structured STATE and PARAM
        % -------------------------------------------------------------
        STATE_LOC          = struct();
        STATE_LOC.c        = c_slice;
        STATE_LOC.e        = Calc_e(pars,c_slice);
        STATE_LOC.E        = E_slice;
        STATE_LOC.mu_e     = mu_slice;
        STATE_LOC.chi      = chi_slice;
        STATE_LOC.omg      = zeros(1,Nmask,Np);
        STATE_LOC.phi      = [];
        STATE_LOC.p        = p_slice;
        STATE_LOC.mask     = ones(1,Nmask,Np);
        STATE_LOC.LE_state = [];

        PARAM_LOC          = PARAM;
        PARAM_LOC.eta      = eta_slice;
        PARAM_LOC.Np       = Np;
        PARAM_LOC.Ne       = Ne;

        % -------------------------------------------------------------
        % Run local equilibrium using structured form
        % -------------------------------------------------------------
        STATE_LOC = LE_Run(STATE_LOC,PARAM_LOC,MODEL);

        % -------------------------------------------------------------
        % Put back only the two pair phases
        % -------------------------------------------------------------
        for ic = 1:length(STATE_LOC.c{io})
            tmp       = c{io}{ic};
            tmp(mask) = STATE_LOC.c{io}{ic}(:);
            c{io}{ic} = tmp;
        end

        for ic = 1:length(STATE_LOC.c{ip})
            tmp       = c{ip}{ic};
            tmp(mask) = STATE_LOC.c{ip}{ic}(:);
            c{ip}{ic} = tmp;
        end

    end
end

%% ------------------------------------------------------------------------
%  Recalculate final map fields after pair initialization
% -------------------------------------------------------------------------
p                = Calc_p(MODEL,phi);
e                = Calc_e(pars,c);
E                = Calc_E_Tot(e,p);

mu_e             = repmat({zeros(ny,nx)},1,Ne);
chi              = repmat({zeros(ny,nx)},Ne,Ne);

%% ------------------------------------------------------------------------
%  Final structured state
% -------------------------------------------------------------------------
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

PARAM.eta        = eta;
PARAM.Np         = Np;
PARAM.Ne         = Ne;

%% ------------------------------------------------------------------------
%  Save structured + legacy variables
% -------------------------------------------------------------------------
save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'pars','E_sc','Lx','Ly','c','E','e','p','phi','eta','mu_e','chi', ...
    'x','dx','nx','y','dy','ny','L_sc','F','Np','Ne')


% plot(GRID.x,STATE.c{1}{1}(3,:),GRID.x,STATE.c{1}{2}(3,:),GRID.x,STATE.c{1}{3}(3,:),GRID.x,STATE.c{1}{4}(3,:));title('c21')


%% ========================================================================
%  Local helper functions
% ========================================================================

function c_slice = Slice_c_Local(c,mask)

Np      = length(c);
Nmask   = nnz(mask);
c_slice = cell(size(c));

for ip = 1:Np
    Nc = length(c{ip});
    c_slice{ip} = cell(size(c{ip}));

    for ic = 1:Nc
        tmp = c{ip}{ic};
        c_slice{ip}{ic} = reshape(tmp(mask),1,Nmask);
    end
end

end

function E_slice = Slice_E_Local(E,mask)

Ne      = length(E);
Nmask   = nnz(mask);
E_slice = cell(size(E));

for ie = 1:Ne
    tmp = E{ie};
    E_slice{ie} = reshape(tmp(mask),1,Nmask);
end

end


function [phi,phase_ID,grain_ID,seed_xy,grain_phase] = Make_RandomPolygon_GrainMap(GRID,Nphase,opt)
%MAKE_RANDOMPOLYGON_GRAINMAP
%
% Generate polygonal grains.
%
% phi(:,:,ig) is one order parameter per grain.
% grain_phase(ig) gives the thermodynamic phase index of grain ig.
%
% phase_ID is the thermodynamic phase map for plotting.
% grain_ID is the grain index map.

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
    opt.phase_mode = 'random';
end

if ~isfield(opt,'phase_fraction') || isempty(opt.phase_fraction)
    opt.phase_fraction = ones(1,Nphase)/Nphase;
end

rng(opt.rng_seed)

x = GRID.x(:).';
y = GRID.y(:);

[X,Y] = meshgrid(x,y);

[ny,nx] = size(X);

xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

Lx = xmax - xmin;
Ly = ymax - ymin;

if isempty(opt.n_grain)

    area_box = max(Lx*Ly,eps);
    area_g   = max(opt.grain_size^2,eps);

    Ngrain = max(Nphase,round(area_box/area_g));

else

    Ngrain = max(Nphase,round(opt.n_grain));

end

% Random grain seeds
sx = xmin + Lx*rand(Ngrain,1);
sy = ymin + Ly*rand(Ngrain,1);

seed_xy = [sx,sy];

% Assign each grid cell to nearest seed
grain_ID = zeros(ny,nx);
bestD2   = inf(ny,nx);

for ig = 1:Ngrain

    D2 = (X - sx(ig)).^2 + (Y - sy(ig)).^2;

    mask = D2 < bestD2;

    bestD2(mask)   = D2(mask);
    grain_ID(mask) = ig;

end

% Assign thermodynamic phase to each grain
grain_phase = zeros(Ngrain,1);

switch lower(opt.phase_mode)

    case 'cycle'

        for ig = 1:Ngrain
            grain_phase(ig) = mod(ig-1,Nphase) + 1;
        end

    case 'fractions'

        frac = opt.phase_fraction(:).';
        frac = frac / sum(frac);

        edge = [0,cumsum(frac)];
        edge(end) = 1;

        r = rand(Ngrain,1);

        for ig = 1:Ngrain
            grain_phase(ig) = find(r(ig) >= edge(1:end-1) & r(ig) <= edge(2:end),1,'first');
        end

    case 'random'

        grain_phase = randi(Nphase,Ngrain,1);

    otherwise

        error('Unknown phase_mode: %s',opt.phase_mode)

end

% Make sure every thermodynamic phase appears at least once
if Ngrain >= Nphase
    grain_phase(1:Nphase) = 1:Nphase;
    grain_phase = grain_phase(randperm(Ngrain));
end

% Thermodynamic phase map
phase_ID = zeros(ny,nx);

for ig = 1:Ngrain
    phase_ID(grain_ID == ig) = grain_phase(ig);
end

% One phi field per grain
phi = zeros(ny,nx,Ngrain);

for ig = 1:Ngrain
    phi(:,:,ig) = double(grain_ID == ig);
end

end