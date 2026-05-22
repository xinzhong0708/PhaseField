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

E_sc             = PHYS.E_sc;
L_sc             = PHYS.L_sc;

eta0             = 2000e10/E_sc;

%% ------------------------------------------------------------------------
%  Thermodynamic data: unique thermodynamic phases
% -------------------------------------------------------------------------
phs_name         = {'Olv','Cpx','Grt','Qtz','Crd','Kya','Opx'};
pars_phase       = Load_Data(phs_name);
Nphase           = length(pars_phase);

%% ------------------------------------------------------------------------
%  Grid
% -------------------------------------------------------------------------
Lx               = 5e-6;
Ly               = 5e-6;

nx               = 100;
ny               = 100;

x                = linspace(0,Lx,nx);
y                = linspace(0,Ly,ny);

dx               = x(2)-x(1);
dy               = y(2)-y(1);

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
%  Initial endmember compositions: unique thermodynamic phases
% -------------------------------------------------------------------------
c_phase          = cell(1,Nphase);

% Olv
c_phase{1}{1}    = 0.0063*ones(ny,nx);
c_phase{1}{2}    = 0.2597*ones(ny,nx);
c_phase{1}{3}    = 0.5133*ones(ny,nx);

% Cpx
c_phase{2}{1}    = 0.0200*ones(ny,nx);
c_phase{2}{2}    = 0.1200*ones(ny,nx);
c_phase{2}{3}    = 0.0300*ones(ny,nx);
c_phase{2}{4}    = 0.6000*ones(ny,nx);

% Grt
c_phase{3}{1}    = 0.4818*ones(ny,nx);
c_phase{3}{2}    = 0.3943*ones(ny,nx);

% Qtz
c_phase{4}{1}    = 1.0000*ones(ny,nx);

% Crd
c_phase{5}{1}    = 1.0000*ones(ny,nx);

% Kya
c_phase{6}{1}    = 1.0000*ones(ny,nx);

% Opx
c_phase{7}{1}    = 0.4000*ones(ny,nx);
c_phase{7}{2}    = 0.3000*ones(ny,nx);
c_phase{7}{3}    = 0.1000*ones(ny,nx);
c_phase{7}{4}    = 0.0500*ones(ny,nx);


%% ------------------------------------------------------------------------
%  Initial phase field: random polygonal grains
% -------------------------------------------------------------------------
MAP_OPT                     = struct();
MAP_OPT.rng_seed            = 5;
MAP_OPT.grain_size          = 7e-7/L_sc;
MAP_OPT.n_grain             = [];
MAP_OPT.phase_mode          = 'random';              % 'cycle', 'random', 'fractions'
MAP_OPT.separate_same_phase = true;
MAP_OPT.phase_fraction      =  ones(1,Nphase)/Nphase;
[phi,phase_ID,grain_ID,seed_xy,grain_phase] = Make_RandomPolygon_GrainMap(GRID,Nphase,MAP_OPT);
Ngrain                      = size(phi,3);
Np                          = Ngrain;                % legacy: number of order parameters/grains

% %% ------------------------------------------------------------------------
% %  Initial phase field: pseudo-2D 1D phase bands along x
% % -------------------------------------------------------------------------
% Ngrain   = Nphase;
% Np       = Ngrain;     % legacy: number of order parameters/grains
% 
% phi      = zeros(GRID.ny,GRID.nx,Ngrain);
% phase_ID = zeros(GRID.ny,GRID.nx);
% grain_ID = zeros(GRID.ny,GRID.nx);
% 
% x_edges  = round(linspace(1,GRID.nx+1,Nphase+1));
% 
% for iph = 1:Nphase
% 
%     ix1 = x_edges(iph);
%     ix2 = x_edges(iph+1) - 1;
% 
%     if iph == Nphase
%         ix2 = GRID.nx;
%     end
% 
%     phi(:,ix1:ix2,iph)      = 1;
%     phase_ID(:,ix1:ix2)     = iph;
%     grain_ID(:,ix1:ix2)     = iph;
% 
% end
% 
% grain_phase = (1:Nphase).';
% seed_xy     = zeros(Ngrain,2);
% for iph = 1:Nphase
%     ix1 = x_edges(iph);
%     ix2 = x_edges(iph+1) - 1;
%     if iph == Nphase
%         ix2 = GRID.nx;
%     end
%     seed_xy(iph,1) = mean(GRID.x(ix1:ix2));
%     seed_xy(iph,2) = mean(GRID.y);
% end

%% ------------------------------------------------------------------------
%  Grain-resolved MODEL
% -------------------------------------------------------------------------
MODEL                    = struct();
MODEL.phs_name           = phs_name;
MODEL.phase_index        = grain_phase(:).';      % grain -> thermodynamic phase
MODEL.pars               = pars_phase(MODEL.phase_index);
eps_phi                  = 1e-14;
MODEL.p_fun              = @(a,phi) phi(:,:,a).^2 ./ (sum(phi.^2,3) + eps_phi);
MODEL.dpdphi             = @(a,b,phi) (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) - 2*phi(:,:,a).*phi(:,:,b).^2 ./ (sum(phi.^2,3) + eps_phi).^2;
pars                     = MODEL.pars;            % grain-sized repeated pars
F                        = MODEL;                 % legacy name

figure(3); clf
pcolor(GRID.x*1e6,GRID.y*1e6,phase_ID)
shading flat
axis equal tight
colorbar
title('Initial random polygonal thermodynamic phase map')
xlabel('x \mum')
ylabel('y \mum')
drawnow

%% ------------------------------------------------------------------------
%  Pair-pair initialization on unique thermodynamic phases only
% -------------------------------------------------------------------------
p_grain          = Calc_p(MODEL,phi);
p_phase          = Collapse_p_By_Phase(p_grain,MODEL.phase_index,Nphase);

e_phase          = Calc_e(pars_phase,c_phase);
E_phase          = Calc_E_Tot(e_phase,p_phase);
Ne               = length(E_phase);

mu_e             = repmat({zeros(ny,nx)},1,Ne);
chi              = repmat({zeros(ny,nx)},Ne,Ne);
eta              = eta0*ones(ny,nx);

PARAM_LE         = struct();
PARAM_LE.Np      = Nphase;
PARAM_LE.Ne      = Ne;
PARAM_LE.eta     = eta;

MODEL_LE         = MODEL;
MODEL_LE.pars    = pars_phase;
MODEL_LE.phase_index = 1:Nphase;


%Calculation

tol_pure         = 1e-12;
p_pair           = 0.05;

for ip = 1:Nphase

    mask = p_phase(:,:,ip) > 1 - tol_pure;

    if ~any(mask(:))
        continue
    end

    Nmask = nnz(mask);

    for io = setdiff(1:Nphase,ip)

        fprintf('Checking pair initialization: phase %d %d\n',ip,io)

        pp          = p_phase;

        pp_ip       = pp(:,:,ip);
        pp_io       = pp(:,:,io);

        pp_ip(mask) = 1 - p_pair;
        pp_io(mask) = p_pair;

        pp(:,:,ip)  = pp_ip;
        pp(:,:,io)  = pp_io;

        e_test      = Calc_e(pars_phase,c_phase);
        E_test      = Calc_E_Tot(e_test,pp);

        p_slice     = zeros(1,Nmask,Nphase);

        for jp = 1:Nphase
            tmp             = pp(:,:,jp);
            p_slice(1,:,jp) = tmp(mask).';
        end

        c_slice     = Slice_c_Local(c_phase,mask);
        E_slice     = Slice_E_Local(E_test,mask);
        eta_slice   = reshape(eta(mask).',1,Nmask);

        mu_slice    = repmat({zeros(1,Nmask)},1,Ne);
        chi_slice   = repmat({zeros(1,Nmask)},Ne,Ne);

        STATE_LOC          = struct();
        STATE_LOC.c        = c_slice;
        STATE_LOC.e        = Calc_e(pars_phase,c_slice);
        STATE_LOC.E        = E_slice;
        STATE_LOC.mu_e     = mu_slice;
        STATE_LOC.chi      = chi_slice;
        STATE_LOC.omg      = zeros(1,Nmask,Nphase);
        STATE_LOC.phi      = [];
        STATE_LOC.p        = p_slice;
        STATE_LOC.mask     = ones(1,Nmask,Nphase);
        STATE_LOC.LE_state = [];

        PARAM_LOC          = PARAM_LE;
        PARAM_LOC.eta      = eta_slice;

        STATE_LOC          = LE_Run(STATE_LOC,PARAM_LOC,MODEL_LE);

        for ic = 1:length(STATE_LOC.c{io})
            tmp             = c_phase{io}{ic};
            tmp(mask)       = STATE_LOC.c{io}{ic}(:);
            c_phase{io}{ic} = tmp;
        end

        for ic = 1:length(STATE_LOC.c{ip})
            tmp             = c_phase{ip}{ic};
            tmp(mask)       = STATE_LOC.c{ip}{ic}(:);
            c_phase{ip}{ic} = tmp;
        end

    end
end

%% ------------------------------------------------------------------------
%  Expand initialized thermodynamic phase compositions to grains
% -------------------------------------------------------------------------
c                = Expand_c_By_Phase(c_phase,MODEL.phase_index);

p                = Calc_p(MODEL,phi);
e                = Calc_e(pars,c);
E                = Calc_E_Tot(e,p);

Ne               = length(E);
mu_e             = repmat({zeros(ny,nx)},1,Ne);
chi              = repmat({zeros(ny,nx)},Ne,Ne);

%% ------------------------------------------------------------------------
%  Final PARAM and STATE
% -------------------------------------------------------------------------
PARAM            = struct();
PARAM.Np         = Np;
PARAM.Ne         = Ne;
PARAM.eta        = eta;

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
%  Save structured + legacy variables
% -------------------------------------------------------------------------
save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'pars','pars_phase','c_phase','phs_name', ...
    'E_sc','Lx','Ly','c','E','e','p','phi','eta','mu_e','chi', ...
    'x','dx','nx','y','dy','ny','L_sc','F','Np','Ne','Nphase','Ngrain', ...
    'phase_ID','grain_ID','grain_phase','seed_xy')


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
        tmp             = c{ip}{ic};
        c_slice{ip}{ic} = reshape(tmp(mask),1,Nmask);
    end

end

end


function E_slice = Slice_E_Local(E,mask)

Ne      = length(E);
Nmask   = nnz(mask);
E_slice = cell(size(E));

for ie = 1:Ne
    tmp         = E{ie};
    E_slice{ie} = reshape(tmp(mask),1,Nmask);
end

end


function p_phase = Collapse_p_By_Phase(p_grain,phase_index,Nphase)

[ny,nx,~] = size(p_grain);
p_phase   = zeros(ny,nx,Nphase);

for iph = 1:Nphase
    grains = find(phase_index == iph);
    if ~isempty(grains)
        p_phase(:,:,iph) = sum(p_grain(:,:,grains),3);
    end
end

end


function c = Expand_c_By_Phase(c_phase,phase_index)

Ngrain = numel(phase_index);
c      = cell(1,Ngrain);

for ig = 1:Ngrain
    iph   = phase_index(ig);
    c{ig} = c_phase{iph};
end

end

function [phi,phase_ID,grain_ID,seed_xy,grain_phase] = Make_RandomPolygon_GrainMap(GRID,Nphase,opt)
%MAKE_RANDOMPOLYGON_GRAINMAP
%
% Generate polygonal grains.
% phi(:,:,ig) is one order parameter per grain.
% grain_phase(ig) gives the thermodynamic phase index of grain ig.
%
% If opt.separate_same_phase = true, adjacent grains are assigned
% different thermodynamic phases when possible.

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

if ~isfield(opt,'separate_same_phase') || isempty(opt.separate_same_phase)
    opt.separate_same_phase = true;
end

rng(opt.rng_seed)

x        = GRID.x(:).';
y        = GRID.y(:);
[X,Y]    = meshgrid(x,y);
[ny,nx]  = size(X);

xmin     = min(x);
xmax     = max(x);
ymin     = min(y);
ymax     = max(y);

Lx       = xmax - xmin;
Ly       = ymax - ymin;

if isempty(opt.n_grain)
    area_box = max(Lx*Ly,eps);
    area_g   = max(opt.grain_size^2,eps);
    Ngrain   = max(Nphase,round(area_box/area_g));
else
    Ngrain   = max(Nphase,round(opt.n_grain));
end

% -------------------------------------------------------------------------
% Random grain seeds
% -------------------------------------------------------------------------
sx       = xmin + Lx*rand(Ngrain,1);
sy       = ymin + Ly*rand(Ngrain,1);
seed_xy  = [sx,sy];

% -------------------------------------------------------------------------
% Voronoi-like nearest-seed grain map
% -------------------------------------------------------------------------
grain_ID = zeros(ny,nx);
bestD2   = inf(ny,nx);

for ig = 1:Ngrain

    D2 = (X - sx(ig)).^2 + (Y - sy(ig)).^2;

    mask = D2 < bestD2;

    bestD2(mask)   = D2(mask);
    grain_ID(mask) = ig;

end

% -------------------------------------------------------------------------
% Grain adjacency graph
% -------------------------------------------------------------------------
adj = false(Ngrain,Ngrain);

A = grain_ID(:,1:end-1);
B = grain_ID(:,2:end);
pair = [A(:),B(:)];
pair = pair(pair(:,1) ~= pair(:,2),:);

for k = 1:size(pair,1)
    i = pair(k,1);
    j = pair(k,2);
    adj(i,j) = true;
    adj(j,i) = true;
end

A = grain_ID(1:end-1,:);
B = grain_ID(2:end,:);
pair = [A(:),B(:)];
pair = pair(pair(:,1) ~= pair(:,2),:);

for k = 1:size(pair,1)
    i = pair(k,1);
    j = pair(k,2);
    adj(i,j) = true;
    adj(j,i) = true;
end

% -------------------------------------------------------------------------
% Assign thermodynamic phase to each grain
% -------------------------------------------------------------------------
if opt.separate_same_phase
    grain_phase = Assign_Phases_Separated(adj,Nphase,opt.phase_mode,opt.phase_fraction);
else
    grain_phase = Assign_Phases_Basic(Ngrain,Nphase,opt.phase_mode,opt.phase_fraction);
end

% -------------------------------------------------------------------------
% Make phase_ID map
% -------------------------------------------------------------------------
phase_ID = zeros(ny,nx);

for ig = 1:Ngrain
    phase_ID(grain_ID == ig) = grain_phase(ig);
end

% -------------------------------------------------------------------------
% One phi field per grain
% -------------------------------------------------------------------------
phi = zeros(ny,nx,Ngrain);

for ig = 1:Ngrain
    phi(:,:,ig) = double(grain_ID == ig);
end

end


function grain_phase = Assign_Phases_Basic(Ngrain,Nphase,phase_mode,phase_fraction)

grain_phase = zeros(Ngrain,1);

switch lower(phase_mode)

    case 'cycle'

        for ig = 1:Ngrain
            grain_phase(ig) = mod(ig-1,Nphase) + 1;
        end

    case 'fractions'

        frac      = phase_fraction(:).';
        frac      = frac / sum(frac);
        edge      = [0,cumsum(frac)];
        edge(end) = 1;
        r         = rand(Ngrain,1);

        for ig = 1:Ngrain
            grain_phase(ig) = find(r(ig) >= edge(1:end-1) & r(ig) <= edge(2:end),1,'first');
        end

    case 'random'

        grain_phase = randi(Nphase,Ngrain,1);

    otherwise

        error('Unknown phase_mode: %s',phase_mode)

end

if Ngrain >= Nphase
    grain_phase(1:Nphase) = 1:Nphase;
    grain_phase = grain_phase(randperm(Ngrain));
end

end


function grain_phase = Assign_Phases_Separated(adj,Nphase,phase_mode,phase_fraction)
%ASSIGN_PHASES_SEPARATED
%
% Greedy graph coloring / phase assignment.
% Adjacent grains avoid the same thermodynamic phase when possible.

Ngrain = size(adj,1);

phase_fraction = phase_fraction(:).';
phase_fraction = phase_fraction / sum(phase_fraction);

target = max(1,round(phase_fraction*Ngrain));
target = target * Ngrain / sum(target);

grain_phase = zeros(Ngrain,1);
count_phase = zeros(1,Nphase);

deg = sum(adj,2);
[~,ord] = sort(deg + 0.01*rand(Ngrain,1),'descend');

for kk = 1:Ngrain

    ig = ord(kk);

    nb = find(adj(ig,:));
    used = unique(grain_phase(nb));
    used = used(used > 0);

    cand = setdiff(1:Nphase,used);

    if isempty(cand)
        cand = 1:Nphase;
    end

    switch lower(phase_mode)

        case 'cycle'

            preferred = mod(ig-1,Nphase) + 1;

            if ismember(preferred,cand)
                iph = preferred;
            else
                iph = cand(1);
            end

        case 'fractions'

            score = count_phase(cand) ./ target(cand);
            score = score + 1e-6*rand(size(score));
            [~,id] = min(score);
            iph = cand(id);

        case 'random'

            iph = cand(randi(numel(cand)));

        otherwise

            error('Unknown phase_mode: %s',phase_mode)

    end

    grain_phase(ig) = iph;
    count_phase(iph) = count_phase(iph) + 1;

end

% -------------------------------------------------------------------------
% Try to make sure every phase appears at least once
% without creating same-phase contacts.
% -------------------------------------------------------------------------
missing = setdiff(1:Nphase,unique(grain_phase).');

for iph = missing

    success = false;

    [~,try_ord] = sort(count_phase(grain_phase),'descend');

    for kk = 1:Ngrain

        ig = try_ord(kk);

        if count_phase(grain_phase(ig)) <= 1
            continue
        end

        nb = find(adj(ig,:));

        if any(grain_phase(nb) == iph)
            continue
        end

        old = grain_phase(ig);

        grain_phase(ig) = iph;
        count_phase(old) = count_phase(old) - 1;
        count_phase(iph) = count_phase(iph) + 1;

        success = true;
        break

    end

    if ~success
        warning('Could not insert phase %d without same-phase contact.',iph)
    end

end

end