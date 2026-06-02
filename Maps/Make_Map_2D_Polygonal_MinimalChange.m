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
eta0             = 5000e10/E_sc;      %Use the same eta in Run_2D

%% ------------------------------------------------------------------------
%  Thermodynamic phases and Gibbs-equilibrium starting values
%  Edit only this input block when changing the phase assemblage.
% -------------------------------------------------------------------------
phs_name         = {'Grt','Olv','Cpx','Fel'};
pars_phase       = Load_Data(phs_name);
Nphase           = length(pars_phase);

%Phase proportions from Gibbs minimization, in the same order as phs_name
phase_prop       = [0.30 0.30 0.20 0.2];
phase_prop       = phase_prop/sum(phase_prop);

%Independent endmember proportions from Gibbs minimization.
c_value          = cell(1,Nphase);
c_value{1}       = [0.38 0.51];                  %Grt
c_value{2}       = [0.25 0.35 0.30 ];            %Olv
c_value{3}       = [0.11 0.26  0 0.02 0.5];      %Cpx
c_value{4}       = [0.45 ];                      %Fel
c_guess          = cell(1,Nphase);

for ip = 1:Nphase
    Nc           = size(pars_phase{ip}.n,1)-1;
    c_guess{ip}  = num2cell(c_value{ip});
end

%% ------------------------------------------------------------------------
%  One reference LE solve at the intended bulk composition
% -------------------------------------------------------------------------
p_ref            = reshape(phase_prop,1,1,Nphase);
e_guess          = Calc_e(pars_phase,c_guess);
E_target         = Calc_E_Tot(e_guess,p_ref);
Ne               = length(E_target);
%LE calculation
[c_ref,mu_ref]   = LE_Calculator(pars_phase,p_ref,c_guess,E_target,eta0,[0.1,2000]);

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

%% ------------------------------------------------------------------------
%  Grid and random polygonal-grain controls
% -------------------------------------------------------------------------
Lx               = 5e-6;
Ly               = 5e-6;

nx               = 80;
ny               = 80;

%Approximate mean equivalent grain diameter.
%Ngrain is calculated from total map area / mean grain area.
grain_size       = 3.0e-6;
rng_seed         = 5;

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
%  Two-dimensional random polygonal grains
%
%  Only the grain geometry is changed relative to the working pseudo-1D
%  initializer.  The compositions c_ref and E_offset above remain those
%  from the same reference LE solution.  Phase fractions are approximate:
%  complete polygonal grains are assigned to phases.
% -------------------------------------------------------------------------
grain_area       = pi/4*(grain_size/L_sc)^2;
Ngrain           = max(Nphase,round(GRID.Lx*GRID.Ly/grain_area));
Np               = Ngrain;

rng(rng_seed,'twister')

[grain_ID,seed_xy] = Generate_Polygonal_Grains_2D(GRID,Ngrain);
[grain_ID,seed_xy] = Remove_Empty_Grains(grain_ID,seed_xy);

Ngrain           = size(seed_xy,1);
Np               = Ngrain;
grain_phase      = Assign_Phases_To_Grains(grain_ID,phase_prop);
grain_phase      = grain_phase(:);

phi              = zeros(ny,nx,Ngrain);
phase_ID         = zeros(ny,nx);

for ig = 1:Ngrain

    mask              = grain_ID == ig;
    phi(:,:,ig)       = mask;
    phase_ID(mask)    = grain_phase(ig);

end

%% ------------------------------------------------------------------------
%  Grain-resolved MODEL
% -------------------------------------------------------------------------
MODEL                    = struct();
MODEL.phs_name           = phs_name;
MODEL.phase_index        = grain_phase(:).';
MODEL.pars               = pars_phase(MODEL.phase_index);

eps_phi                  = 0;
MODEL.p_fun              = @(a,phi) phi(:,:,a).^2 ./ (sum(phi.^2,3) + eps_phi);
MODEL.dpdphi             = @(a,b,phi) (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) - 2*phi(:,:,a).*phi(:,:,b).^2 ./ (sum(phi.^2,3) + eps_phi).^2;

pars                     = MODEL.pars;
F                        = MODEL;

figure(3); clf
imagesc(x*1e6,y*1e6,phase_ID)
set(gca,'YDir','normal')
axis image tight
colorbar
title('Initial 2-D random polygonal phase map')
xlabel('x \mum');ylabel('y \mum')


%% ------------------------------------------------------------------------
%  Construct the penalty-consistent conserved composition field
% -------------------------------------------------------------------------
c                = Expand_c_By_Phase(c_phase,MODEL.phase_index);
p                = Calc_p(MODEL,phi);
e                = Calc_e(pars,c);
E                = Calc_E_Tot(e,p);
eta              = eta0*ones(ny,nx);

%The penalty formulation requires E - Emix = mu/eta at coexistence.
%Adding the same offset at every node makes each pure phase and each
%future diffuse mixture consistent with the same reference mu.
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

p_phase          = Collapse_p_By_Phase(STATE.p,MODEL.phase_index,Nphase);
phase_prop_map   = zeros(1,Nphase);

for ip = 1:Nphase
    tmp                = p_phase(:,:,ip);
    phase_prop_map(ip) = mean(tmp(:));
end

omega_phase      = zeros(ny,nx,Nphase);

for ip = 1:Nphase
    ig                  = find(MODEL.phase_index == ip,1,'first');
    omega_phase(:,:,ip) = STATE.omg(:,:,ig);
end

for ip = 1:Nphase-1
    for jp = ip+1:Nphase
        domg = omega_phase(:,:,ip) - omega_phase(:,:,jp);
        fprintf('%s - %s: max|domg| = %.8e, mean|domg| = %.8e\n', ...
                phs_name{ip},phs_name{jp}, max(abs(domg(:))),mean(abs(domg(:))))
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
save('Map2d.mat','PHYS','GRID','MODEL','PARAM','STATE', ...
     'E_sc','L_sc','eta','pars','Np','Ngrain','grain_ID','phase_ID', ...
     'grain_phase','seed_xy','grain_size','rng_seed','phase_prop_map')

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


function [grain_ID,seed_xy] = Generate_Polygonal_Grains_2D(GRID,Ngrain)
%GENERATE_POLYGONAL_GRAINS_2D Random Voronoi polygonal grain map.

[X,Y]    = meshgrid(GRID.x,GRID.y);
seed_xy  = [GRID.Lx*rand(Ngrain,1), GRID.Ly*rand(Ngrain,1)];

dist2    = inf(GRID.ny,GRID.nx);
grain_ID = zeros(GRID.ny,GRID.nx);

for ig = 1:Ngrain

    d2   = (X-seed_xy(ig,1)).^2+(Y-seed_xy(ig,2)).^2;
    mask = d2 < dist2;

    dist2(mask)    = d2(mask);
    grain_ID(mask) = ig;

end

end


function [grain_ID,seed_xy] = Remove_Empty_Grains(grain_ID,seed_xy)
%REMOVE_EMPTY_GRAINS Remove seeds without a represented grid cell.

grain_old = unique(grain_ID(:));

if length(grain_old) == size(seed_xy,1)
    return
end

id_new              = zeros(size(seed_xy,1),1);
id_new(grain_old)   = 1:length(grain_old);
grain_ID            = reshape(id_new(grain_ID(:)),size(grain_ID));
seed_xy             = seed_xy(grain_old,:);

end


function grain_phase = Assign_Phases_To_Grains(grain_ID,phase_prop)
%ASSIGN_PHASES_TO_GRAINS Assign complete grains to approximate phase areas.

Ngrain              = max(grain_ID(:));
Nphase              = length(phase_prop);
grain_area          = accumarray(grain_ID(:),1,[Ngrain,1]);
target_area         = phase_prop(:).'*numel(grain_ID);
phase_area          = zeros(1,Nphase);
grain_phase         = zeros(Ngrain,1);

[~,igrain]          = sort(grain_area,'descend');

for ii = 1:Ngrain

    ig              = igrain(ii);
    err             = zeros(1,Nphase);

    for ip = 1:Nphase

        area_try     = phase_area;
        area_try(ip) = area_try(ip)+grain_area(ig);
        err(ip)      = sum((area_try-target_area).^2);

    end

    [~,ip]          = min(err);
    grain_phase(ig) = ip;
    phase_area(ip)  = phase_area(ip)+grain_area(ig);

end

%When there are enough grains, retain at least one grain for every phase.
for ip = find(phase_prop > 0 & accumarray(grain_phase,1,[Nphase,1]).' == 0)

    count           = accumarray(grain_phase,1,[Nphase,1]).';
    donor           = find(count > 1,1,'first');
    id_donor        = find(grain_phase == donor);
    [~,imin]        = min(grain_area(id_donor));
    grain_phase(id_donor(imin)) = ip;

end

end
