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
eta0             = 4000e10/E_sc;      %Use the same eta in Run_2D

%% ------------------------------------------------------------------------
%  Thermodynamic phases and Gibbs-equilibrium starting values
%  Edit only this input block when changing the phase assemblage.
% -------------------------------------------------------------------------
phs_name         = {'Grt','Cpx','Fel','Olv'};
pars_phase       = Load_Data(phs_name);
Nphase           = length(pars_phase);

%Phase proportions from Gibbs minimization, in the same order as phs_name
phase_prop       = [0.3 0.2 0.3 0.2];
phase_prop       = phase_prop/sum(phase_prop);

%Independent endmember proportions from Gibbs minimization.
c_value          = cell(1,Nphase);
c_value{1}       = [0.30 0.60 ];                 %Grt
c_value{2}       = [0.04 0.28 0 0.04 0.44];      %Cpx
c_value{3}       = [0.34];                       %Fel
c_value{4}       = [0.05 0.40 0.40];             %Olv
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
%  Grid
% -------------------------------------------------------------------------
Lx               = 5e-6;
Ly               = 5e-6;

nx               = 200;
ny               = 4;

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
%  Pseudo-2D repeating 1D phase bands along x
%
%  Default cyclic order:
%     Nphase = 2, Nrepeat = 2:  1 | 2 | 1 | 2
%     Nphase = 3, Nrepeat = 2:  1 | 2 | 3 | 1 | 2 | 3
%
%  band_order may be edited manually to choose particular interfaces.
%  Each phase keeps its total requested phase_prop even when it appears
%  multiple times in band_order.
% -------------------------------------------------------------------------
Nrepeat          = 1;
band_order       = repmat(1:Nphase,1,Nrepeat);

Noccur           = accumarray(band_order(:),1,[Nphase,1]).';

band_fraction    = phase_prop(band_order)./Noccur(band_order);
grain_phase      = band_order;
Ngrain           = numel(grain_phase);
Np               = Ngrain;

x_edges          = 1 + round(nx*[0 cumsum(band_fraction)]);
x_edges(end)     = nx + 1;

phi              = zeros(ny,nx,Ngrain);
phase_ID         = zeros(ny,nx);
grain_ID         = zeros(ny,nx);
seed_xy          = zeros(Ngrain,2);

for ig = 1:Ngrain

    iph                      = grain_phase(ig);
    ix1                      = x_edges(ig);
    ix2                      = x_edges(ig+1)-1;

    phi(:,ix1:ix2,ig)       = 1;
    phase_ID(:,ix1:ix2)     = iph;
    grain_ID(:,ix1:ix2)     = ig;

    seed_xy(ig,1)           = mean(x(ix1:ix2));
    seed_xy(ig,2)           = mean(y);
end

grain_phase      = grain_phase(:);

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
pcolor(x*1e6,y*1e6,phase_ID)
shading flat;axis equal tight;colorbar
title('Initial repeating 1D thermodynamic phase bands')
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
save('Map2d.mat','PHYS','GRID','MODEL','PARAM','STATE', 'E_sc','L_sc','eta', 'pars', 'Np')

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
