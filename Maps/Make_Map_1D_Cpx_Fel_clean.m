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
eta0             = 1000e10/E_sc;

%% ------------------------------------------------------------------------
%  Thermodynamic phases and target bulk equilibrium
% -------------------------------------------------------------------------
phs_name         = {'Grt','Cpx','Fel'};
pars_phase       = Load_Data(phs_name);
Nphase           = length(pars_phase);

%Phase proportions expected from Gibbs minimization.
%Replace this line with your linprog phase proportions when needed.
phase_prop       = [1/3 1/3 1/3];
phase_prop       = phase_prop/sum(phase_prop);

%Independent endmember proportions from Gibbs minimization.
%These are used as starting compositions for the finite-eta LE solve.
c_guess          = cell(1,Nphase);

% Grt
c_guess{1}{1}    =  0.42;
c_guess{1}{2}    =  0.48;

%Cpx
c_guess{2}{1}    =  0.11;
c_guess{2}{2}    =  0.26;
c_guess{2}{3}    = -0.00;
c_guess{2}{4}    =  0.00;
c_guess{2}{5}    =  0.52;

%Feldspar
c_guess{3}{1}    =  0.48;

%% ------------------------------------------------------------------------
%  One reference LE solve at the intended bulk composition
% -------------------------------------------------------------------------
p_ref            = reshape(phase_prop,1,1,Nphase);
e_guess          = Calc_e(pars_phase,c_guess);
E_target         = Calc_E_Tot(e_guess,p_ref);
Ne               = length(E_target);

[c_ref,mu_ref]   = LE_Calculator(pars_phase,p_ref,c_guess,E_target,eta0,[0.2,1000]);

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
        omega_ref(ip) = omega_ref(ip) - e_ref{ip}{ie} .* mu_ref{ie};
    end

end

fprintf('\nReference LE state at prescribed phase proportions:\n')
for ip = 1:Nphase
    fprintf('%s endmember proportion:\n',phs_name{ip})
    disp(cell2mat(c_ref{ip}(:)).')
end
fprintf('max|omega - omega(1)| = %.8e\n',max(abs(omega_ref-omega_ref(1))))

fprintf('\nFinite-eta E offset: calculated residual and mu/eta should agree:\n')
for ie = 1:Ne
    fprintf('element %d: E_offset = %.8e, mu/eta = %.8e\n', ...
            ie,E_offset{ie},mu_ref{ie}/eta0)
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
%  Nrepeat = 1: Cpx | Fel; Nrepeat = 2: Cpx | Fel | Cpx | Fel.
% -------------------------------------------------------------------------
Nrepeat          = 1;
grain_phase      = repmat(1:Nphase,1,Nrepeat);
band_fraction    = repmat(phase_prop/Nrepeat,1,Nrepeat);

Ngrain           = numel(grain_phase);
Np               = Ngrain;

phi              = zeros(ny,nx,Ngrain);
phase_ID         = zeros(ny,nx);
grain_ID         = zeros(ny,nx);
seed_xy          = zeros(Ngrain,2);

x_edges          = 1 + round(nx*[0 cumsum(band_fraction)]);
x_edges(end)     = nx + 1;

for ig = 1:Ngrain

    iph          = grain_phase(ig);
    ix1          = x_edges(ig);
    ix2          = x_edges(ig+1)-1;

    phi(:,ix1:ix2,ig)     = 1;
    phase_ID(:,ix1:ix2)   = iph;
    grain_ID(:,ix1:ix2)   = ig;

    seed_xy(ig,1)         = mean(x(ix1:ix2));
    seed_xy(ig,2)         = mean(y);

end

grain_phase      = grain_phase(:);

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
title('Initial repeating 1D thermodynamic phase bands')
xlabel('x \mum')
ylabel('y \mum')
drawnow

%% ------------------------------------------------------------------------
%  Construct penalty-consistent conserved composition field
% -------------------------------------------------------------------------
c                = Expand_c_By_Phase(c_phase,MODEL.phase_index);
p                = Calc_p(MODEL,phi);
e                = Calc_e(pars,c);
eta              = eta0*ones(ny,nx);
E                = Calc_E_Tot(e,p);

%Add the same finite-eta residual everywhere. This prevents the first LE
%call from shifting c simply because the penalty formulation needs
%E - Emix = mu/eta at coexistence.
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
                max(abs(domg(:))),mean(abs(domg(:))));
    end
end

E_mean           = zeros(Ne,1);
E_target_vec     = zeros(Ne,1);

for ie = 1:Ne
    E_mean(ie)       = mean(STATE.E{ie}(:));
    E_target_vec(ie) = E_target{ie};
end

fprintf('max|mean(E map) - E target| = %.8e\n',max(abs(E_mean-E_target_vec)))

%Keep legacy variables consistent with STATE after LE evaluation.
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
    'pars','pars_phase','c_phase','c_guess','c_ref','mu_ref','E_target', ...
    'E_offset','phase_prop','phs_name','E_sc','Lx','Ly','c','E','e','p', ...
    'phi','eta','mu_e','chi','x','dx','nx','y','dy','ny','L_sc','F', ...
    'Np','Ne','Nphase','Ngrain','phase_ID','grain_ID','grain_phase', ...
    'seed_xy','Nrepeat')

%% ========================================================================
%  Local helper function
% ========================================================================
function c = Expand_c_By_Phase(c_phase,phase_index)

Ngrain = numel(phase_index);
c      = cell(1,Ngrain);

for ig = 1:Ngrain
    c{ig} = c_phase{phase_index(ig)};
end

end
