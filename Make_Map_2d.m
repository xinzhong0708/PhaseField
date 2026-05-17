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

eta0             = 1000e10/E_sc;

%% ------------------------------------------------------------------------
%  Thermodynamic data
% -------------------------------------------------------------------------
pars             = Load_Data({'Olv','Cpx','Grt','Qtz','Crd','Spl'});
% pars             = Load_Data({'Qtz','Crd'});

Np               = length(pars);

%% ------------------------------------------------------------------------
%  Model functions
% -------------------------------------------------------------------------
MODEL            = struct();
MODEL.pars       = pars;

eps_phi          = 1e-14;

MODEL.p_fun      = @(a,phi) ...
    phi(:,:,a).^2 ./ (sum(phi.^2,3) + eps_phi);

MODEL.dpdphi     = @(a,b,phi) ...
    (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) ...
    - 2*phi(:,:,a).*phi(:,:,b).^2 ./ (sum(phi.^2,3) + eps_phi).^2;

% legacy name, for older functions
F                = MODEL;

%% ------------------------------------------------------------------------
%  Grid
% -------------------------------------------------------------------------
Lx               = 5e-6;
Ly               = 5e-6;

nx               = 60*6;
ny               = 4;

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

% Spl
c{6}{1}          = 0.5000*ones(ny,nx);
c{6}{2}          = 0.1000*ones(ny,nx);
c{6}{3}          = 0.1000*ones(ny,nx);



% % Qtz
% c{1}{1}          = 1.0000*ones(ny,nx);
% 
% % Crd
% c{2}{1}          = 1.0000*ones(ny,nx);



%% ------------------------------------------------------------------------
%  Initial phase field
% -------------------------------------------------------------------------
phi              = zeros(ny,nx,Np);

md               = nx/6;
phi(:,1:md,1)              = 1;
phi(:,1*md+1:2*md,2)       = 1;
phi(:,2*md+1:3*md,3)       = 1;
phi(:,3*md+1:4*md,4)       = 1;
phi(:,4*md+1:5*md,5)       = 1;
phi(:,5*md+1:6*md,6)       = 1;

% md               = nx/2;
% phi(:,1:md,1)              = 1;
% phi(:,1*md+1:2*md,2)       = 1;
 

% md               = nx/2;
% phi(:,1:md,1)              = 1;
% phi(:,1*md+1:2*md,2)       = 1;
% phi(:,2*md+1:3*md,3)       = 1;
% phi(:,3*md+1:4*md,4)       = 1;


% md               = nx/2;
% phi(   1:md,    1:md  ,1)  =  1;
% phi(md+1:nx,    1:md  ,2)  =  1;
% phi(1:md-10, md+1:nx  ,3)  =  1;
% phi(md-9:nx, md+1:nx  ,4)  =  1;


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