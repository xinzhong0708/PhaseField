%% Create_Map2D_OlivineMelt_Metadata.m
% Generate Map2d.mat for olivine seeds in melt.
%
% Active phases are fixed here:
%
%   phs_name = {'Olivine','Melt'}
%
% Thermodynamic / numerical metadata are read from:
%
%   MetaData.xlsx
%
% Grain orientation is saved as:
%
%   PARAM.theta_grain(ig)
%
% Olivine grains receive random orientations in [0,pi).
% Melt orientation is set to zero.

clear; clf

addpath('..\bin')
addpath('..\')
addpath('..\Thermo')
addpath('..\Thermo\Solutions')

%% ========================================================================
%  User controls
% ========================================================================

metadata_file = 'MetaData.xlsx';

% Active phases for this map
phs_name = {'Olivine','Melt'};

% Domain
Lx = 2e-5;
Ly = 2e-5;
nx = 90;
ny = 90;

% Olivine seeds in melt
N_olivine_seed = 1;
olivine_frac   = 0.005;
min_gap_fac    = 3.50;      % minimum center distance relative to 2*radius
rng_seed       = 1;

% Grain orientation controls
theta_rng_seed   = 1002;
theta_grain_user = [];

%% ========================================================================
%  Read metadata
% ========================================================================

META   = Read_Map_Metadata(metadata_file,phs_name);

T      = META.T;
P      = META.P;
Cname  = META.Cname;
solmod = META.solmod;

PHYS        = struct();
PHYS.E_sc   = META.E_sc;
PHYS.L_sc   = META.L_sc;
PHYS.t_sc   = META.t_sc;
PHYS.vref   = META.vref;
PHYS.M_phs  = META.M_phs;

E_sc        = PHYS.E_sc;
L_sc        = PHYS.L_sc;
vref        = PHYS.vref;
eta0        = 4000e10/E_sc;

phase_prop = [olivine_frac,1-olivine_frac];

% Initial independent endmember compositions.
% Prefer values from phase sheets:
%   Initial c    value1 value2 ...
%
% Fallback values are only placeholders.
c_value = cell(1,numel(phs_name));
c_value{1} = Get_Phase_Init_C(META,'Olivine',[0.45 0.45]);
c_value{2} = Get_Phase_Init_C(META,'Melt'   ,[0.35]);

%% ========================================================================
%  Build thermodynamics and reference LE
% ========================================================================

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


%% ========================================================================
%  Grid
% ========================================================================

x_SI  = linspace(0,Lx,nx);
y_SI  = linspace(0,Ly,ny);
dx_SI = x_SI(2)-x_SI(1);
dy_SI = y_SI(2)-y_SI(1);

x  = x_SI/L_sc;
y  = y_SI/L_sc;
dx = dx_SI/L_sc;
dy = dy_SI/L_sc;

GRID       = struct();
GRID.x     = x;
GRID.y     = y;
GRID.dx    = dx;
GRID.dy    = dy;
GRID.nx    = nx;
GRID.ny    = ny;
GRID.Lx    = Lx/L_sc;
GRID.Ly    = Ly/L_sc;

% Store SI copy for checking / postprocessing
GRID.x_SI  = x_SI;
GRID.y_SI  = y_SI;
GRID.dx_SI = dx_SI;
GRID.dy_SI = dy_SI;
GRID.Lx_SI = Lx;
GRID.Ly_SI = Ly;

% This tells Read_PFM_Metadata not to scale GRID again.
GRID.is_scaled = 1;

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
%  Generate olivine seeds in melt
% ========================================================================

rng(rng_seed,'twister')

[X,Y] = meshgrid(GRID.x,GRID.y);

ol_area       = olivine_frac*GRID.Lx*GRID.Ly;
seed_radius   = sqrt(ol_area/(N_olivine_seed*pi));
seed_radius_m = seed_radius*L_sc;

seed_xy       = zeros(N_olivine_seed+1,2);

for iseed = 1:N_olivine_seed

    accepted = 0;

    for itry = 1:10000

        trial_xy = [seed_radius + (GRID.Lx-2*seed_radius)*rand, ...
                    seed_radius + (GRID.Ly-2*seed_radius)*rand];

        if iseed == 1

            accepted = 1;

        else

            d2 = sum((seed_xy(1:iseed-1,:) - trial_xy).^2,2);
            accepted = all(sqrt(d2) > min_gap_fac*2*seed_radius);

        end

        if accepted == 1
            seed_xy(iseed,:) = trial_xy;
            break
        end
    end

    if accepted == 0
        error('Could not place non-overlapping olivine seeds. Reduce N_olivine_seed or olivine_frac.')
    end
end

% Melt grain uses domain center only as a marker
melt_grain              = N_olivine_seed + 1;
seed_xy(melt_grain,:)   = [mean(GRID.x),mean(GRID.y)];

seed_xy(1,1) = GRID.Lx/2;
seed_xy(1,2) = GRID.Ly/2;

Ngrain = N_olivine_seed + 1;
Np     = Ngrain;

grain_phase = zeros(Ngrain,1);
grain_phase(1:N_olivine_seed) = 1;     % Olivine
grain_phase(melt_grain)       = 2;     % Melt

grain_ID = melt_grain*ones(ny,nx);
phase_ID = 2*ones(ny,nx);

for ig = 1:N_olivine_seed

    mask = (X-seed_xy(ig,1)).^2 + (Y-seed_xy(ig,2)).^2 <= seed_radius^2;

    grain_ID(mask) = ig;
    phase_ID(mask) = 1;

end

phi = zeros(ny,nx,Ngrain);

for ig = 1:Ngrain
    phi(:,:,ig) = grain_ID == ig;
end

% Recompute realised phase fraction after rasterization
phase_prop_geom = Compute_Phase_Fraction_From_Map(phase_ID,Nphase);
phase_prop_map  = phase_prop_geom;

%% ========================================================================
%  Grain orientations
% ========================================================================

theta_grain = zeros(1,Ngrain);

if isempty(theta_grain_user)

    rng_state = rng;

    if ~isempty(theta_rng_seed)
        rng(theta_rng_seed,'twister')
    end

    theta_grain(1:N_olivine_seed) = 0*pi*rand(1,N_olivine_seed);
    theta_grain(melt_grain)       = 0;

    rng(rng_state)

else

    if numel(theta_grain_user) ~= Ngrain
        error('theta_grain_user must be [] or length Ngrain.')
    end

    theta_grain = theta_grain_user(:).';

end

%% ========================================================================
%  Grain-resolved MODEL
% ========================================================================

MODEL             = struct();
MODEL.phs_name    = phs_name;
MODEL.phase_index = grain_phase(:).';
MODEL.pars        = pars_phase(MODEL.phase_index);

MODEL.Cname       = Cname;
MODEL.solmod      = solmod;
MODEL.T           = T;
MODEL.P           = P;
MODEL.E_sc        = E_sc;
MODEL.vref        = vref;
MODEL.eta         = eta0;

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
title('Initial olivine seeds in melt')
xlabel('x (\mum)')
ylabel('y (\mum)')
hold on
plot(seed_xy(1:N_olivine_seed,1)*L_sc*1e6, ...
     seed_xy(1:N_olivine_seed,2)*L_sc*1e6,'k.','MarkerSize',8)
hold off
drawnow

fprintf('\nInitial map: olivine seeds in melt\n')
fprintf('Number of olivine grains = %d\n',N_olivine_seed)
fprintf('Olivine seed radius      = %.4e m\n',seed_radius_m)
fprintf('\nPhase       requested        realised area fraction      difference\n')

for ip = 1:Nphase
    fprintf('%-12s %.8f         %.8f                 %+.3e\n', ...
        phs_name{ip},phase_prop(ip),phase_prop_geom(ip), ...
        phase_prop_geom(ip)-phase_prop(ip))
end

fprintf('\nGrain orientations theta_grain:\n')
disp(theta_grain)

%% ========================================================================
%  Construct locally penalty-consistent conserved composition field
% ========================================================================

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
chi  = repmat({zeros(ny,nx)},Ne,Ne);

%% ========================================================================
%  PARAM and STATE
% ========================================================================

PARAM            = struct();
PARAM.Np         = Np;
PARAM.Ne         = Ne;
PARAM.eta        = eta;
PARAM.use_WScale = 0;
PARAM.T          = T;
PARAM.P          = P;
PARAM.theta_grain = theta_grain;

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
%  Completeness check before LE
% ========================================================================

if numel(MODEL.phase_index) ~= Ngrain
    error('MODEL.phase_index length does not match Ngrain.')
end

if numel(MODEL.pars) ~= Ngrain
    error('MODEL.pars length does not match Ngrain.')
end

if size(STATE.phi,3) ~= Ngrain
    error('STATE.phi third dimension does not match Ngrain.')
end

if size(STATE.p,3) ~= Ngrain
    error('STATE.p third dimension does not match Ngrain.')
end

if numel(PARAM.theta_grain) ~= Ngrain
    error('PARAM.theta_grain length does not match Ngrain.')
end

if PARAM.Np ~= Ngrain
    error('PARAM.Np does not match Ngrain.')
end

if PARAM.Ne ~= Ne
    error('PARAM.Ne does not match Ne.')
end

%% ========================================================================
%  Initial LE check
% ========================================================================

STATE_INI = STATE;

PARAM.LE_mode = 'LE';
STATE         = LE_Run_Mode_New(STATE,PARAM,MODEL);

p_phase        = Collapse_p_By_Phase(STATE.p,MODEL.phase_index,Nphase);
phase_prop_map = zeros(1,Nphase);

for ip = 1:Nphase
    tmp = p_phase(:,:,ip);
    phase_prop_map(ip) = mean(tmp(:));
end

omega_phase = zeros(ny,nx,Nphase);

for ip = 1:Nphase
    grains = find(MODEL.phase_index == ip);

    if ~isempty(grains)
        omega_phase(:,:,ip) = mean(STATE.omg(:,:,grains),3);
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
phi  = STATE.phi;

%% ========================================================================
%  Save
% ========================================================================

save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'E_sc','L_sc','vref','eta','pars','Np','Ne','Nphase','Ngrain', ...
    'T','P','Cname','solmod', ...
    'phs_name','phase_prop','phase_prop_geom','phase_prop_map', ...
    'grain_ID','phase_ID','grain_phase','seed_xy','theta_grain', ...
    'theta_grain_user','theta_rng_seed', ...
    'N_olivine_seed','olivine_frac','seed_radius','rng_seed', ...
    'E_target','E_offset','E_bulk_shift','c_ref','mu_ref','omega_ref', ...
    'c','e','E','p','mu_e','chi','phi')

fprintf('\nSaved Map2d.mat\n')
fprintf('Saved full-name phases:\n')
disp(MODEL.phs_name)
fprintf('Saved PARAM.theta_grain:\n')
disp(PARAM.theta_grain)

%% ========================================================================
%  Local helper functions
% ========================================================================

function META = Read_Map_Metadata(xlsx_file,phs_name)

Cmain = readcell(xlsx_file,'Sheet','Main');
Cpt   = readcell(xlsx_file,'Sheet','PTt_Path');

META = struct();

elem_list = Get_Row_List(Cmain,'Independent c (Si dependent)');

META.Cname  = elem_list;
META.solmod = Get_Value(Cmain,'Solution name','solution_models_PFM');
META.vref   = Get_Value(Cmain,'Molar volume',2e-5);
META.E_sc   = Get_Value(Cmain,'Energy scale',1e9);
META.L_sc   = Get_Value(Cmain,'Length scale',1);
META.t_sc   = Get_Value(Cmain,'Time scale',1);

if isempty(META.Cname)
    error('Read_Map_Metadata: element list is missing in Main sheet.')
end

if ~any(strcmp(META.Cname,'Si'))
    META.Cname{end+1} = 'Si';
end

if ~any(strcmp(META.Cname,'O'))
    META.Cname{end+1} = 'O';
end

[META.T,META.P] = Read_First_PT(Cpt);

Nphase = numel(phs_name);
Ne     = numel(elem_list);

META.M_phs = zeros(Nphase,Ne);

for ip = 1:Nphase

    Cph = readcell(xlsx_file,'Sheet',phs_name{ip});
    META.phase(ip).name             = phs_name{ip};
    META.phase(ip).interface_energy = Get_Value(Cph,'Interface energy',0.5);
    META.phase(ip).aniso_mode       = Get_Value(Cph,'Anisotropy mode','iso');
    META.phase(ip).init_c           = Get_Row_List_Number(Cph,'Initial c');

end

end


function [T0,P0] = Read_First_PT(Cpt)

header_row = [];

for i = 1:size(Cpt,1)

    row_txt = cell(1,size(Cpt,2));

    for j = 1:size(Cpt,2)
        row_txt{j} = strtrim(Cell_String(Cpt{i,j}));
    end

    if any(strcmpi(row_txt,'P (Pa)')) && ...
       any(strcmpi(row_txt,'T (K)'))  && ...
       any(strcmpi(row_txt,'t (s)'))
        header_row = i;
        break
    end
end

if isempty(header_row)
    error('Read_First_PT: PTt_Path header row not found.')
end

colP = [];
colT = [];

for j = 1:size(Cpt,2)

    txt = strtrim(Cell_String(Cpt{header_row,j}));

    if strcmpi(txt,'P (Pa)')
        colP = j;
    elseif strcmpi(txt,'T (K)')
        colT = j;
    end
end

for i = header_row+1:size(Cpt,1)

    if Is_Empty_Cell(Cpt{i,colP}) || Is_Empty_Cell(Cpt{i,colT})
        continue
    end

    P0 = Cell_Number(Cpt{i,colP});
    T0 = Cell_Number(Cpt{i,colT});
    return
end

error('Read_First_PT: no valid P-T point found.')

end


function c0 = Get_Phase_Init_C(META,phase_name,default)

c0 = default;

for ip = 1:numel(META.phase)

    if strcmp(META.phase(ip).name,phase_name)

        if isfield(META.phase(ip),'init_c') && ~isempty(META.phase(ip).init_c)
            c0 = META.phase(ip).init_c;
        end

        return
    end
end

end


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


function M_diag = Read_Phase_Mobility(Cph,elem_list)

mob_row = [];

for i = 1:size(Cph,1)

    if strcmpi(strtrim(Cell_String(Cph{i,1})),'Mobility')
        mob_row = i;
        break
    end
end

if isempty(mob_row)
    error('Read_Phase_Mobility: Mobility table is missing.')
end

Ne     = numel(elem_list);
M_diag = zeros(1,Ne);

for ie = 1:Ne

    elem = elem_list{ie};

    col = [];
    row = [];

    for j = 2:size(Cph,2)
        if strcmpi(strtrim(Cell_String(Cph{mob_row,j})),elem)
            col = j;
            break
        end
    end

    for i = mob_row+1:size(Cph,1)

        if Is_Empty_Cell(Cph{i,1})
            continue
        end

        if strcmpi(strtrim(Cell_String(Cph{i,1})),'Facet')
            break
        end

        if strcmpi(strtrim(Cell_String(Cph{i,1})),elem)
            row = i;
            break
        end
    end

    if isempty(row) || isempty(col)
        error('Read_Phase_Mobility: missing diagonal mobility for element %s.',elem)
    end

    M_diag(ie) = Cell_Number(Cph{row,col});

end

end


function val = Get_Value(C,key,default)

val = default;
key = strtrim(key);

for i = 1:size(C,1)

    key_i = strtrim(Cell_String(C{i,1}));

    if strcmpi(key_i,key)

        if size(C,2) < 2 || Is_Empty_Cell(C{i,2})
            return
        end

        raw = C{i,2};

        if isnumeric(raw)
            val = raw;
        elseif ischar(raw) || isstring(raw)

            tmp = str2double(raw);

            if isnan(tmp)
                val = char(raw);
            else
                val = tmp;
            end

        else
            val = raw;
        end

        return
    end
end

end


function list = Get_Row_List(C,key)

list = {};
key  = strtrim(key);

for i = 1:size(C,1)

    key_i = strtrim(Cell_String(C{i,1}));

    if strcmpi(key_i,key)

        for j = 2:size(C,2)

            if ~Is_Empty_Cell(C{i,j})
                list{end+1} = strtrim(Cell_String(C{i,j})); %#ok<AGROW>
            end
        end

        return
    end
end

end


function list = Get_Row_List_Number(C,key)

list = [];
key  = strtrim(key);

for i = 1:size(C,1)

    key_i = strtrim(Cell_String(C{i,1}));

    if strcmpi(key_i,key)

        for j = 2:size(C,2)

            if ~Is_Empty_Cell(C{i,j})
                list(end+1) = Cell_Number(C{i,j}); %#ok<AGROW>
            end
        end

        return
    end
end

end


function s = Cell_String(x)

if Is_Empty_Cell(x)
    s = '';
elseif ischar(x)
    s = x;
elseif isstring(x)
    s = char(x);
elseif isnumeric(x)
    s = num2str(x);
else
    try
        s = char(string(x));
    catch
        s = '';
    end
end

end


function x = Cell_Number(c)

if isnumeric(c)
    x = c;
elseif ischar(c) || isstring(c)
    x = str2double(c);
else
    x = NaN;
end

if isnan(x)
    error('Expected numeric value.')
end

end


function tf = Is_Empty_Cell(x)

tf = false;

if isempty(x)
    tf = true;
    return
end

try
    if ismissing(x)
        tf = true;
        return
    end
catch
end

if ischar(x) && isempty(strtrim(x))
    tf = true;
    return
end

if isstring(x) && strlength(strtrim(x)) == 0
    tf = true;
    return
end

if isnumeric(x) && isscalar(x) && isnan(x)
    tf = true;
    return
end

end