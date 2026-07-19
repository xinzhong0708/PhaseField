%% Make_Map_2D_Bands_PairwiseLE.m
% Clean pseudo-1D symmetric band map with pairwise LE initialization.
%
% Default layout:
%   Kyanite | Orthopyroxene | Quartz | Garnet | Quartz | Orthopyroxene | Kyanite
%
% The global all-phase reference LE is deliberately avoided.  Instead, each
% phase is initialized from pairwise LE problems:
%
%   A + eps*B,  A + eps*C,  A + eps*D, ...
%
% where A is the major/present phase and the other phase is a small tunable
% amount.  The median of the returned A-compositions is used as the initial
% c for phase A.  This avoids sending many phases to LE_Calculator at once.

clear; clf

addpath('..\bin')
addpath('..\')
addpath('..\Thermo')
addpath('..\Thermo\Solutions')

%% ========================================================================
%  User controls
% ========================================================================

metadata_file = '..\Metadata.xlsx';
sync_PT_from_metadata  = 1;
sync_eta_from_metadata = 1;

T = 950 + 273.15;       % K
P = 0.4e9;              % Pa

Cname  = {'Fe' 'Mg' 'Ca' 'Al' 'Si' 'O'};
solmod = 'solution_models_PFM';

% Phase order used by thermodynamic arrays.
% To switch back to Cpx, change Orthopyroxene to Clinopyroxene and set c_value{1}.
phs_name = {'Orthopyroxene','Garnet','Kyanite','Quartz'};
% phs_name = {'Clinopyroxene','Garnet','Kyanite','Quartz'};

% Requested phase proportions in the same order as phs_name.
% Values are normalized below.  With nx = 320 these give exact symmetric bands:
%   Kya 50 | Opx 75 | Qtz 25 | Grt 20 | Qtz 25 | Opx 75 | Kya 50
phase_prop = [0.40 0.04 0.10 0.20];
phase_prop = phase_prop/sum(phase_prop);

% Initial independent endmember compositions.
% [] means automatic middle guess based on the number of independent variables.
c_value = cell(1,numel(phs_name));
c_value{1} = [0.30 0.10 0.36 0.2]; 
c_value{2} = [0.45 0.45];         % Garnet
c_value{3} = 1.0;                 % Kyanite
c_value{4} = 1.0;                 % Quartz

% Domain.  nx = 320 is chosen so the default normalized proportions are exact.
Lx = 100e-6;
Ly = 100e-6;
nx = 320;
ny = 4;

% Symmetric band layout by phase name.
band_name = {'Kyanite','Orthopyroxene','Quartz','Garnet', ...
             'Quartz','Orthopyroxene','Kyanite'};
% band_name = {'Kyanite','Clinopyroxene','Quartz','Garnet', ...
%              'Quartz','Clinopyroxene','Kyanite'};

% Pairwise LE reference controls.
pair_p_minor = 2e-2;       % small amount of phase B in A+B reference
pair_alpha   = 0.1;
pair_iter    = 2000;

% Penalty scaling.
PHYS        = struct();
PHYS.E_sc   = 1e9;
PHYS.L_sc   = 1;
PHYS.vref   = 2e-5;
    
E_sc = PHYS.E_sc;
L_sc = PHYS.L_sc;
vref = PHYS.vref;
eta0 = 5000e10/E_sc;

% Initial eta.  The map is sharp one-hot by default, so uniform eta is safest.
init_eta_mode = 'uniform';

% Optional initial LE check.  Keep off for the map generator; Run_2D will do LE.
do_initial_LE_check = 0;

%% ========================================================================
%  Metadata P-T / eta
% ========================================================================

if sync_PT_from_metadata == 1
    [T,P] = Read_First_PT_From_Metadata(metadata_file,T,P);
end

if sync_eta_from_metadata == 1
    eta0 = Read_Eta_From_Metadata(metadata_file,eta0,E_sc);
end

%% ========================================================================
%  Thermodynamics and pairwise reference compositions
% ========================================================================

pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref);
Nphase     = numel(pars_phase);

c_guess = Build_c_Guess(pars_phase,c_value);
e_guess = Calc_e(pars_phase,c_guess);
Ne      = numel(Calc_E_Tot(e_guess,reshape(phase_prop,1,1,Nphase)));

[c_ref,mu_ref,pair_diag] = Build_Pairwise_LE_Reference( ...
    pars_phase,c_guess,phase_prop,eta0,pair_p_minor,[pair_alpha pair_iter]);

e_ref = Calc_e(pars_phase,c_ref);

%% ========================================================================
%  Grid and symmetric band map
% ========================================================================

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

band_phase = Phase_Names_To_Index(band_name,phs_name);
band_width = Symmetric_Band_Widths(phase_prop,band_phase,nx);

Ngrain = numel(band_phase);
Np     = Ngrain;

phi       = zeros(ny,nx,Ngrain);
phase_ID  = zeros(ny,nx);
grain_ID  = zeros(ny,nx);
seed_xy   = zeros(Ngrain,2);

grain_phase = band_phase(:);
ix1 = 1;

for ig = 1:Ngrain

    ix2 = ix1 + band_width(ig) - 1;

    phi(:,ix1:ix2,ig)      = 1;
    phase_ID(:,ix1:ix2)    = grain_phase(ig);
    grain_ID(:,ix1:ix2)    = ig;
    seed_xy(ig,1)          = mean(x(ix1:ix2));
    seed_xy(ig,2)          = mean(y);

    ix1 = ix2 + 1;

end

Check_Initial_Phi_Map(phi,grain_ID)
phase_prop_geom = Compute_Phase_Fraction_From_Map(phase_ID,Nphase);

theta_grain = zeros(1,Ngrain);

%% ========================================================================
%  MODEL
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

eps_phi = 1e-14;
MODEL.p_fun = @(a,phi) phi(:,:,a).^2 ./ (sum(phi.^2,3) + eps_phi);
MODEL.dpdphi = @(a,b,phi) ...
    (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) ...
    - 2*phi(:,:,a).*phi(:,:,b).^2 ./ (sum(phi.^2,3) + eps_phi).^2;

pars = MODEL.pars;
F    = MODEL; %#ok<NASGU>

%% ========================================================================
%  Initial STATE
% ========================================================================

p = Calc_p(MODEL,phi);

eta_bulk = eta0*ones(ny,nx);

switch lower(init_eta_mode)
    case 'uniform'
        eta_init = eta_bulk;
    otherwise
        error('Only uniform init_eta_mode is used in this clean band initializer.')
end

c_phase = cell(1,Nphase);
for ip = 1:Nphase
    c_phase{ip} = cell(size(c_ref{ip}));
    for ic = 1:numel(c_ref{ip})
        c_phase{ip}{ic} = c_ref{ip}{ic}*ones(ny,nx);
    end
end

c = Expand_c_By_Phase(c_phase,MODEL.phase_index);
e = Calc_e(pars,c);
E = Calc_E_Tot(e,p);

mu_e = cell(1,Ne);
for ie = 1:Ne
    mu_e{ie} = mu_ref{ie}*ones(ny,nx);
    E{ie}    = E{ie} + mu_e{ie}./eta_init;
end

chi = repmat({zeros(ny,nx)},Ne,Ne);
omg = CalcOmegaLocal(pars,c,e,mu_e,ny,nx,Ne,Ngrain);

PARAM            = struct();
PARAM.Np         = Np;
PARAM.Ne         = Ne;
PARAM.eta        = eta_bulk;
PARAM.eta_bulk   = eta_bulk;
PARAM.eta_init   = eta_init;
PARAM.use_WScale = 0;
PARAM.T          = T;
PARAM.P          = P;
PARAM.theta_grain = theta_grain;
PARAM.init_eta_mode = init_eta_mode;
PARAM.pairwise_LE_init = 1;
PARAM.pair_p_minor = pair_p_minor;
PARAM.pair_alpha   = pair_alpha;
PARAM.pair_iter    = pair_iter;

STATE          = struct();
STATE.c        = c;
STATE.e        = e;
STATE.E        = E;
STATE.mu_e     = mu_e;
STATE.chi      = chi;
STATE.omg      = omg;
STATE.omg_raw  = omg;
STATE.phi      = phi;
STATE.p        = p;
STATE.mask     = ones(ny,nx,Np);
STATE.LE_state = [];

%% ========================================================================
%  Optional initial LE check
% ========================================================================

STATE_INI = STATE;

if do_initial_LE_check == 1

    PARAM_LE = PARAM;
    PARAM_LE.LE_mode = 'LE';
    PARAM_LE.eta = eta_init;

    if exist('LE_Run_Mode','file') == 2
        STATE = LE_Run_Mode(STATE,PARAM_LE,MODEL);
    else
        STATE = LE_Run_Mode_New(STATE,PARAM_LE,MODEL);
    end

end

p_phase = Collapse_p_By_Phase(STATE.p,MODEL.phase_index,Nphase);
phase_prop_ref = zeros(1,Nphase);

for ip = 1:Nphase
    tmp = p_phase(:,:,ip);
    phase_prop_ref(ip) = mean(tmp(:));
end

phase_prop_ref = phase_prop_ref/sum(phase_prop_ref);
phase_prop_map = phase_prop_ref;
phase_prop_le_ref = phase_prop;

%% ========================================================================
%  Plot and print diagnostics
% ========================================================================

figure(3); clf
imagesc(x*1e6,y*1e6,phase_ID)
set(gca,'YDir','normal')
axis image tight
colorbar
title('Symmetric pseudo-1D phase bands')
xlabel('x (\mum)')
ylabel('y (\mum)')
drawnow

fprintf('\nSymmetric band map saved as Map2d.mat\n')
fprintf('Band order:\n')
for ig = 1:Ngrain
    fprintf('  grain %d: %-16s width = %d pixels\n', ...
        ig,phs_name{grain_phase(ig)},band_width(ig))
end

fprintf('\nPhase             requested        realised area fraction      difference\n')
for ip = 1:Nphase
    fprintf('%-16s    %.8f         %.8f                 %+.3e\n', ...
        phs_name{ip},phase_prop(ip),phase_prop_geom(ip), ...
        phase_prop_geom(ip)-phase_prop(ip))
end

fprintf('\nPairwise LE diagnostics: failed pairs = %d of %d\n', ...
    nnz(~pair_diag.ok(:)),numel(pair_diag.ok))

if do_initial_LE_check == 1
    dc_max = Max_C_Diff(STATE.c,STATE_INI.c);
    fprintf('max|c after optional LE check - c before| = %.8e\n',dc_max)
else
    fprintf('Initial all-phase LE check skipped. Run_2D will perform the first LE.\n')
end

%% ========================================================================
%  Save
% ========================================================================

save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'E_sc','L_sc','vref','eta0','eta_bulk','eta_init', ...
    'pars','Np','Ne','Nphase','Ngrain', ...
    'T','P','Cname','solmod', ...
    'phs_name','phase_prop','phase_prop_geom','phase_prop_ref', ...
    'phase_prop_le_ref','phase_prop_map', ...
    'band_name','band_phase','band_width','grain_ID','phase_ID', ...
    'grain_phase','seed_xy','theta_grain', ...
    'pair_p_minor','pair_alpha','pair_iter','pair_diag', ...
    'c_ref','mu_ref','c_guess', ...
    'c','e','E','p','mu_e','chi','phi')

fprintf('\nSaved Map2d.mat\n')

%% ========================================================================
%  Local helper functions
% ========================================================================

function [T,P] = Read_First_PT_From_Metadata(metadata_file,T_default,P_default)

T = T_default;
P = P_default;

if exist(metadata_file,'file') ~= 2
    if exist('Metadata.xlsx','file') == 2
        metadata_file = 'Metadata.xlsx';
    elseif exist('MetaData.xlsx','file') == 2
        metadata_file = 'MetaData.xlsx';
    end
end

if exist(metadata_file,'file') ~= 2
    warning('Metadata file not found. Using map P-T controls.')
    return
end

try
    Cpt = readcell(metadata_file,'Sheet','PTt_Path');
catch
    warning('Could not read PTt_Path. Using map P-T controls.')
    return
end

header_row = [];
for i = 1:size(Cpt,1)
    row_txt = cell(1,size(Cpt,2));
    for j = 1:size(Cpt,2)
        row_txt{j} = strtrim(Cell_String_Local(Cpt{i,j}));
    end
    if any(strcmpi(row_txt,'P (Pa)')) && any(strcmpi(row_txt,'T (K)')) && any(strcmpi(row_txt,'t (s)'))
        header_row = i;
        break
    end
end

if isempty(header_row)
    warning('PTt_Path header not found. Using map P-T controls.')
    return
end

colP = find(strcmpi(strtrim(string(Cpt(header_row,:))),'P (Pa)'),1);
colT = find(strcmpi(strtrim(string(Cpt(header_row,:))),'T (K)'),1);
colt = find(strcmpi(strtrim(string(Cpt(header_row,:))),'t (s)'),1);

pt_rows = [];
for i = header_row+1:size(Cpt,1)
    if Is_Empty_Cell_Local(Cpt{i,colP}) || Is_Empty_Cell_Local(Cpt{i,colT}) || Is_Empty_Cell_Local(Cpt{i,colt})
        continue
    end
    pt_rows(end+1,:) = [Cell_Number_Local(Cpt{i,colt}), ...
                        Cell_Number_Local(Cpt{i,colT}), ...
                        Cell_Number_Local(Cpt{i,colP})]; %#ok<AGROW>
end

if isempty(pt_rows)
    warning('PTt_Path contains no valid rows. Using map P-T controls.')
    return
end

[~,ord] = sort(pt_rows(:,1));
pt_rows = pt_rows(ord,:);
T = pt_rows(1,2);
P = pt_rows(1,3);

fprintf('Map P-T from metadata: T = %.12g K, P = %.12g Pa\n',T,P)

end

function eta0 = Read_Eta_From_Metadata(metadata_file,eta0,E_sc)

if exist(metadata_file,'file') ~= 2
    if exist('Metadata.xlsx','file') == 2
        metadata_file = 'Metadata.xlsx';
    elseif exist('MetaData.xlsx','file') == 2
        metadata_file = 'MetaData.xlsx';
    end
end

if exist(metadata_file,'file') ~= 2
    warning('Metadata file not found. Using map eta control.')
    return
end

try
    Cmain = readcell(metadata_file,'Sheet','Main');
catch
    warning('Could not read metadata Main sheet. Using map eta control.')
    return
end

eta_SI = Get_Metadata_Value_Any(Cmain,{'Penalty eta','eta'},[]);
if ~isempty(eta_SI) && isnumeric(eta_SI) && isfinite(eta_SI)
    eta0 = eta_SI/E_sc;
end

end

function val = Get_Metadata_Value_Any(C,keys,default)

val = default;

for ik = 1:numel(keys)
    key = strtrim(keys{ik});
    for i = 1:size(C,1)
        if strcmpi(strtrim(Cell_String_Local(C{i,1})),key)
            if size(C,2) < 2 || Is_Empty_Cell_Local(C{i,2})
                return
            end
            raw = C{i,2};
            if isnumeric(raw)
                val = raw;
            else
                tmp = str2double(string(raw));
                if isnan(tmp)
                    val = raw;
                else
                    val = tmp;
                end
            end
            return
        end
    end
end

end

function pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref)

Nphase = numel(phs_name);
pars_phase = cell(1,Nphase);

for ip = 1:Nphase

    td = init_thermo(phs_name(ip),Cname,solmod);
    g0 = cell2mat(tl_g0(T,P,td));
    n  = td.n_em(:,1:end-1);

    pars            = td;
    pars.n          = n;
    pars.P          = P;
    pars.T          = T;
    pars.g0         = g0;
    pars.E_sc       = E_sc;
    pars.vref       = vref;
    pars.phase_name = phs_name(ip);

    pars_phase{ip} = pars;

end

end

function c_guess = Build_c_Guess(pars_phase,c_value)

Nphase = numel(pars_phase);
c_guess = cell(1,Nphase);

for ip = 1:Nphase

    if ip <= numel(c_value) && ~isempty(c_value{ip})
        v = c_value{ip};
    else
        Nc = size(pars_phase{ip}.n,2);
        if Nc <= 0
            v = 1.0;
        else
            v = 0.5*ones(1,Nc);
        end
    end

    c_guess{ip} = num2cell(v);

end

end

function [c_ref,mu_ref,diag] = Build_Pairwise_LE_Reference( ...
    pars_phase,c_guess,phase_prop,eta_ref,pminor,level)

Nphase = numel(pars_phase);
e_guess = Calc_e(pars_phase,c_guess);
E_tmp   = Calc_E_Tot(e_guess,reshape(phase_prop,1,1,Nphase));
Ne      = numel(E_tmp);

store_c = cell(1,Nphase);
for ip = 1:Nphase
    store_c{ip} = cell(1,numel(c_guess{ip}));
    for ic = 1:numel(c_guess{ip})
        store_c{ip}{ic} = [];
    end
end

store_mu = cell(1,Ne);
for ie = 1:Ne
    store_mu{ie} = [];
end

diag.ok = false(Nphase,Nphase);
diag.msg = strings(Nphase,Nphase);

for ia = 1:Nphase

    for ib = 1:Nphase

        if ia == ib
            continue
        end

        pars_pair = pars_phase([ia ib]);
        c_pair    = {c_guess{ia},c_guess{ib}};
        p_pair    = zeros(1,1,2);
        p_pair(1,1,1) = 1 - pminor;
        p_pair(1,1,2) = pminor;

        e_pair = Calc_e(pars_pair,c_pair);
        E_pair = Calc_E_Tot(e_pair,p_pair);

        try
            [c_eq,mu_eq] = LE_Calculator(pars_pair,p_pair,c_pair,E_pair,eta_ref,level);
            diag.ok(ia,ib) = true;
            diag.msg(ia,ib) = "ok";
        catch ME
            c_eq = c_pair;
            mu_eq = [];
            diag.ok(ia,ib) = false;
            diag.msg(ia,ib) = string(ME.message);
            warning('Pairwise LE failed for phase %d + phase %d: %s',ia,ib,ME.message)
        end

        for ic = 1:numel(c_guess{ia})
            store_c{ia}{ic}(end+1) = c_eq{1}{ic}; %#ok<AGROW>
        end

        if ~isempty(mu_eq)
            for ie = 1:Ne
                store_mu{ie}(end+1) = mu_eq{ie}; %#ok<AGROW>
            end
        end

    end
end

c_ref = c_guess;

for ip = 1:Nphase
    for ic = 1:numel(c_guess{ip})
        vals = store_c{ip}{ic};
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            c_ref{ip}{ic} = median(vals);
        end
    end
end

mu_ref = cell(1,Ne);
for ie = 1:Ne
    vals = store_mu{ie};
    vals = vals(isfinite(vals));
    if isempty(vals)
        mu_ref{ie} = 0;
    else
        mu_ref{ie} = median(vals);
    end
end

end

function idx = Phase_Names_To_Index(names,phs_name)

idx = zeros(1,numel(names));

for i = 1:numel(names)
    id = find(strcmpi(phs_name,names{i}),1);
    if isempty(id)
        error('Band phase %s is not listed in phs_name.',names{i})
    end
    idx(i) = id;
end

end

function band_width = Symmetric_Band_Widths(phase_prop,band_phase,nx)

% This assumes the layout [A B C D C B A].
if numel(band_phase) ~= 7 || any(band_phase ~= fliplr(band_phase))
    error('Symmetric_Band_Widths expects a 7-band symmetric layout.')
end

A = band_phase(1);
B = band_phase(2);
C = band_phase(3);
D = band_phase(4);

wA = max(1,round(nx*phase_prop(A)/2));
wB = max(1,round(nx*phase_prop(B)/2));
wC = max(1,round(nx*phase_prop(C)/2));
wD = nx - 2*(wA+wB+wC);

if wD < 1
    error('Central band has non-positive width. Increase nx or reduce side phase proportions.')
end

band_width = [wA wB wC wD wC wB wA];

end

function omg = CalcOmegaLocal(pars,c,e,mu_e,ny,nx,Ne,Np)

omg = zeros(ny,nx,Np);

for ip = 1:Np
    A           = PhaseThermo(pars{ip},c{ip});
    omg(:,:,ip) = reshape(A.g,ny,[]);
    for ie = 1:Ne
        omg(:,:,ip) = omg(:,:,ip) - e{ip}{ie}.*mu_e{ie};
    end
end

end

function c = Expand_c_By_Phase(c_phase,phase_index)

Ngrain = numel(phase_index);
c = cell(1,Ngrain);

for ig = 1:Ngrain
    c{ig} = c_phase{phase_index(ig)};
end

end

function p_phase = Collapse_p_By_Phase(p_grain,phase_index,Nphase)

[ny,nx,~] = size(p_grain);
p_phase = zeros(ny,nx,Nphase);

for ip = 1:Nphase
    grains = find(phase_index == ip);
    if ~isempty(grains)
        p_phase(:,:,ip) = sum(p_grain(:,:,grains),3);
    end
end

end

function phase_prop_geom = Compute_Phase_Fraction_From_Map(phase_ID,Nphase)

phase_prop_geom = zeros(1,Nphase);
N = numel(phase_ID);

for ip = 1:Nphase
    phase_prop_geom(ip) = nnz(phase_ID == ip)/N;
end

end

function Check_Initial_Phi_Map(phi,grain_ID)

s = sum(phi,3);
if any(abs(s(:)-1) > 0)
    error('Check_Initial_Phi_Map: initial phi is not one-hot.')
end

[~,gid] = max(phi,[],3);
if any(gid(:) ~= grain_ID(:))
    error('Check_Initial_Phi_Map: phi and grain_ID are inconsistent.')
end

end

function d = Max_C_Diff(c1,c0)

d = 0;

for ig = 1:numel(c1)
    for ic = 1:numel(c1{ig})
        d = max(d,max(abs(c1{ig}{ic}(:)-c0{ig}{ic}(:))));
    end
end

end

function val = Cell_Number_Local(x)

if isnumeric(x)
    val = x;
else
    val = str2double(string(x));
end

if ~isfinite(val)
    error('Metadata contains a nonnumeric P-T-time value.')
end

end

function s = Cell_String_Local(x)

if ismissing(x)
    s = '';
elseif ischar(x)
    s = x;
elseif isstring(x)
    s = char(x);
elseif isnumeric(x)
    s = num2str(x);
else
    s = '';
end

end

function tf = Is_Empty_Cell_Local(x)

if isempty(x)
    tf = true;
elseif ismissing(x)
    tf = true;
elseif ischar(x) || isstring(x)
    tf = strlength(string(x)) == 0;
else
    tf = false;
end

end
