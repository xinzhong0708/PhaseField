clear; clf

addpath('..\bin')
addpath('..\')
addpath('..\Thermo')
addpath('..\Thermo\Solutions')

%% ========================================================================
%  User controls
% ========================================================================

metadata_file          = '..\Metadata.xlsx';
sync_PT_from_metadata  = 1;
sync_eta_from_metadata = 1;

if sync_PT_from_metadata == 1
    [T,P] = Read_First_PT_From_Metadata(metadata_file,1,1);
end

Cname  = {'Fe' 'Mg' 'Ca' 'Al' 'Si' 'O'};
solmod = 'solution_models_PFM';

% Thermodynamic phase order used everywhere in MODEL.phase_index.
phs_name = {'Garnet','Clinopyroxene','Kyanite','Quartz','Orthopyroxene'};

% Requested phase proportions in the phs_name order above.
% They do not need to sum to one.
phase_prop = [0.05 0.30  0.15 0.22 0.2];
phase_prop = phase_prop/sum(phase_prop);

% Initial independent endmember compositions in the phs_name order above.
c_value = cell(1,numel(phs_name));
c_value{1} = [0.50  0.45]; 
c_value{2} = [0.25  0.35  0.25  0.20];         
c_value{3} = [1];         
c_value{4} = [1];         
c_value{5} = [0.25  0.35  0.25  0.20];         

% Symmetric band order.  The middle band is the Garnet seed.
band_phase_name = {'Quartz','Orthopyroxene','Clinopyroxene','Kyanite','Garnet','Kyanite','Clinopyroxene','Orthopyroxene','Quartz'};

% Domain.  nx=320 gives exact integer widths for the default phase_prop:
% kya/cpx/qtz/grt/qtz/cpx/kya = 50/75/25/20/25/75/50.
Lx = 500e-6;
Ly = 500e-6;
nx = 300;
ny = 4;

% Scaling / penalty
PHYS        = struct();
PHYS.E_sc   = 1e9;
PHYS.L_sc   = 1;
PHYS.vref   = 2e-5;

E_sc        = PHYS.E_sc;
L_sc        = PHYS.L_sc;
vref        = PHYS.vref;
eta0        = 5000e10/E_sc;

% Initial eta used to build E.
init_eta_mode        = 'interface_damped';  % 'uniform' or 'interface_damped'
init_eta_damp_factor = 0.1;
init_eta_q2          = 4;
init_eta_p02         = 3e-3;
init_eta_q3          = 4;
init_eta_p03         = 3e-3;
init_eta_nsmooth     = 1;
init_eta_halo_cut    = 1e-3;

% Reference LE controls
init_ref_phase_prop_mode = 'requested';     % 'requested' or 'realised_p'
init_ref_eta_mode        = 'bulk';          % 'bulk' or 'harmonic_init'

% Initial nonlinear LE projection before saving
do_initial_LE = 1;

if sync_eta_from_metadata == 1
    [eta0,init_eta_damp_factor] = Read_Init_Eta_From_Metadata( ...
        metadata_file,eta0,init_eta_damp_factor,E_sc);
end

%% ========================================================================
%  Thermodynamics
% ========================================================================

pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref);
Nphase     = numel(pars_phase);

c_guess = cell(1,Nphase);
for ip = 1:Nphase
    c_guess{ip} = num2cell(c_value{ip});
end

e_guess = Calc_e(pars_phase,c_guess);
E_guess = Calc_E_Tot(e_guess,reshape(phase_prop,1,1,Nphase));
Ne      = numel(E_guess);

%% ========================================================================
%  Grid
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

%% ========================================================================
%  Symmetric phase/grain map
% ========================================================================

grain_phase = Phase_Names_To_ID(phs_name,band_phase_name);
Ngrain      = numel(grain_phase);
Np          = Ngrain;

[band_width,band_fraction] = Local_Symmetric_Band_Widths( ...
    phase_prop,grain_phase,nx);

phi      = zeros(ny,nx,Ngrain);
phase_ID = zeros(ny,nx);
grain_ID = zeros(ny,nx);
seed_xy  = zeros(Ngrain,2);

ix1 = 1;
for ig = 1:Ngrain

    ix2 = ix1 + band_width(ig) - 1;
    iph = grain_phase(ig);

    phi(:,ix1:ix2,ig)   = 1;
    phase_ID(:,ix1:ix2) = iph;
    grain_ID(:,ix1:ix2) = ig;

    seed_xy(ig,1) = mean(x(ix1:ix2));
    seed_xy(ig,2) = mean(y);

    ix1 = ix2 + 1;

end

Check_Initial_Phi_Map(phi,grain_ID)

phase_prop_geom = Compute_Phase_Fraction_From_Map(phase_ID,Nphase);

% Zero orientation for pseudo-1D symmetric bands.
theta_grain = zeros(1,Ngrain);

fprintf('Band order:\n')
for ig = 1:Ngrain
    fprintf('  grain %d: %-16s width = %d cells, x = %.4e m\n', ...
        ig,phs_name{grain_phase(ig)},band_width(ig),seed_xy(ig,1)*L_sc)
end

fprintf('\nPhase             requested        realised area fraction      difference\n')
for ip = 1:Nphase
    fprintf('%-16s    %.8f         %.8f                 %+.3e\n', ...
        phs_name{ip},phase_prop(ip),phase_prop_geom(ip), ...
        phase_prop_geom(ip)-phase_prop(ip))
end

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

MODEL.p_fun = @(a,phi) ...
    phi(:,:,a).^2 ./ (sum(phi.^2,3) + eps_phi);

MODEL.dpdphi = @(a,b,phi) ...
    (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) ...
    - 2*phi(:,:,a).*phi(:,:,b).^2 ./ (sum(phi.^2,3) + eps_phi).^2;

pars = MODEL.pars;
F    = MODEL; %#ok<NASGU>

%% ========================================================================
%  Plot map
% ========================================================================

figure(3); clf
imagesc(x*1e6,y*1e6,phase_ID)
set(gca,'YDir','normal')
axis image tight
colorbar
title('Symmetric pseudo-1D map: kya cpx qtz grt qtz cpx kya')
xlabel('x (\mum)')
ylabel('y (\mum)')
drawnow

%% ========================================================================
%  Penalty-consistent initial E, c, mu
% ========================================================================

p = Calc_p(MODEL,phi);

p_phase_ref    = Collapse_p_By_Phase(p,MODEL.phase_index,Nphase);
phase_prop_ref = zeros(1,Nphase);

for ip = 1:Nphase
    tmp = p_phase_ref(:,:,ip);
    phase_prop_ref(ip) = mean(tmp(:));
end

phase_prop_ref = phase_prop_ref/sum(phase_prop_ref);

eta_bulk = eta0*ones(ny,nx);

switch lower(init_eta_mode)

    case 'uniform'

        eta_init = eta_bulk;

    case 'interface_damped'

        if exist('Eta_Damping_SmoothHalo','file') == 2
            eta_init = Eta_Damping_SmoothHalo( ...
                p_phase_ref,eta_bulk,init_eta_damp_factor*eta_bulk, ...
                init_eta_q2,init_eta_p02,init_eta_q3,init_eta_p03, ...
                init_eta_damp_factor*eta_bulk,init_eta_nsmooth, ...
                init_eta_damp_factor*eta_bulk,init_eta_halo_cut);
        else
            warning('Eta_Damping_SmoothHalo not found. Using uniform eta.')
            eta_init = eta_bulk;
        end

    otherwise

        error('Unknown init_eta_mode: %s',init_eta_mode)

end

switch lower(init_ref_phase_prop_mode)
    case 'requested'
        phase_prop_le_ref = phase_prop;
    case 'realised_p'
        phase_prop_le_ref = phase_prop_ref;
    otherwise
        error('Unknown init_ref_phase_prop_mode: %s',init_ref_phase_prop_mode)
end

phase_prop_le_ref = phase_prop_le_ref/sum(phase_prop_le_ref);

switch lower(init_ref_eta_mode)
    case 'bulk'
        eta_ref = eta0;
    case 'harmonic_init'
        eta_ref = 1/mean(1./eta_init(:));
    otherwise
        error('Unknown init_ref_eta_mode: %s',init_ref_eta_mode)
end

p_ref    = reshape(phase_prop_le_ref,1,1,Nphase);
E_target = Calc_E_Tot(e_guess,p_ref);

[c_ref,mu_ref] = LE_Calculator(pars_phase,p_ref,c_guess,E_target,eta_ref,[0.1,2000]);

e_ref     = Calc_e(pars_phase,c_ref);
E_mix_ref = Calc_E_Tot(e_ref,p_ref);

E_offset_ref = cell(1,Ne);
E_offset     = cell(1,Ne);
omega_ref    = zeros(1,Nphase);

for ie = 1:Ne
    E_offset_ref{ie} = E_target{ie} - E_mix_ref{ie};
    E_offset{ie}     = mu_ref{ie}./eta_init;
end

for ip = 1:Nphase

    omega_ref(ip) = PhaseG(pars_phase{ip},c_ref{ip});

    for ie = 1:Ne
        omega_ref(ip) = omega_ref(ip) - e_ref{ip}{ie}.*mu_ref{ie};
    end

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

for ie = 1:Ne
    E{ie} = E{ie} + E_offset{ie};
end

E_bulk_shift = zeros(1,Ne);

mu_e = cell(1,Ne);
for ie = 1:Ne
    mu_e{ie} = mu_ref{ie}*ones(ny,nx);
end

chi = repmat({zeros(ny,nx)},Ne,Ne);

fprintf('\nInitial eta mode = %s, eta range = %.4e to %.4e, eta_ref = %.4e\n', ...
    init_eta_mode,min(eta_init(:)),max(eta_init(:)),eta_ref)
fprintf('Reference LE phase proportions (%s):',init_ref_phase_prop_mode)
fprintf(' %.8f',phase_prop_le_ref)
fprintf('\n')

%% ========================================================================
%  PARAM and STATE
% ========================================================================

PARAM            = struct();
PARAM.Np         = Np;
PARAM.Ne         = Ne;
PARAM.eta        = eta_bulk;
PARAM.eta_bulk   = eta_bulk;
PARAM.eta_init   = eta_init;
PARAM.use_WScale = 0;
PARAM.T          = T;
PARAM.P          = P;

PARAM.init_eta_mode            = init_eta_mode;
PARAM.init_eta_damp_factor     = init_eta_damp_factor;
PARAM.init_eta_synced_E        = 1;
PARAM.init_ref_phase_prop_mode = init_ref_phase_prop_mode;
PARAM.init_ref_eta_mode        = init_ref_eta_mode;
PARAM.theta_grain              = theta_grain;

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

Check_State_Size(STATE,MODEL,PARAM,Ngrain,Ne)

%% ========================================================================
%  Initial LE projection
% ========================================================================

STATE_INI = STATE;

if do_initial_LE == 1

    PARAM_LE         = PARAM;
    PARAM_LE.LE_mode = 'LE';
    PARAM_LE.eta     = eta_init;

    if exist('LE_Run_Mode','file') == 2
        STATE = LE_Run_Mode(STATE,PARAM_LE,MODEL);
    elseif exist('LE_Run_Mode_New','file') == 2
        STATE = LE_Run_Mode_New(STATE,PARAM_LE,MODEL);
    else
        error('Neither LE_Run_Mode nor LE_Run_Mode_New was found.')
    end

end

p_phase        = Collapse_p_By_Phase(STATE.p,MODEL.phase_index,Nphase);
phase_prop_map = zeros(1,Nphase);

for ip = 1:Nphase
    tmp = p_phase(:,:,ip);
    phase_prop_map(ip) = mean(tmp(:));
end

dc_max = 0;

for ig = 1:Ngrain
    for ic = 1:numel(STATE.c{ig})
        dc_max = max(dc_max, ...
            max(abs(STATE.c{ig}{ic}(:)-STATE_INI.c{ig}{ic}(:))));
    end
end

mu_jump = zeros(1,Ne);
E_mu_consistency = zeros(1,Ne);
E_mix_check = Calc_E_Tot(STATE.e,STATE.p);

for ie = 1:Ne
    mu_jump(ie) = max(STATE.mu_e{ie}(:)) - min(STATE.mu_e{ie}(:));
    E_mu_consistency(ie) = max(abs( ...
        eta_init(:).*(STATE.E{ie}(:)-E_mix_check{ie}(:)) - STATE.mu_e{ie}(:)));
end

fprintf('\nAfter initial LE projection:\n')
fprintf('max|c after first LE - c before LE| = %.8e\n',dc_max)
fprintf('max|eta_init*(E-Emix)-mu_e|        = %.8e\n',max(E_mu_consistency))
fprintf('initial mu_e jump after LE:\n')
disp(mu_jump)

% Keep legacy variables consistent.
c    = STATE.c;
e    = STATE.e;
E    = STATE.E;
p    = STATE.p;
mu_e = STATE.mu_e;
chi  = STATE.chi;
phi  = STATE.phi;
eta  = eta_bulk;

%% ========================================================================
%  Save
% ========================================================================

map_mode        = 'bands';
grain_size      = NaN;
grain_size_real = NaN;
Ngrain_user     = [];
rng_seed        = [];
periodic_map    = 0;

save('Map2d.mat', ...
    'PHYS','GRID','MODEL','PARAM','STATE', ...
    'E_sc','L_sc','vref','eta','eta_bulk','eta_init', ...
    'pars','Np','Ne','Nphase','Ngrain', ...
    'T','P','Cname','solmod', ...
    'phs_name','phase_prop','phase_prop_geom','phase_prop_ref', ...
    'phase_prop_le_ref','phase_prop_map', ...
    'map_mode','grain_ID','phase_ID','grain_phase','seed_xy','theta_grain', ...
    'grain_size','grain_size_real','Ngrain_user','rng_seed','periodic_map', ...
    'band_phase_name','band_width','band_fraction', ...
    'init_ref_phase_prop_mode','init_ref_eta_mode', ...
    'init_eta_mode','init_eta_damp_factor','init_eta_q2','init_eta_p02', ...
    'init_eta_q3','init_eta_p03','init_eta_nsmooth','init_eta_halo_cut', ...
    'eta_ref','E_target','E_offset','E_offset_ref','E_bulk_shift', ...
    'c_ref','mu_ref','omega_ref', ...
    'c','e','E','p','mu_e','chi','phi')

fprintf('\nSaved Map2d.mat\n')
fprintf('Saved symmetric band order:\n')
disp(band_phase_name)

disp([mean(STATE.E{1},'all') mean(STATE.E{2},'all') mean(STATE.E{3},'all')  mean(STATE.E{4},'all')])


%% ========================================================================
%  Local helper functions
% ========================================================================

function phase_id = Phase_Names_To_ID(phs_name,phase_name)

phase_id = zeros(1,numel(phase_name));

for i = 1:numel(phase_name)

    id = find(strcmpi(phs_name,phase_name{i}),1,'first');

    if isempty(id)
        error('Phase "%s" is not in phs_name.',phase_name{i})
    end

    phase_id(i) = id;

end

end


function [band_width,band_fraction] = Local_Symmetric_Band_Widths(phase_prop,band_order,nx)

Nband = numel(band_order);

if mod(Nband,2) ~= 1
    error('Local_Symmetric_Band_Widths requires an odd number of bands.')
end

ic = (Nband+1)/2;

for i = 1:ic-1
    if band_order(i) ~= band_order(Nband+1-i)
        error('Band order is not symmetric.')
    end
end

Nphase = numel(phase_prop);
Noccur = accumarray(band_order(:),1,[Nphase,1]).';

band_fraction = phase_prop(band_order)./Noccur(band_order);
band_fraction = band_fraction/sum(band_fraction);

target = nx*band_fraction;
npair  = ic-1;

best_err = inf;
best_w   = [];

for wc = 1:nx-2*npair

    if mod(nx-wc,2) ~= 0
        continue
    end

    side = (nx-wc)/2;

    if side < npair
        continue
    end

    wleft = Allocate_Integer_Widths(target(1:npair),side);
    w     = [wleft,wc,fliplr(wleft)];

    err = sum((w-target).^2);

    if err < best_err
        best_err = err;
        best_w   = w;
    end

end

if isempty(best_w)
    error('Could not build symmetric band widths.')
end

band_width = best_w;

end


function width = Allocate_Integer_Widths(target,total_width)

target = target(:).';
target = max(target,eps);
raw    = total_width*target/sum(target);

width = floor(raw);
width(raw > 0 & width < 1) = 1;

while sum(width) < total_width

    res = raw - width;
    [~,ord] = sort(res,'descend');

    for k = 1:numel(ord)
        width(ord(k)) = width(ord(k)) + 1;
        if sum(width) == total_width
            break
        end
    end

end

while sum(width) > total_width

    res = raw - width;
    [~,ord] = sort(res,'ascend');

    changed = 0;

    for k = 1:numel(ord)

        ii = ord(k);

        if width(ii) > 1
            width(ii) = width(ii) - 1;
            changed = 1;
        end

        if sum(width) == total_width
            break
        end
    end

    if changed == 0
        error('Allocate_Integer_Widths failed.')
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


function Check_State_Size(STATE,MODEL,PARAM,Ngrain,Ne)

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

end


function [T,P] = Read_First_PT_From_Metadata(metadata_file,T_default,P_default)

T = T_default;
P = P_default;

metadata_file = Resolve_Metadata_File(metadata_file);

if exist(metadata_file,'file') ~= 2
    warning('Metadata file was not found. Using map P-T controls.')
    return
end

try
    Cpt = readcell(metadata_file,'Sheet','PTt_Path');
catch
    warning('Could not read PTt_Path. Using map P-T controls.')
    return
end

header = [];
for i = 1:size(Cpt,1)

    row = strings(1,size(Cpt,2));

    for j = 1:size(Cpt,2)
        row(j) = string(Cell_String_Local(Cpt{i,j}));
    end

    if any(strcmpi(row,'P (Pa)')) && ...
       any(strcmpi(row,'T (K)'))  && ...
       any(strcmpi(row,'t (s)'))
        header = i;
        break
    end

end

if isempty(header)
    warning('PTt_Path header not found. Using map P-T controls.')
    return
end

colP = Find_Header(Cpt,header,'P (Pa)');
colT = Find_Header(Cpt,header,'T (K)');
colt = Find_Header(Cpt,header,'t (s)');

pt = [];

for i = header+1:size(Cpt,1)

    if Is_Empty_Cell_Local(Cpt{i,colP}) || ...
       Is_Empty_Cell_Local(Cpt{i,colT}) || ...
       Is_Empty_Cell_Local(Cpt{i,colt})
        continue
    end

    pt(end+1,:) = [ ...
        Cell_Number_Local(Cpt{i,colt}), ...
        Cell_Number_Local(Cpt{i,colT}), ...
        Cell_Number_Local(Cpt{i,colP})]; %#ok<AGROW>

end

if isempty(pt)
    warning('PTt_Path contains no valid rows. Using map P-T controls.')
    return
end

[~,ord] = sort(pt(:,1));
pt = pt(ord,:);

T = pt(1,2);
P = pt(1,3);

fprintf('Map P-T from %s: T = %.12g K, P = %.12g Pa\n',metadata_file,T,P)

end


function [eta0,init_eta_damp_factor] = Read_Init_Eta_From_Metadata( ...
    metadata_file,eta0,init_eta_damp_factor,E_sc)

metadata_file = Resolve_Metadata_File(metadata_file);

if exist(metadata_file,'file') ~= 2
    warning('Metadata file was not found. Using map eta controls.')
    return
end

try
    Cmain = readcell(metadata_file,'Sheet','Main');
catch
    warning('Could not read metadata Main sheet. Using map eta controls.')
    return
end

eta_SI = Get_Metadata_Value(Cmain,{'Penalty eta','eta'},[]);

if ~isempty(eta_SI) && isnumeric(eta_SI) && isfinite(eta_SI)
    eta0 = eta_SI/E_sc;
end

tmp = Get_Metadata_Value(Cmain, ...
    {'interface damping factor','Interface damping factor'}, ...
    init_eta_damp_factor);

if isnumeric(tmp) && isfinite(tmp)
    init_eta_damp_factor = min(max(tmp,eps),1);
end

end


function metadata_file = Resolve_Metadata_File(metadata_file)

if exist(metadata_file,'file') == 2
    return
end

if exist('Metadata.xlsx','file') == 2
    metadata_file = 'Metadata.xlsx';
elseif exist('MetaData.xlsx','file') == 2
    metadata_file = 'MetaData.xlsx';
end

end


function col = Find_Header(C,header,key)

col = [];

for j = 1:size(C,2)

    if strcmpi(strtrim(Cell_String_Local(C{header,j})),key)
        col = j;
        return
    end

end

error('Find_Header: key "%s" not found.',key)

end


function val = Get_Metadata_Value(C,keys,default)

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

end


function val = Cell_Number_Local(x)

if isnumeric(x)
    val = x;
elseif ischar(x) || isstring(x)
    val = str2double(x);
else
    val = NaN;
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
