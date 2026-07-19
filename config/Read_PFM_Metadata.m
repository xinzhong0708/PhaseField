function [PHYS,NUM,PARAM,MODEL,GRID] = Read_PFM_Metadata(xlsx_file,GRID,MODEL,STATE,eta,PARAM)
%READ_PFM_METADATA Read PFM metadata from Excel and scale to code units.
%
% Excel metadata is assumed to be in SI / physical units.
%
% Scaling:
%
%   x_code      = x_SI / L_sc
%   t_code      = t_SI / t_sc
%   sigma_code  = sigma_SI / (E_sc*L_sc)
%   kappa_code  = kappa_SI / (E_sc*L_sc^2)
%   eta_code    = eta_SI / E_sc
%
% Required sheets:
%   Main
%   PTt_Path
%   one sheet for each phase in MODEL.phs_name
%
% Mobility table is interpreted as constant elemental mobility in SI.
% It is scaled here and stored in PHYS.M_phs.

if nargin < 6 || isempty(PARAM)
    PARAM = struct();
end

% ------------------------------------------------------------
% Read sheets
% ------------------------------------------------------------
Cmain      = readcell(xlsx_file,'Sheet','Main');
Cpt        = readcell(xlsx_file,'Sheet','PTt_Path');
sheet_list = sheetnames(xlsx_file);

if ~isfield(MODEL,'phs_name') || isempty(MODEL.phs_name)
    error('Read_PFM_Metadata: MODEL.phs_name is missing. Active phases must come from the map.')
end

% ------------------------------------------------------------
% Main metadata
% ------------------------------------------------------------
elem_list = Get_Row_List(Cmain,'Independent c (Si dependent)');
solmod    = Get_Value_Any(Cmain,{'Solution name','Solution model'},'solution_models_PFM');
vref      = Get_Value_Any(Cmain,{'Molar volume','vref'},2e-5);

if isempty(elem_list)
    error('Read_PFM_Metadata: element list is missing in Main sheet.')
end

% ------------------------------------------------------------
% Check phase sheets
% ------------------------------------------------------------
Nphase = numel(MODEL.phs_name);

for ip = 1:Nphase

    phase_name = MODEL.phs_name{ip};

    if ~any(strcmp(sheet_list,phase_name))
        error('Read_PFM_Metadata: missing phase sheet "%s".',phase_name)
    end

end

% ------------------------------------------------------------
% PHYS: scales
% ------------------------------------------------------------
PHYS             = struct();

PHYS.E_sc        = Get_Value_Any(Cmain,{'Energy scale','Energy scaling'},1e9);
PHYS.L_sc        = Get_Value_Any(Cmain,{'Length scale','Length scaling'},1e-6);
PHYS.t_sc        = Get_Value_Any(Cmain,{'Time scale','Time scaling'},1);
PHYS.vref        = vref;
PHYS.dceq        = Get_Value_Any(Cmain,{'dceq','dc eq'},0.5);

% ------------------------------------------------------------
% Scale GRID once
% ------------------------------------------------------------
if ~isfield(GRID,'is_scaled') || GRID.is_scaled == 0

    GRID.x_SI  = GRID.x;
    GRID.y_SI  = GRID.y;
    GRID.dx_SI = GRID.dx;
    GRID.dy_SI = GRID.dy;

    if isfield(GRID,'Lx')
        GRID.Lx_SI = GRID.Lx;
    end

    if isfield(GRID,'Ly')
        GRID.Ly_SI = GRID.Ly;
    end

    GRID.x  = GRID.x /PHYS.L_sc;
    GRID.y  = GRID.y /PHYS.L_sc;
    GRID.dx = GRID.dx/PHYS.L_sc;
    GRID.dy = GRID.dy/PHYS.L_sc;

    if isfield(GRID,'Lx')
        GRID.Lx = GRID.Lx/PHYS.L_sc;
    end

    if isfield(GRID,'Ly')
        GRID.Ly = GRID.Ly/PHYS.L_sc;
    end

    GRID.is_scaled = 1;

end

% ------------------------------------------------------------
% Interface parameters from SI to code units
% ------------------------------------------------------------
thick_fac        = Get_Value_Any(Cmain,{'Interface thick factor','Interface thickness factor'},3);
PHYS.l           = thick_fac*GRID.dx;

kappa_SI         = Get_Value_Any(Cmain,{'4th order kappa','kappa'},0);
PHYS.kappa_SI    = kappa_SI;
PHYS.kappa       = kappa_SI/(PHYS.E_sc*PHYS.L_sc^2);

sigma_phase_SI   = zeros(1,Nphase);

for ip = 1:Nphase

    Cph = readcell(xlsx_file,'Sheet',MODEL.phs_name{ip});
    sigma_phase_SI(ip) = Get_Value_Any(Cph,{'Interface energy','Interface energy J/m2','sigma_J_m2'},0.5);

end

PHYS.sigma_phase_SI = sigma_phase_SI;
PHYS.sigma_phase    = sigma_phase_SI/(PHYS.E_sc*PHYS.L_sc);

% Current solver still uses one scalar sigma.
% For now use the mean phase value.
PHYS.sigma      = mean(PHYS.sigma_phase);

PHYS.m          = 6*PHYS.sigma/PHYS.l;
PHYS.kap        = 3/4*PHYS.sigma*PHYS.l;

% ------------------------------------------------------------
% Global W/kappa interface ramp geometry
% ------------------------------------------------------------
zero_grid = Get_Value_Any(Cmain, {'Interface zero kappa grid','Interface zero kappa grid'},2);
ramp_grid = Get_Value_Any(Cmain, {'Interface ramp grid','Interface kappa ramp grid'},6);
zero_grid = max(0,round(zero_grid));
ramp_grid = max(1,round(ramp_grid));
PARAM.Interface_zero_kapp_grid = zero_grid;
PARAM.Interface_ramp_grid      = ramp_grid;
PARAM.ramp_zero_width = zero_grid;
PARAM.ramp_width      = ramp_grid;

% ------------------------------------------------------------
% Eta
% ------------------------------------------------------------
eta_SI = Get_Value_Any(Cmain,{'Penalty eta','eta'},[]);

if isempty(eta_SI)

    % eta from Map2d.mat is normally already scaled.  Prefer an explicit
    % bulk eta if the map also stores an initialization-specific eta map.
    if isfield(PARAM,'eta_bulk') && ~isempty(PARAM.eta_bulk)
        PHYS.eta = PARAM.eta_bulk;
    else
        PHYS.eta = eta;
    end

else

    PHYS.eta = eta_SI/PHYS.E_sc*ones(GRID.ny,GRID.nx);

end

% ------------------------------------------------------------
% NUM: all time values from SI to code units
% ------------------------------------------------------------
NUM                     = struct();

NUM.dt_phy              = Get_Value_Any(Cmain,{'dt initial','Initial dt'},1e-2)/PHYS.t_sc;
NUM.dt_max              = Get_Value_Any(Cmain,{'dt max','Maximum dt'},100)/PHYS.t_sc;
NUM.dt_min              = Get_Value_Any(Cmain,{'dt min','Minimum dt'},1e-15)/PHYS.t_sc;

NUM.t_phy               = 0;
NUM.time                = 0;

NUM.dE_target           = Get_Value_Any(Cmain,{'dE target'},3e-2);
NUM.dp_target           = Get_Value_Any(Cmain,{'dp target'},3e-2);
NUM.dmu_target          = Get_Value_Any(Cmain,{'dmu target'},1e2);

NUM.dt_good_count       = Get_Value_Any(Cmain,{'dt good count'},0);
NUM.dt_grow_after       = Get_Value_Any(Cmain,{'dt grow after'},8);
NUM.dt_grow_fac         = Get_Value_Any(Cmain,{'dt grow factor'},1.15);
NUM.dt_shrink_fac       = Get_Value_Any(Cmain,{'dt shrink factor'},0.5);
NUM.err_grow            = Get_Value_Any(Cmain,{'error grow'},0.25);

NUM.phi_mask_cut        = Get_Value_Any(Cmain,{'phi mask cutoff'},1e-8);
NUM.phi_mask_thick      = Get_Value_Any(Cmain,{'phi mask thickness grid'},3);
NUM.norm_phi            = Get_Value_Any(Cmain,{'normalize phi or not'},1);
NUM.cut_phi             = Get_Value_Any(Cmain,{'cut phi or not'},0);
NUM.norm_E              = Get_Value_Any(Cmain,{'normalize E or not'},1);
NUM.int_damp            = Get_Value_Any(Cmain,{'interface damping factor'},0.1);
NUM.use_Jphi            = Get_Value_Any(Cmain,{'use J_phi or not'},0);

NUM.CHLE_p_cut          = Get_Value_Any(Cmain,{'CH LE p cutoff'},1e-8);
NUM.CHLE_band_thick     = Get_Value_Any(Cmain,{'CH LE band thickness'},16);
NUM.CHLE_res_rel        = [];

NUM.linear_solver       = Get_Value_Any(Cmain,{'Linear solver method'},'bicgstab_ilu');
NUM.linear_tol          = Get_Value_Any(Cmain,{'Linear solver tolerance'},1e-8);
NUM.linear_maxit        = Get_Value_Any(Cmain,{'Linear solver max iteration'},500);
NUM.ilu_reuse           = Get_Value_Any(Cmain,{'Linear solver ilu reuse cache'},1);
NUM.ilu_reuse_steps     = Get_Value_Any(Cmain,{'Linear solver reuse steps'},5);
NUM.ilu_rebuild         = Get_Value_Any(Cmain,{'Linear solver ilu rebuild or not'},1);
NUM.direct_fallback     = Get_Value_Any(Cmain,{'Linear solver fallback'},1);

NUM.ilu_reuse_ACCH          = Get_Value_Any(Cmain,{'Ilu preconditioner reuse cache'},1);
NUM.ilu_reuse_steps_ACCH    = Get_Value_Any(Cmain,{'Ilu preconditioner reuse steps'},10);
NUM.ACCH_mask_update        = Get_Value_Any(Cmain,{'AC-CH solver mask update'},NUM.ilu_reuse_steps_ACCH);
NUM.ilu_cache_check_pattern = Get_Value_Any(Cmain,{'Ilu cache check pattern'},0);

% ------------------------------------------------------------
% P-T-t path
% ------------------------------------------------------------
[PARAM.PT.t_path,PARAM.PT.T_path,PARAM.PT.P_path] = Read_PT_Path(Cpt,PHYS);

NUM.t_tot = max(PARAM.PT.t_path);

% ------------------------------------------------------------
% MODEL
% ------------------------------------------------------------
MODEL.solmod = solmod;
MODEL.E_sc   = PHYS.E_sc;
MODEL.vref   = PHYS.vref;

if ~isfield(MODEL,'Cname') || isempty(MODEL.Cname)

    MODEL.Cname = elem_list;

    if ~any(strcmp(MODEL.Cname,'Si'))
        MODEL.Cname{end+1} = 'Si';
    end

    if ~any(strcmp(MODEL.Cname,'O'))
        MODEL.Cname{end+1} = 'O';
    end

end

eps_phi = 1e-14;

MODEL.dgdphi = @(phi) 2*PHYS.m*phi.*(phi - 1).^2 + PHYS.m*phi.^2.*(2.*phi - 2);

MODEL.p_fun  = @(a,phi) phi(:,:,a).^2./(sum(phi.^2,3) + eps_phi);

MODEL.dpdphi = @(a,b,phi) (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) - 2*phi(:,:,a).*phi(:,:,b).^2./(sum(phi.^2,3) + eps_phi).^2;

% ------------------------------------------------------------
% PARAM
% ------------------------------------------------------------
PARAM.Np            = length(STATE.c);
PARAM.Ne            = length(STATE.E);
PARAM.eta           = PHYS.eta;
PARAM.use_WScale    = 0;

PARAM.use_CS_chi    = Get_Value_Any(Cmain,{'Use convex split'},1);
PARAM.CS_chi_floor  = Get_Value_Any(Cmain,{'Convex split floor'},1e-9);
PARAM.CS_chi_cap    = Get_Value_Any(Cmain,{'Convex split cap'},1e-2);

PARAM.use_kappa_c   = Get_Value_Any(Cmain,{'Use 4th order term for c'},1);
PARAM.L_fac         = Get_Value_Any(Cmain,{'L scaling factor'},0.5);
PARAM.LE_mode       = Get_Value_Any(Cmain,{'LE solver (LE or GP)'},'LE');

PARAM.M_L_floor_fac      = Get_Value_Any(Cmain,{'M L floor factor'},1e-2);
PARAM.M_L_interface_only = Get_Value_Any(Cmain,{'M L interface only'},1);
PARAM.M_L_p_cut          = Get_Value_Any(Cmain,{'M L p cutoff'},1e-8);

PARAM.dceq          = PHYS.dceq;

PARAM.update_PT_every = Get_Value_Any(Cmain,{'Update PT every','PT update every'},1);

PARAM.update_PT_every = round(PARAM.update_PT_every);
if PARAM.update_PT_every < 1
    PARAM.update_PT_every = 1;
end

% Temporary initial values. Compute_M_And_L updates these during the run.
PARAM.L             = zeros(GRID.ny,GRID.nx);
PARAM.Lm            = zeros(GRID.ny,GRID.nx);
PARAM.LK            = zeros(GRID.ny,GRID.nx);

if isfield(MODEL,'pars')
    PARAM.kappa_phase = PHYS.kappa .* cellfun(@(x) size(x.n,1) > 1, MODEL.pars);
end

if ~isfield(PARAM,'theta_grain') || numel(PARAM.theta_grain) ~= PARAM.Np
    PARAM.theta_grain = zeros(1,PARAM.Np);
end

% ------------------------------------------------------------
% Phase mobility and anisotropy metadata
% ------------------------------------------------------------
Ne = PARAM.Ne;

if numel(elem_list) ~= Ne
    error('Read_PFM_Metadata: number of Excel mobility elements does not match STATE.E / PARAM.Ne.')
end

PHYS.M_phs_raw      = zeros(Nphase,Ne);
PHYS.M_phs          = zeros(Nphase,Ne);

PARAM.aniso_mode    = cell(1,Nphase);
PARAM.aniso_phase   = [];
PARAM.WScale_phase_factor = ones(1,Nphase);
facet_all           = struct([]);

for ip = 1:Nphase

    phase_name = MODEL.phs_name{ip};
    Cph        = readcell(xlsx_file,'Sheet',phase_name);
    % Interface excess-energy damping factor.
    % 1 = original excess energy
    % 0 = excess energy fully removed at the strongest interface/tail limit
    wfac       = Get_Value_Any(Cph,{'Interface damp W factor'},1);
    wfac       = min(max(wfac,0),1);
    PARAM.WScale_phase_factor(ip) = wfac;

    % Mobility in Excel is the raw mobility used by the PF solver.
    % This mobility is the coefficient in:
    %   dE/dt = div( M * grad(mu_e) )
    M_raw      = Read_Phase_Mobility(Cph,elem_list);

    PHYS.M_phs_raw(ip,:) = M_raw;
    PHYS.M_phs(ip,:)     = M_raw/PHYS.L_sc^2*PHYS.t_sc*PHYS.E_sc;

    mode = lower(Get_Value_Any(Cph,{'Anisotropy mode'},'iso'));
    PARAM.aniso_mode{ip} = mode;

    if ~strcmpi(mode,'iso')
        PARAM.aniso_phase(end+1) = ip; %#ok<AGROW>
    end

    facet_all(ip).phase_name = phase_name;
    facet_all(ip).mode       = mode;
    facet_all(ip).nfold      = Get_Value_Any(Cph,{'Anisotropy n fold'},[]);
    facet_all(ip).q          = Get_Value_Any(Cph,{'Anisotropy q'},0.25);
    facet_all(ip).amin       = Get_Value_Any(Cph,{'Anisotropy min'},0.2);
    facet_all(ip).amax       = Get_Value_Any(Cph,{'Anisotropy max'},10);

    [hkl,theta_deg,A_weight] = Read_Facet_Table(Cph);

    facet_all(ip).hkl       = hkl;
    facet_all(ip).theta_deg = theta_deg;
    facet_all(ip).theta     = theta_deg*pi/180;
    facet_all(ip).A_weight  = A_weight;
    facet_all(ip).A         = A_weight;

    if ~isempty(A_weight)
        facet_all(ip).A_ref = max(A_weight);
    else
        facet_all(ip).A_ref = 1.0;
    end

    facet_all(ip).sigma_ref = 1.0;

end

PHYS.M0 = max(PHYS.M_phs(:));

if PHYS.M0 <= 0
    error('Read_PFM_Metadata: PHYS.M0 is zero. Check mobility tables.')
end

PHYS.L = 4*PHYS.m/3/PHYS.kap/(PHYS.dceq^2/PHYS.M0);

PARAM.use_aniso = 0;

if ~isempty(PARAM.aniso_phase)
    PARAM.use_aniso = 1;
    PARAM.facet     = facet_all;
end

% ------------------------------------------------------------
% Scaling diagnostics
% ------------------------------------------------------------
fprintf('\nMetadata scaling check:\n')
fprintf('GRID.dx              = %.6e\n',GRID.dx)
fprintf('PHYS.l              = %.6e\n',PHYS.l)
fprintf('PHYS.sigma          = %.6e\n',PHYS.sigma)
fprintf('PHYS.kappa          = %.6e\n',PHYS.kappa)
fprintf('PHYS.M0             = %.6e\n',PHYS.M0)
fprintf('PHYS.m              = %.6e\n',PHYS.m)
fprintf('PHYS.kap            = %.6e\n',PHYS.kap)
fprintf('NUM.dt_phy          = %.6e\n',NUM.dt_phy)
fprintf('NUM.t_tot           = %.6e\n',NUM.t_tot)

end

%% ========================================================================
%  Local helper functions
% ========================================================================

function val = Get_Value_Any(C,keys,default)

val = default;

for ik = 1:numel(keys)

    key = strtrim(keys{ik});

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


function [t_path,T_path,P_path] = Read_PT_Path(Cpt,PHYS)

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
    error('Read_PFM_Metadata: PTt_Path header row not found.')
end

colP = [];
colT = [];
colt = [];

for j = 1:size(Cpt,2)

    txt = strtrim(Cell_String(Cpt{header_row,j}));

    if strcmpi(txt,'P (Pa)')
        colP = j;
    elseif strcmpi(txt,'T (K)')
        colT = j;
    elseif strcmpi(txt,'t (s)')
        colt = j;
    end

end

P_path = [];
T_path = [];
t_path = [];

for i = header_row+1:size(Cpt,1)

    if Is_Empty_Cell(Cpt{i,colP}) || Is_Empty_Cell(Cpt{i,colT}) || Is_Empty_Cell(Cpt{i,colt})
        continue
    end

    P_path(end+1) = Cell_Number(Cpt{i,colP}); %#ok<AGROW>
    T_path(end+1) = Cell_Number(Cpt{i,colT}); %#ok<AGROW>
    t_path(end+1) = Cell_Number(Cpt{i,colt}); %#ok<AGROW>

end

if isempty(t_path)
    error('Read_PFM_Metadata: PTt_Path contains no valid points.')
end

[t_path,ord] = sort(t_path);
T_path       = T_path(ord);
P_path       = P_path(ord);

if any(diff(t_path) <= 0)
    error('Read_PFM_Metadata: PTt_Path time must be strictly increasing.')
end

t_path = t_path/PHYS.t_sc;

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
    error('Read_PFM_Metadata: Mobility table is missing.')
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
        error('Read_PFM_Metadata: missing diagonal mobility for element %s.',elem)
    end

    M_diag(ie) = Cell_Number(Cph{row,col});

end

% Warn if off-diagonal values are nonzero.
for i = mob_row+1:size(Cph,1)

    if Is_Empty_Cell(Cph{i,1})
        continue
    end

    if strcmpi(strtrim(Cell_String(Cph{i,1})),'Facet')
        break
    end

    row_elem = strtrim(Cell_String(Cph{i,1}));

    for j = 2:size(Cph,2)

        col_elem = strtrim(Cell_String(Cph{mob_row,j}));

        if Is_Empty_Cell(Cph{i,j}) || isempty(col_elem)
            continue
        end

        val = Cell_Number(Cph{i,j});

        if ~strcmpi(row_elem,col_elem) && abs(val) > 0
            warning('Read_PFM_Metadata: off-diagonal mobility %s-%s is nonzero but ignored for now.', ...
                row_elem,col_elem)
        end

    end
end

end


function [hkl,theta_deg,A_weight] = Read_Facet_Table(Cph)

hkl       = {};
theta_deg = [];
A_weight  = [];

facet_row = [];

for i = 1:size(Cph,1)

    if strcmpi(strtrim(Cell_String(Cph{i,1})),'Facet')
        facet_row = i;
        break
    end

end

if isempty(facet_row)
    return
end

for i = facet_row+1:size(Cph,1)

    if Is_Empty_Cell(Cph{i,1})
        continue
    end

    hkl{end+1}       = strtrim(Cell_String(Cph{i,1})); %#ok<AGROW>
    theta_deg(end+1) = Cell_Number(Cph{i,2}); %#ok<AGROW>
    A_weight(end+1)  = Cell_Number(Cph{i,3}); %#ok<AGROW>

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
    error('Read_PFM_Metadata: expected numeric value.')
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
