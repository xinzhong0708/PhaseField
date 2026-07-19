function [PHYS,NUM,PARAM,MODEL,GRID] = Read_PFM_Metadata(xlsx_file,GRID,MODEL,STATE,eta,PARAM)
%READ_PFM_METADATA Read PFM metadata from Excel and scale to code units.
%
% Excel metadata is assumed to use SI / physical units.
%
% Required sheets:
%   Main
%   PTt_Path
%   Interface
%   one sheet for each phase in MODEL.phs_name
%
% The Interface sheet contains the pairwise interface-energy matrix:
%
%                  Phase-1  Phase-2 ...
%   Phase-1        sigma11  sigma12 ...
%   Phase-2        sigma21  sigma22 ...
%
% Values are in J/m^2. Rows and columns are matched by phase name, so their
% order does not need to match MODEL.phs_name.
%
% Pairwise scaling:
%   sigma_code = sigma_SI/(E_sc*L_sc)
%   m_pair     = 6*sigma_pair/l
%   kap_pair   = 3/4*sigma_pair*l
%
% The diagonal matrix entry is the interface energy between two different
% grains of the same thermodynamic phase. The self-pair of one grain with
% itself is set to zero later in Compute_M_And_L.
%
% MODEL.dgdphi retains one scalar reference m. The actual pairwise barrier
% coefficient used by the AC solver is PARAM.Lm_AC, built by
% Compute_M_And_L. Calc_S_AllenCahn divides out the reference m.

if nargin < 6 || isempty(PARAM)
    PARAM = struct();
end

% -------------------------------------------------------------------------
% Read sheets
% -------------------------------------------------------------------------
Cmain      = readcell(xlsx_file,'Sheet','Main');
Cpt        = readcell(xlsx_file,'Sheet','PTt_Path');
Cint       = readcell(xlsx_file,'Sheet','Interface');
sheet_list = sheetnames(xlsx_file);

if ~isfield(MODEL,'phs_name') || isempty(MODEL.phs_name)
    error(['Read_PFM_Metadata: MODEL.phs_name is missing. ', ...
           'Active phases must come from the map.'])
end

MODEL.phs_name = Cellstr_Row(MODEL.phs_name);
Nphase         = numel(MODEL.phs_name);

% -------------------------------------------------------------------------
% Main metadata
% -------------------------------------------------------------------------
elem_list = Get_Row_List(Cmain,'Independent c (Si dependent)');
solmod    = Get_Value_Any(Cmain,{'Solution name','Solution model'}, ...
                          'solution_models_PFM');
vref      = Get_Value_Any(Cmain,{'Molar volume','vref'},2e-5);

if isempty(elem_list)
    error('Read_PFM_Metadata: element list is missing in Main sheet.')
end

% -------------------------------------------------------------------------
% Check required sheets
% -------------------------------------------------------------------------
required_sheet = [{'Main','PTt_Path','Interface'}, MODEL.phs_name];

for is = 1:numel(required_sheet)

    if ~any(strcmpi(sheet_list,required_sheet{is}))
        error('Read_PFM_Metadata: missing sheet "%s".',required_sheet{is})
    end

end

% -------------------------------------------------------------------------
% PHYS: scales
% -------------------------------------------------------------------------
PHYS             = struct();

PHYS.E_sc        = Get_Value_Any(Cmain,{'Energy scale','Energy scaling'},1e9);
PHYS.L_sc        = Get_Value_Any(Cmain,{'Length scale','Length scaling'},1e-6);
PHYS.t_sc        = Get_Value_Any(Cmain,{'Time scale','Time scaling'},1);
PHYS.vref        = vref;
PHYS.dceq        = Get_Value_Any(Cmain,{'dceq','dc eq'},0.5);

if ~isscalar(PHYS.E_sc) || PHYS.E_sc <= 0 || ~isfinite(PHYS.E_sc)
    error('Read_PFM_Metadata: Energy scale must be positive.')
end
if ~isscalar(PHYS.L_sc) || PHYS.L_sc <= 0 || ~isfinite(PHYS.L_sc)
    error('Read_PFM_Metadata: Length scale must be positive.')
end
if ~isscalar(PHYS.t_sc) || PHYS.t_sc <= 0 || ~isfinite(PHYS.t_sc)
    error('Read_PFM_Metadata: Time scale must be positive.')
end

% -------------------------------------------------------------------------
% Scale GRID once
% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
% Interface thickness and composition-gradient kappa
% -------------------------------------------------------------------------
thick_fac = Get_Value_Any(Cmain, ...
    {'Interface thick factor','Interface thickness factor'},3);

if ~isscalar(thick_fac) || thick_fac <= 0 || ~isfinite(thick_fac)
    error('Read_PFM_Metadata: interface thickness factor must be positive.')
end

PHYS.l = thick_fac*GRID.dx;

kappa_SI      = Get_Value_Any(Cmain,{'4th order kappa','kappa'},0);
PHYS.kappa_SI = kappa_SI;
PHYS.kappa    = kappa_SI/(PHYS.E_sc*PHYS.L_sc^2);

% -------------------------------------------------------------------------
% Pairwise interface energy
% -------------------------------------------------------------------------
sigma_pair_SI = Read_Interface_Matrix(Cint,MODEL.phs_name);

% Validate symmetry before scaling.
sigma_scale = max(1,max(abs(sigma_pair_SI(:))));
sym_err     = max(abs(sigma_pair_SI - sigma_pair_SI.'),[],'all');

if sym_err > 1e-10*sigma_scale
    error(['Read_PFM_Metadata: Interface matrix is not symmetric. ', ...
           'Maximum mismatch = %.6e J/m^2.'],sym_err)
end

sigma_pair_SI = 0.5*(sigma_pair_SI + sigma_pair_SI.');

if any(~isfinite(sigma_pair_SI(:))) || any(sigma_pair_SI(:) <= 0)
    error(['Read_PFM_Metadata: all Interface matrix entries must be ', ...
           'finite and strictly positive.'])
end

PHYS.interface_phase_names = MODEL.phs_name;
PHYS.sigma_pair_SI         = sigma_pair_SI;
PHYS.sigma_pair            = sigma_pair_SI/(PHYS.E_sc*PHYS.L_sc);

PHYS.m_pair   = 6*PHYS.sigma_pair/PHYS.l;
PHYS.kap_pair = (3/4)*PHYS.sigma_pair*PHYS.l;

% A positive reference is required only because MODEL.dgdphi contains m.
% The reference cancels in Calc_S_AllenCahn and does not replace pairwise m.
upper_mask        = triu(true(Nphase));
sigma_unique      = PHYS.sigma_pair(upper_mask);
sigma_unique      = sigma_unique(isfinite(sigma_unique) & sigma_unique > 0);
PHYS.sigma        = median(sigma_unique);
PHYS.sigma_ref    = PHYS.sigma;
PHYS.sigma_ref_SI = PHYS.sigma*PHYS.E_sc*PHYS.L_sc;
PHYS.m            = 6*PHYS.sigma_ref/PHYS.l;
PHYS.kap          = (3/4)*PHYS.sigma_ref*PHYS.l;

% -------------------------------------------------------------------------
% Global W/kappa interface-ramp geometry
% -------------------------------------------------------------------------
zero_grid = Get_Value_Any(Cmain, ...
    {'Interface zero kappa grid','Interface zero kappa grid'},2);
ramp_grid = Get_Value_Any(Cmain, ...
    {'Interface ramp grid','Interface kappa ramp grid'},6);

zero_grid = max(0,round(zero_grid));
ramp_grid = max(1,round(ramp_grid));

PARAM.Interface_zero_kapp_grid = zero_grid;
PARAM.Interface_ramp_grid      = ramp_grid;
PARAM.ramp_zero_width          = zero_grid;
PARAM.ramp_width               = ramp_grid;

% -------------------------------------------------------------------------
% Eta
% -------------------------------------------------------------------------
eta_SI = Get_Value_Any(Cmain,{'Penalty eta','eta'},[]);

if isempty(eta_SI)

    if isfield(PARAM,'eta_bulk') && ~isempty(PARAM.eta_bulk)
        PHYS.eta = PARAM.eta_bulk;
    else
        PHYS.eta = eta;
    end

else

    PHYS.eta = eta_SI/PHYS.E_sc*ones(GRID.ny,GRID.nx);

end

% -------------------------------------------------------------------------
% NUM
% -------------------------------------------------------------------------
NUM                     = struct();

NUM.dt_phy              = Get_Value_Any(Cmain, ...
    {'dt initial','Initial dt'},1e-2)/PHYS.t_sc;
NUM.dt_max              = Get_Value_Any(Cmain, ...
    {'dt max','Maximum dt'},100)/PHYS.t_sc;
NUM.dt_min              = Get_Value_Any(Cmain, ...
    {'dt min','Minimum dt'},1e-15)/PHYS.t_sc;

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
NUM.phi_mask_thick      = Get_Value_Any(Cmain, ...
    {'phi mask thickness grid'},3);
NUM.norm_phi            = Get_Value_Any(Cmain,{'normalize phi or not'},1);
NUM.cut_phi             = Get_Value_Any(Cmain,{'cut phi or not'},0);
NUM.norm_E              = Get_Value_Any(Cmain,{'normalize E or not'},1);
NUM.int_damp            = Get_Value_Any(Cmain, ...
    {'interface damping factor'},0.1);
NUM.use_Jphi            = Get_Value_Any(Cmain,{'use J_phi or not'},0);

NUM.CHLE_p_cut          = Get_Value_Any(Cmain,{'CH LE p cutoff'},1e-8);
NUM.CHLE_band_thick     = Get_Value_Any(Cmain, ...
    {'CH LE band thickness'},16);
NUM.CHLE_res_rel        = [];

NUM.linear_solver       = Get_Value_Any(Cmain, ...
    {'Linear solver method'},'bicgstab_ilu');
NUM.linear_tol          = Get_Value_Any(Cmain, ...
    {'Linear solver tolerance'},1e-8);
NUM.linear_maxit        = Get_Value_Any(Cmain, ...
    {'Linear solver max iteration'},500);
NUM.ilu_reuse           = Get_Value_Any(Cmain, ...
    {'Linear solver ilu reuse cache'},1);
NUM.ilu_reuse_steps     = Get_Value_Any(Cmain, ...
    {'Linear solver reuse steps'},5);
NUM.ilu_rebuild         = Get_Value_Any(Cmain, ...
    {'Linear solver ilu rebuild or not'},1);
NUM.direct_fallback     = Get_Value_Any(Cmain, ...
    {'Linear solver fallback'},1);

NUM.ilu_reuse_ACCH          = Get_Value_Any(Cmain, ...
    {'Ilu preconditioner reuse cache'},1);
NUM.ilu_reuse_steps_ACCH    = Get_Value_Any(Cmain, ...
    {'Ilu preconditioner reuse steps'},10);
NUM.ACCH_mask_update        = Get_Value_Any(Cmain, ...
    {'AC-CH solver mask update'},NUM.ilu_reuse_steps_ACCH);
NUM.ilu_cache_check_pattern = Get_Value_Any(Cmain, ...
    {'Ilu cache check pattern'},0);

% -------------------------------------------------------------------------
% P-T-t path
% -------------------------------------------------------------------------
[PARAM.PT.t_path,PARAM.PT.T_path,PARAM.PT.P_path] = ...
    Read_PT_Path(Cpt,PHYS);

NUM.t_tot = max(PARAM.PT.t_path);

% -------------------------------------------------------------------------
% MODEL
% -------------------------------------------------------------------------
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

% Reference double-well derivative. Pairwise physical m is supplied later
% through PARAM.Lm_AC and the reference is divided out in Calc_S_AllenCahn.
MODEL.dgdphi = @(phi) ...
    2*PHYS.m*phi.*(phi - 1).^2 + ...
    PHYS.m*phi.^2.*(2.*phi - 2);

MODEL.p_fun  = @(a,phi) ...
    phi(:,:,a).^2./(sum(phi.^2,3) + eps_phi);

MODEL.dpdphi = @(a,b,phi) ...
    (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) - ...
    2*phi(:,:,a).*phi(:,:,b).^2./ ...
    (sum(phi.^2,3) + eps_phi).^2;

% -------------------------------------------------------------------------
% PARAM
% -------------------------------------------------------------------------
PARAM.Np         = length(STATE.c);
PARAM.Ne         = length(STATE.E);
PARAM.eta        = PHYS.eta;
PARAM.use_WScale = 0;

PARAM.use_CS_chi   = Get_Value_Any(Cmain,{'Use convex split'},1);
PARAM.CS_chi_floor = Get_Value_Any(Cmain,{'Convex split floor'},1e-9);
PARAM.CS_chi_cap   = Get_Value_Any(Cmain,{'Convex split cap'},1e-2);

PARAM.use_kappa_c = Get_Value_Any(Cmain, ...
    {'Use 4th order term for c'},1);
PARAM.L_fac       = Get_Value_Any(Cmain,{'L scaling factor'},0.5);
PARAM.LE_mode     = Get_Value_Any(Cmain, ...
    {'LE solver (LE or GP)'},'LE');

PARAM.M_L_floor_fac      = Get_Value_Any(Cmain, ...
    {'M L floor factor'},1e-2);
PARAM.M_L_interface_only = Get_Value_Any(Cmain, ...
    {'M L interface only'},1);
PARAM.M_L_p_cut          = Get_Value_Any(Cmain, ...
    {'M L p cutoff'},1e-8);

PARAM.interface_pair_p_cut = Get_Value_Any(Cmain, ...
    {'Interface pair p cutoff','Pair interface p cutoff'},1e-12);

PARAM.dceq = PHYS.dceq;

PARAM.update_PT_every = Get_Value_Any(Cmain, ...
    {'Update PT every','PT update every'},1);
PARAM.update_PT_every = round(PARAM.update_PT_every);

if PARAM.update_PT_every < 1
    PARAM.update_PT_every = 1;
end

% Pairwise phase-level coefficients in code units.
PARAM.interface_phase_names = MODEL.phs_name;
PARAM.sigma_pair_phase_SI   = PHYS.sigma_pair_SI;
PARAM.sigma_pair_phase      = PHYS.sigma_pair;
PARAM.m_pair_phase          = PHYS.m_pair;
PARAM.kap_pair_phase        = PHYS.kap_pair;
PARAM.m_dgdphi_ref          = PHYS.m;
PARAM.m_AC_base             = PHYS.m;

% Temporary fields. Compute_M_And_L replaces them each timestep.
PARAM.L  = zeros(GRID.ny,GRID.nx);
PARAM.Lm = zeros(GRID.ny,GRID.nx);
PARAM.LK = zeros(GRID.ny,GRID.nx);

PARAM.L_AC  = zeros(GRID.ny,GRID.nx,PARAM.Np);
PARAM.Lm_AC = zeros(GRID.ny,GRID.nx,PARAM.Np);
PARAM.LK_AC = zeros(GRID.ny,GRID.nx,PARAM.Np);

if isfield(MODEL,'pars')
    PARAM.kappa_phase = PHYS.kappa .* ...
        cellfun(@(x) size(x.n,1) > 1,MODEL.pars);
end

if ~isfield(PARAM,'theta_grain') || numel(PARAM.theta_grain) ~= PARAM.Np
    PARAM.theta_grain = zeros(1,PARAM.Np);
end

% -------------------------------------------------------------------------
% Phase mobility and anisotropy metadata
% -------------------------------------------------------------------------
Ne = PARAM.Ne;

if numel(elem_list) ~= Ne
    error(['Read_PFM_Metadata: number of Excel mobility elements does ', ...
           'not match STATE.E / PARAM.Ne.'])
end

PHYS.M_phs_raw = zeros(Nphase,Ne);
PHYS.M_phs     = zeros(Nphase,Ne);

PARAM.aniso_mode = cell(1,Nphase);
PARAM.aniso_phase = [];
PARAM.WScale_phase_factor = ones(1,Nphase);
facet_all = struct([]);

for ip = 1:Nphase

    phase_name = MODEL.phs_name{ip};
    Cph        = readcell(xlsx_file,'Sheet',phase_name);

    wfac = Get_Value_Any(Cph,{'Interface damp W factor'},1);
    wfac = min(max(wfac,0),1);
    PARAM.WScale_phase_factor(ip) = wfac;

    M_raw = Read_Phase_Mobility(Cph,elem_list);

    PHYS.M_phs_raw(ip,:) = M_raw;
    PHYS.M_phs(ip,:)     = ...
        M_raw/PHYS.L_sc^2*PHYS.t_sc*PHYS.E_sc;

    mode = lower(Get_Value_Any(Cph,{'Anisotropy mode'},'iso'));
    PARAM.aniso_mode{ip} = mode;

    if ~strcmpi(mode,'iso')
        PARAM.aniso_phase(end+1) = ip; %#ok<AGROW>
    end

    facet_all(ip).phase_name = phase_name;
    facet_all(ip).mode       = mode;
    facet_all(ip).nfold      = Get_Value_Any(Cph, ...
        {'Anisotropy n fold'},[]);
    facet_all(ip).q          = Get_Value_Any(Cph, ...
        {'Anisotropy q'},0.25);
    facet_all(ip).amin       = Get_Value_Any(Cph, ...
        {'Anisotropy min'},0.2);
    facet_all(ip).amax       = Get_Value_Any(Cph, ...
        {'Anisotropy max'},10);

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

% -------------------------------------------------------------------------
% Diagnostics
% -------------------------------------------------------------------------
fprintf('\nMetadata scaling check:\n')
fprintf('GRID.dx                    = %.6e\n',GRID.dx)
fprintf('PHYS.l                    = %.6e\n',PHYS.l)
fprintf('sigma_pair_SI min/max     = %.6e  %.6e J/m2\n', ...
    min(PHYS.sigma_pair_SI,[],'all'), ...
    max(PHYS.sigma_pair_SI,[],'all'))
fprintf('sigma reference           = %.6e J/m2\n',PHYS.sigma_ref_SI)
fprintf('m_pair min/max            = %.6e  %.6e\n', ...
    min(PHYS.m_pair,[],'all'),max(PHYS.m_pair,[],'all'))
fprintf('kap_pair min/max          = %.6e  %.6e\n', ...
    min(PHYS.kap_pair,[],'all'),max(PHYS.kap_pair,[],'all'))
fprintf('PHYS.kappa                = %.6e\n',PHYS.kappa)
fprintf('PHYS.M0                   = %.6e\n',PHYS.M0)
fprintf('NUM.dt_phy                = %.6e\n',NUM.dt_phy)
fprintf('NUM.t_tot                 = %.6e\n',NUM.t_tot)

end


%% =========================================================================
% Local helper functions
% =========================================================================
function sigma_pair = Read_Interface_Matrix(Cint,phase_names)
%READ_INTERFACE_MATRIX Read pairwise sigma by matching phase names.

phase_names = Cellstr_Row(phase_names);
Nphase      = numel(phase_names);

header_row = [];
col_idx    = zeros(1,Nphase);

% Find a row that contains every requested phase name.
for ir = 1:size(Cint,1)

    col_try = zeros(1,Nphase);

    for ip = 1:Nphase

        % Pair-column names must be in data columns, not the row-label column.
        for jc = 2:size(Cint,2)

            if strcmpi(strtrim(Cell_String(Cint{ir,jc})),phase_names{ip})
                col_try(ip) = jc;
                break
            end

        end

    end

    if all(col_try > 0) && numel(unique(col_try)) == Nphase
        header_row = ir;
        col_idx    = col_try;
        break
    end

end

if isempty(header_row)
    error(['Read_PFM_Metadata: Interface sheet does not contain a ', ...
           'header row with every active MODEL.phs_name.'])
end

row_idx = zeros(1,Nphase);

for ip = 1:Nphase

    matches = [];

    for ir = header_row+1:size(Cint,1)

        % Row labels are expected in the columns before the first data column.
        for jc = 1:max(1,min(col_idx)-1)

            if strcmpi(strtrim(Cell_String(Cint{ir,jc})),phase_names{ip})
                matches(end+1) = ir; %#ok<AGROW>
                break
            end

        end

    end

    matches = unique(matches);

    if isempty(matches)
        error('Read_PFM_Metadata: Interface row "%s" is missing.', ...
              phase_names{ip})
    elseif numel(matches) > 1
        error('Read_PFM_Metadata: Interface row "%s" is duplicated.', ...
              phase_names{ip})
    end

    row_idx(ip) = matches;

end

sigma_pair = zeros(Nphase,Nphase);

for ip = 1:Nphase
    for jp = 1:Nphase

        c = Cint{row_idx(ip),col_idx(jp)};

        if Is_Empty_Cell(c)
            error(['Read_PFM_Metadata: Interface value for %s-%s ', ...
                   'is empty.'],phase_names{ip},phase_names{jp})
        end

        sigma_pair(ip,jp) = Cell_Number(c);

    end
end

end


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

    if Is_Empty_Cell(Cpt{i,colP}) || ...
       Is_Empty_Cell(Cpt{i,colT}) || ...
       Is_Empty_Cell(Cpt{i,colt})
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
    col  = [];
    row  = [];

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
        error('Read_PFM_Metadata: missing diagonal mobility for element %s.', ...
              elem)
    end

    M_diag(ie) = Cell_Number(Cph{row,col});

end

% Warn if off-diagonal mobility values are present.
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
            warning(['Read_PFM_Metadata: off-diagonal mobility %s-%s ', ...
                     'is nonzero but ignored for now.'], ...
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


function c = Cellstr_Row(c)

if isstring(c)
    c = cellstr(c);
elseif ischar(c)
    c = {c};
end

c = c(:).';

for i = 1:numel(c)
    c{i} = strtrim(Cell_String(c{i}));
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

if ~isscalar(x) || isnan(x)
    error('Read_PFM_Metadata: expected one numeric value.')
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
