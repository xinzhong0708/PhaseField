function [PHYS,NUM,PARAM,MODEL] = Read_PFM_Metadata(xlsx_file,GRID,MODEL,STATE,eta,PARAM)
%READ_PFM_METADATA Read PFM metadata from Excel.
%
% Required sheets:
%   Main
%   PTt_Path
%   one sheet for each phase in MODEL.phs_name
%
% Mobility table is interpreted directly as code mobility M.
% No D-to-M conversion is done here.

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
solmod    = Get_Value(Cmain,'Solution name','solution_models_PFM');
vref      = Get_Value(Cmain,'Molar volume',2e-5);

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
% PHYS
% ------------------------------------------------------------
PHYS             = struct();

PHYS.E_sc        = Get_Value(Cmain,'Energy scale',1e9);
PHYS.L_sc        = Get_Value(Cmain,'Length scale',1e-6);
PHYS.t_sc        = Get_Value(Cmain,'Time scale',1);
PHYS.vref        = vref;

PHYS.chi_ref     = Get_Value(Cmain,'chi_ref',1e-2);
PHYS.dceq        = Get_Value(Cmain,'dceq',0.5);

thick_fac        = Get_Value(Cmain,'Interface thick factor',3);
PHYS.l           = thick_fac*GRID.dx/PHYS.L_sc;

kappa_phys       = Get_Value(Cmain,'4th order kappa',0);
PHYS.kappa       = kappa_phys/(PHYS.E_sc*PHYS.L_sc^2);

sigma_phase_Jm2  = zeros(1,Nphase);

for ip = 1:Nphase

    Cph = readcell(xlsx_file,'Sheet',MODEL.phs_name{ip});
    sigma_phase_Jm2(ip) = Get_Value(Cph,'Interface energy',0.5);

end

PHYS.sigma_phase_Jm2 = sigma_phase_Jm2;

% Current solver still uses one scalar sigma.
% For now use the mean phase value.
PHYS.sigma      = mean(sigma_phase_Jm2)/(PHYS.E_sc*PHYS.L_sc);

PHYS.m          = 6*PHYS.sigma/PHYS.l;
PHYS.kap        = 3/4*PHYS.sigma*PHYS.l;
PHYS.eta        = MODEL.eta;

% ------------------------------------------------------------
% NUM
% ------------------------------------------------------------
NUM                     = struct();

NUM.dt_phy              = Get_Value(Cmain,'dt initial',1e-2)/PHYS.t_sc;
NUM.dt_max              = Get_Value(Cmain,'dt max',100)/PHYS.t_sc;
NUM.dt_min              = Get_Value(Cmain,'dt min',1e-15)/PHYS.t_sc;

NUM.t_phy               = 0;
NUM.time                = 0;

NUM.dE_target           = Get_Value(Cmain,'dE target',3e-2);
NUM.dp_target           = Get_Value(Cmain,'dp target',3e-2);
NUM.dmu_target          = Get_Value(Cmain,'dmu target',1e2);

NUM.dt_good_count       = Get_Value(Cmain,'dt good count',0);
NUM.dt_grow_after       = Get_Value(Cmain,'dt grow after',8);
NUM.dt_grow_fac         = Get_Value(Cmain,'dt grow factor',1.15);
NUM.dt_shrink_fac       = Get_Value(Cmain,'dt shrink factor',0.5);
NUM.err_grow            = Get_Value(Cmain,'error grow',0.25);

NUM.phi_mask_cut        = Get_Value(Cmain,'phi mask cutoff',1e-8);
NUM.phi_mask_thick      = Get_Value(Cmain,'phi mask thickness grid',3);
NUM.norm_phi            = Get_Value(Cmain,'normalize phi or not',1);
NUM.cut_phi             = Get_Value(Cmain,'cut phi or not',0);
NUM.norm_E              = Get_Value(Cmain,'normalize E or not',1);
NUM.int_damp            = Get_Value(Cmain,'interface damping factor',0.1);
NUM.use_Jphi            = Get_Value(Cmain,'use J_phi or not',0);

NUM.CHLE_p_cut          = Get_Value(Cmain,'CH LE p cutoff',1e-8);
NUM.CHLE_band_thick     = Get_Value(Cmain,'CH LE band thickness',16);
NUM.CHLE_res_rel        = [];

NUM.linear_solver       = Get_Value(Cmain,'Linear solver method','bicgstab_ilu');
NUM.linear_tol          = Get_Value(Cmain,'Linear solver tolerance',1e-8);
NUM.linear_maxit        = Get_Value(Cmain,'Linear solver max iteration',500);
NUM.ilu_reuse           = Get_Value(Cmain,'Linear solver ilu reuse cache',1);
NUM.ilu_reuse_steps     = Get_Value(Cmain,'Linear solver reuse steps',5);
NUM.ilu_rebuild         = Get_Value(Cmain,'Linear solver ilu rebuild or not',1);
NUM.direct_fallback     = Get_Value(Cmain,'Linear solver fallback',1);

NUM.ilu_reuse_ACCH          = Get_Value(Cmain,'Ilu preconditioner reuse cache',1);
NUM.ilu_reuse_steps_ACCH    = Get_Value(Cmain,'Ilu preconditioner reuse steps',10);
NUM.ACCH_mask_update        = Get_Value(Cmain,'AC-CH solver mask update',NUM.ilu_reuse_steps_ACCH);
NUM.ilu_cache_check_pattern = Get_Value(Cmain,'Ilu cache check pattern',0);

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

MODEL.p_fun  = @(a,phi) ...
    phi(:,:,a).^2./(sum(phi.^2,3) + eps_phi);

MODEL.dpdphi = @(a,b,phi) ...
    (a==b)*2*phi(:,:,b)./(sum(phi.^2,3) + eps_phi) ...
    - 2*phi(:,:,a).*phi(:,:,b).^2./(sum(phi.^2,3) + eps_phi).^2;

% ------------------------------------------------------------
% PARAM
% ------------------------------------------------------------
PARAM.Np            = length(STATE.c);
PARAM.Ne            = length(STATE.E);
PARAM.eta           = MODEL.eta;
PARAM.use_WScale    = 0;

PARAM.use_CS_chi    = Get_Value(Cmain,'Use convex split',1);
PARAM.CS_chi_floor  = Get_Value(Cmain,'Convex split floor',1e-9);
PARAM.CS_chi_cap    = Get_Value(Cmain,'Convex split cap',1e-2);

PARAM.use_kappa_c   = Get_Value(Cmain,'Use 4th order term for c',1);
PARAM.L_fac         = Get_Value(Cmain,'L scaling factor',0.5);
PARAM.LE_mode       = Get_Value(Cmain,'LE solver (LE or GP)','LE');

PARAM.M_L_floor_fac      = Get_Value(Cmain,'M L floor factor',1e-2);
PARAM.M_L_interface_only = Get_Value(Cmain,'M L interface only',1);
PARAM.M_L_p_cut          = Get_Value(Cmain,'M L p cutoff',1e-8);

PARAM.dceq          = PHYS.dceq;

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

PHYS.M_phs          = zeros(Nphase,Ne);
PARAM.aniso_mode    = cell(1,Nphase);
PARAM.aniso_phase   = [];
facet_all           = struct([]);

for ip = 1:Nphase

    phase_name = MODEL.phs_name{ip};
    Cph        = readcell(xlsx_file,'Sheet',phase_name);

    % Mobility is loaded directly as code mobility M.
    PHYS.M_phs(ip,:) = Read_Phase_Mobility(Cph,elem_list);

    mode = lower(Get_Value(Cph,'Anisotropy mode','iso'));
    PARAM.aniso_mode{ip} = mode;

    if ~strcmpi(mode,'iso')
        PARAM.aniso_phase(end+1) = ip; %#ok<AGROW>
    end

    facet_all(ip).phase_name = phase_name;
    facet_all(ip).mode       = mode;
    facet_all(ip).nfold      = Get_Value(Cph,'Anisotropy n fold',[]);
    facet_all(ip).q          = Get_Value(Cph,'Anisotropy q',0.25);
    facet_all(ip).amin       = Get_Value(Cph,'Anisotropy min',0.2);
    facet_all(ip).amax       = Get_Value(Cph,'Anisotropy max',10);

    [hkl,theta_deg,A_weight] = Read_Facet_Table(Cph);

    facet_all(ip).hkl       = hkl;
    facet_all(ip).theta_deg = theta_deg;
    facet_all(ip).theta     = theta_deg*pi/180;
    facet_all(ip).A_weight  = A_weight;

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

end

%% ========================================================================
%  Local helper functions
% ========================================================================

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