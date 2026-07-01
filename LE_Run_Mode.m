function [STATE] = LE_Run_Mode(STATE,PARAM,MODEL)
%This function perform local equilibrium calculation
%Two modes are available
%1) 'LE': Classical KKS type fixed E mode with fixed p and E
%2) 'GP': Grand potential based local equilibrium with fixed mu_e
%
%Stability update in this version:
%There is no ghost return-map. Tiny phases can be active at a low p, but
%their composition is written back with a p-weighted relaxation and a dc cap.
%This avoids low-p c/omega spikes while keeping phases thermodynamically visible.

%Choose the mode
mode       =  'LE';
if isfield(PARAM,'LE_mode')
    mode   = PARAM.LE_mode;
end

%p mode on and off
%These thresholds control thermodynamic phase activation for LE. They can be
%lower than before because low-p c is no longer written back abruptly.
sc         =  0.01;
p_tail     =  sc*1e-2;
p_full     =  sc*1e-1;
p_on       =  sc*5e-2;
p_off      =  sc*2e-2;

if isfield(PARAM,'LE_p_tail'), p_tail = PARAM.LE_p_tail; end
if isfield(PARAM,'LE_p_full'), p_full = PARAM.LE_p_full; end
if isfield(PARAM,'LE_p_on'),   p_on   = PARAM.LE_p_on;   end
if isfield(PARAM,'LE_p_off'),  p_off  = PARAM.LE_p_off;  end

%Maximal number of active phases for LE/GP
Pmax       =  3;
if isfield(PARAM,'LE_Pmax'), Pmax = PARAM.LE_Pmax; end

%Local equilibrium parameters
alpha_LE   = [0.6 0.4 0.3 0.2];
iter_LE    = [100 100 100 100];
iter_GP    = [100 100 100 100];

if isfield(PARAM,'LE_alpha_LE'), alpha_LE = PARAM.LE_alpha_LE; end
if isfield(PARAM,'LE_iter_LE'),  iter_LE  = PARAM.LE_iter_LE;  end
if isfield(PARAM,'LE_iter_GP'),  iter_GP  = PARAM.LE_iter_GP;  end

%Unpack fields and prepare the calculation
ny         =  size(STATE.p,1);
nx         =  size(STATE.p,2);
N          =  ny*nx;
p          =  UnpackP(STATE.p);
c          =  UnpackC(STATE.c);
E          =  UnpackCell(STATE.E);
mu_e       =  UnpackCell(STATE.mu_e);
chi        =  UnpackChi(STATE.chi);
pars       =  MODEL.pars;

%Prepare eta
if isscalar(PARAM.eta)
    eta    =  PARAM.eta*ones(1,N);
else
    eta    =  reshape(PARAM.eta,1,[]);
end

%Keep grain-resolved fields before collapse
p_grain    =  p;
c_grain    =  c;

%Find the active phase with smoothed thermodynamic p and hysteresis
if isfield(STATE,'LE_state') && ~isempty(STATE.LE_state)
    LE_state_old = STATE.LE_state;
else
    LE_state_old = struct();
end


%Phase index of the grains
%phase_index(ig) is the thermodynamic phase ID of grain ig.
phase_index = MODEL.phase_index(:)';

%Build phase-level p only for active-set detection.
[phase_id,grain_to_phase,pars_phase,p_phase] = BuildPhaseInfo(pars,p_grain,phase_index);

%Number of thermodynamic phases after grouping grains
Np   = size(p_phase,3);
Ne   = numel(E);
Pmax = min(Pmax,Np);

%Representative pars for each thermodynamic phase
pars = pars_phase;

%Damping of the interface excess energy, if used.
if isfield(PARAM,'w_scale_phase') && ~isempty(PARAM.w_scale_phase)
    pars_inter = Apply_WScale_Map(pars,PARAM.w_scale_phase,ny,nx);
else
    pars_inter = pars;
end

%Elemental composition of grain-resolved fields.
%This is used to subtract inactive grains from E_loc.
e_grain = Calc_e(MODEL.pars,c_grain);


%Find active thermodynamic phases using phase-collapsed p
if isfield(STATE,'LE_state') && ~isempty(STATE.LE_state)
    LE_state_old = STATE.LE_state;
else
    LE_state_old = struct();
end

[p_le,p_th,active] = BuildActiveSet_Tail(p_phase,LE_state_old,Pmax,p_tail,p_full,p_on,p_off);

n_active_map       = reshape(sum(active,2),ny,nx);

%Group grids by identical active thermodynamic phase set
[active_sets,~,set_id] = unique(active,'rows');
set_id_map             = reshape(set_id,ny,nx);
LE_diag                = InitLEDiag(size(active_sets,1),N);


%MAIN PART OF THE LE
for iset = 1:size(active_sets,1)

    %Active thermodynamic phases for this grid group
    ph_act = find(active_sets(iset,:));
    k      = numel(ph_act);

    %Linear grid ids belonging to this active set
    ids = find(set_id == iset).';

    if isempty(ids)
        continue
    end

    %Build local collapsed problem only for these ids
    if isfield(PARAM,'use_WScale') && PARAM.use_WScale
        pars_use = pars_inter;
    else
        pars_use = pars;
    end

    [pars_loc,c_loc,p_loc,E_loc,eta_loc] = BuildLocalLEBlock(pars_use,c_grain,p_grain,e_grain,E,eta,phase_index,phase_id,ph_act,ids);

    %Choose damping based on number of active phases
    kk = min(k,numel(alpha_LE));
    aa = alpha_LE(kk);
    diag_loc = InitLocalLEDiag(numel(ids));

    %MODE 1: KKS TYPE LOCAL EQUILIBRIUM WITH E AND P FIXED
    if strcmpi(mode,'LE')
        mm                         = iter_LE(kk);
        [c_loc,mu_loc,chi_loc,diag_loc] = LE_Calculator(pars_loc,p_loc,c_loc,E_loc,eta_loc,[aa,mm]);
    %MODE 2: GRAND POTENTIAL METHOD WITH FIXED P AND MU_E
    elseif strcmpi(mode,'GP')
        if k == 1
            mm                     = iter_LE(kk);
            [c_loc,mu_loc,chi_loc,diag_loc] = LE_Calculator(pars_loc,p_loc,c_loc,E_loc,eta_loc,[aa,mm]);
        else
            mm                     = iter_GP(kk);
            mu_loc                 = SliceCellIds(mu_e,ids);
            [c_loc,chi_loc]        = GP_Calculator(pars_loc,p_loc,c_loc,mu_loc,E_loc,eta_loc,aa,mm);
        end
    else
        error('LE_Run_Mode: unknown LE_mode "%s".',mode)
    end

    LE_diag = AccumulateLEDiag(LE_diag,diag_loc,ids,iset);

    %Write returned phase-local c back to the real grain-resolved storage.
    c_grain = AssignLocalCToGrains(c_grain,p_grain,phase_index,phase_id,ph_act,ids,c_loc,PARAM,p_tail,p_full);

    %mu_e and chi are grid quantities, not grain quantities
    mu_e    = AssignCellIds(mu_e,ids,mu_loc);
    chi     = AssignChiIds(chi,ids,chi_loc);

end


% -------------------------------------------------------------------------
% Pack output
% -------------------------------------------------------------------------
c_out       = PackC(c_grain,ny);
E           = PackCell(E,ny);
mu_e        = PackCell(mu_e,ny);
chi         = PackChi(chi,ny);

% -------------------------------------------------------------------------
% Raw thermodynamics for e and omega
% -------------------------------------------------------------------------
% c_grain is still unpacked as row-vector fields, which is what PhaseG uses.
% c_out is packed and used for Calc_e.
%
% Here omega is intentionally calculated with MODEL.pars, not pars_inter.
% Therefore W scaling, if used in LE, does not directly change the final
% thermodynamic driving force.
e_out       = Calc_e(MODEL.pars,c_out);
omg_raw     = CalcOmegaLocal(MODEL.pars,c_grain,e_out,mu_e,ny,nx,Ne,numel(c_out));
omg         = omg_raw;

% -------------------------------------------------------------------------
% Save LE state memory
% -------------------------------------------------------------------------
LE_state                 = struct();
LE_state.mode            = mode;
LE_state.active          = active;
LE_state.p_le            = p_le;
LE_state.p_th            = p_th;
LE_state.phase_index     = phase_index;
LE_state.phase_id        = phase_id;
LE_state.grain_to_phase  = grain_to_phase;
LE_state.n_active_map    = n_active_map;
LE_state.set_id_map      = set_id_map;
LE_state.p_tail          = p_tail;
LE_state.p_full          = p_full;
LE_state.p_on            = p_on;
LE_state.p_off           = p_off;
LE_state.Pmax            = Pmax;
LE_state.diag            = LE_diag;
if isfield(PARAM,'LE_c_write_pmin'),  LE_state.c_write_pmin  = PARAM.LE_c_write_pmin;  end
if isfield(PARAM,'LE_c_full_pmin'),   LE_state.c_full_pmin   = PARAM.LE_c_full_pmin;   end
if isfield(PARAM,'LE_c_write_dcmax'), LE_state.c_write_dcmax = PARAM.LE_c_write_dcmax; end

% -------------------------------------------------------------------------
% Update STATE
% -------------------------------------------------------------------------
STATE.c        = c_out;
STATE.e        = e_out;
STATE.E        = E;
STATE.mu_e     = mu_e;
STATE.chi      = chi;
STATE.omg      = omg;
STATE.omg_raw  = omg_raw;
STATE.LE_state = LE_state;
STATE.LE_diag  = LE_diag;
end




% ALL HELPER FUNCTIONS BELOW
% =========================================================================
% LE diagnostics
% =========================================================================
function diag = InitLEDiag(n_sets,n_nodes)
diag = struct();
diag.failed = false;
diag.n_sets = n_sets;
diag.n_nodes = n_nodes;
diag.n_failed_sets = 0;
diag.n_failed_nodes = 0;
diag.failed_set_ids = [];
diag.failed_node_sample = [];
diag.max_dc = 0;
diag.max_cchg = 0;
diag.max_iter = 0;
diag.min_alpha_stage = inf;
diag.messages = {};
end


function diag = InitLocalLEDiag(n_nodes)
diag = struct();
diag.converged = true;
diag.failed = false;
diag.n_nodes = n_nodes;
diag.alpha_stage = NaN;
diag.n_iter = 0;
diag.max_dc = 0;
diag.max_cchg = 0;
diag.message = '';
end


function diag = AccumulateLEDiag(diag,diag_loc,ids,iset)
if isfield(diag_loc,'max_dc') && isfinite(diag_loc.max_dc)
    diag.max_dc = max(diag.max_dc,diag_loc.max_dc);
end
if isfield(diag_loc,'max_cchg') && isfinite(diag_loc.max_cchg)
    diag.max_cchg = max(diag.max_cchg,diag_loc.max_cchg);
end
if isfield(diag_loc,'n_iter') && isfinite(diag_loc.n_iter)
    diag.max_iter = max(diag.max_iter,diag_loc.n_iter);
end
if isfield(diag_loc,'alpha_stage') && isfinite(diag_loc.alpha_stage)
    if isfinite(diag.min_alpha_stage)
        diag.min_alpha_stage = min(diag.min_alpha_stage,diag_loc.alpha_stage);
    else
        diag.min_alpha_stage = diag_loc.alpha_stage;
    end
end

if isfield(diag_loc,'failed') && diag_loc.failed
    diag.failed = true;
    diag.n_failed_sets = diag.n_failed_sets + 1;
    diag.n_failed_nodes = diag.n_failed_nodes + numel(ids);
    diag.failed_set_ids(end+1) = iset;

    n_keep = max(0,200-numel(diag.failed_node_sample));
    if n_keep > 0
        diag.failed_node_sample = [diag.failed_node_sample ids(1:min(n_keep,numel(ids)))]; %#ok<AGROW>
    end

    if isfield(diag_loc,'message') && ~isempty(diag_loc.message)
        diag.messages{end+1} = diag_loc.message; %#ok<AGROW>
    end
end

if isinf(diag.min_alpha_stage)
    diag.min_alpha_stage = NaN;
end
end


% =========================================================================
% Omega
% =========================================================================
function omg_col = CalcOmegaLocal(pars,c,e_col,mu_e,ny,nx,Ne,Np)
omg_col = zeros(ny,nx,Np);
for ip = 1:Np
    omg_col(:,:,ip) = reshape(PhaseG(pars{ip},c{ip}),ny,[]);
    for ie = 1:Ne
        omg_col(:,:,ip) = omg_col(:,:,ip) - e_col{ip}{ie}.*mu_e{ie};
    end
end
end


% =========================================================================
% Pack / unpack
% =========================================================================
function p = UnpackP(p)
[ny,nx,np] = size(p);
p          = reshape(p,1,ny*nx,np);
end


function c = UnpackC(c)
for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        c{ip}{ic} = reshape(c{ip}{ic},1,[]);
    end
end
end


function f = UnpackCell(f)
for i = 1:numel(f)
    f{i} = reshape(f{i},1,[]);
end
end


function ch = UnpackChi(ch)
for i = 1:size(ch,1)
    for j = 1:size(ch,2)
        ch{i,j} = reshape(ch{i,j},1,[]);
    end
end
end


function c = PackC(c,ny)
for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        c{ip}{ic} = reshape(c{ip}{ic},ny,[]);
    end
end
end


function f = PackCell(f,ny)
for i = 1:numel(f)
    f{i} = reshape(f{i},ny,[]);
end
end


function ch = PackChi(ch,ny)
for i = 1:size(ch,1)
    for j = 1:size(ch,2)
        ch{i,j} = reshape(ch{i,j},ny,[]);
    end
end
end


% =========================================================================
% WScale map
% =========================================================================
function pars_cur = Apply_WScale_Map(pars,w_scale_phase,ny,nx)
pars_cur = pars;
Np       = numel(pars);
for ip = 1:Np
    ws = w_scale_phase(:,:,ip);
    if isscalar(ws)
        pars_cur{ip}.w_scale = ws;
    else
        pars_cur{ip}.w_scale = reshape(ws,1,ny*nx);
    end
end
end



% =========================================================================
% Build thermodynamic phase information from grains
% =========================================================================
function [phase_id,grain_to_phase,pars_phase,p_phase] = BuildPhaseInfo(pars_grain,p_grain,phase_index)
%BUILDPHASEINFO Build phase-level p and grain-to-phase map.
%
% This does not collapse c.
% It only sums p over grains with the same thermodynamic phase ID.

Ngrain = size(p_grain,3);
N      = size(p_grain,2);

if numel(phase_index) ~= Ngrain
    error('BuildPhaseInfo: phase_index length does not match number of grains.')
end

phase_id = unique(phase_index,'stable');
Nphase   = numel(phase_id);

grain_to_phase = zeros(1,Ngrain);
pars_phase     = cell(1,Nphase);
p_phase        = zeros(1,N,Nphase);

for iph = 1:Nphase

    grains = find(phase_index == phase_id(iph));

    if isempty(grains)
        error('BuildPhaseInfo: phase %d has no grains.',phase_id(iph))
    end

    grain_to_phase(grains) = iph;

    %Use first grain as representative thermodynamic parameter.
    pars_phase{iph} = pars_grain{grains(1)};

    %Collapse p only.
    for ig = grains
        p_phase(:,:,iph) = p_phase(:,:,iph) + p_grain(:,:,ig);
    end
end

if any(grain_to_phase < 1)
    error('BuildPhaseInfo: some grains were not assigned.')
end

end



% =========================================================================
% Build local LE problem for one active-set group
% =========================================================================
function [pars_loc,c_loc,p_loc,E_loc,eta_loc] = BuildLocalLEBlock( ...
    pars_phase,c_grain,p_grain,e_grain,E,eta, ...
    phase_index,phase_id,ph_act,ids)

Nloc = numel(ids);
k    = numel(ph_act);
Ne   = numel(E);

pars_loc = cell(1,k);
c_loc    = cell(1,k);
p_loc    = zeros(1,Nloc,k);

%Full E at this local block
E_loc = SliceCellIds(E,ids);

%Original thermodynamic phase IDs that are active in this block
active_phase_id = phase_id(ph_act);

% ------------------------------------------------------------
% Build local p_loc and c_loc for active phases
% ------------------------------------------------------------
for ia = 1:k

    iph = ph_act(ia);
    pid = phase_id(iph);

    grains = find(phase_index == pid);
    ig0    = grains(1);

    pars_loc{ia} = pars_phase{iph};

    % If node-wise W scale is present, slice it to the current local ids.
    if isfield(pars_loc{ia},'w_scale') && ~isempty(pars_loc{ia}.w_scale)

        ws = reshape(pars_loc{ia}.w_scale,1,[]);

        if isscalar(ws)

            pars_loc{ia}.w_scale = ws;

        elseif numel(ws) == Nloc

            pars_loc{ia}.w_scale = ws;

        elseif numel(ws) >= max(ids)

            pars_loc{ia}.w_scale = ws(ids);

        else

            error('BuildLocalLEBlock: w_scale has incompatible size.')

        end
    end

    Nc = numel(c_grain{ig0});
    c_loc{ia} = cell(1,Nc);

    %Collapsed local p of this thermodynamic phase
    p_sum = zeros(1,Nloc);

    for ig = grains
        pg = reshape(p_grain(:,:,ig),1,[]);
        p_sum = p_sum + pg(ids);
    end

    p_loc(1,:,ia) = p_sum;

    %p-weighted c of this thermodynamic phase
    good = p_sum > eps;

    for ic = 1:Nc

        num = zeros(1,Nloc);

        for ig = grains
            pg = reshape(p_grain(:,:,ig),1,[]);
            cg = reshape(c_grain{ig}{ic},1,[]);
            num = num + pg(ids).*cg(ids);
        end

        cg0 = reshape(c_grain{ig0}{ic},1,[]);
        tmp = cg0(ids);

        tmp(good) = num(good)./p_sum(good);

        c_loc{ia}{ic} = tmp;

    end
end

% ------------------------------------------------------------
% Subtract inactive grain contributions from E_loc
% ------------------------------------------------------------
inactive_grains = find(~ismember(phase_index,active_phase_id));

Efix = zeros(Ne,Nloc);

for ig = inactive_grains

    pg = reshape(p_grain(:,:,ig),1,[]);
    pg = pg(ids);

    for ie = 1:Ne
        eg = reshape(e_grain{ig}{ie},1,[]);
        Efix(ie,:) = Efix(ie,:) + pg.*eg(ids);
    end
end

for ie = 1:Ne
    E_loc{ie} = E_loc{ie} - Efix(ie,:);
end

eta_loc = eta(ids);

end


% =========================================================================
% Assign local phase c back to grain-resolved storage
% =========================================================================
function c_grain = AssignLocalCToGrains(c_grain,p_grain,phase_index,phase_id,ph_act,ids,c_loc,PARAM,p_tail,p_full)
%ASSIGNLOCALCTOGRAINS Write local LE compositions back safely.
%
%The LE solve may include phases with very small p so that their chemical
%driving force is not lost. However, the composition of such low-p phases is
%weakly constrained and can jump. Therefore the write-back is p-weighted and
%optionally capped:
%
%   p <= LE_c_write_pmin : keep old c
%   p >= LE_c_full_pmin  : full write
%   between them         : smooth partial write
%
%This is the replacement for the unstable ghost return-map.

p_write = p_tail;
p_fullw = p_full;
dcmax   = inf;

if isfield(PARAM,'LE_c_write_pmin')
    p_write = PARAM.LE_c_write_pmin;
end
if isfield(PARAM,'LE_c_full_pmin')
    p_fullw = PARAM.LE_c_full_pmin;
end
if isfield(PARAM,'LE_c_write_dcmax')
    dcmax = PARAM.LE_c_write_dcmax;
end

p_write = max(p_write,0);
p_fullw = max(p_fullw,p_write+eps);

for ia = 1:numel(ph_act)

    iph = ph_act(ia);
    pid = phase_id(iph);

    grains = find(phase_index == pid);

    for ig = grains

        pg = reshape(p_grain(:,:,ig),1,[]);
        pp = pg(ids);

        w  = (pp - p_write)./(p_fullw - p_write);
        w  = min(max(w,0),1);
        w  = w.^2.*(3 - 2*w);

        take = w > 0;

        if ~any(take)
            continue
        end

        id_take = ids(take);
        ww      = w(take);

        for ic = 1:numel(c_grain{ig})

            A = reshape(c_grain{ig}{ic},1,[]);
            B = reshape(c_loc{ia}{ic},1,[]);

            dA = B(take) - A(id_take);

            if isfinite(dcmax) && dcmax > 0
                dA = min(max(dA,-dcmax),dcmax);
            end

            A(id_take) = A(id_take) + ww.*dA;
            c_grain{ig}{ic} = A;

        end
    end
end

end


function f = SliceCellIds(f,ids)
f2 = cell(size(f));
for i = 1:numel(f)
    x = reshape(f{i},1,[]);
    f2{i} = x(ids);
end
f = f2;
end


function f = AssignCellIds(f,ids,f_sub)
for i = 1:numel(f)
    x = reshape(f{i},1,[]);
    x(ids) = reshape(f_sub{i},1,[]);
    f{i} = x;
end
end


function ch = AssignChiIds(ch,ids,ch_sub)
Ne = size(ch,1);
for i = 1:Ne
    for j = 1:Ne
        x = reshape(ch{i,j},1,[]);
        x(ids) = reshape(ch_sub{i,j},1,[]);
        ch{i,j} = x;
    end
end
end



function [p_le,p_th,active] = BuildActiveSet_Tail(p,LE_state,Pmax,p_tail,p_full,p_on,p_off)

N  = size(p,2);
Np = size(p,3);

%Absolute smoothed phase amount used for active-set decision.
%This is not normalized.
p_abs = CalcThermoWeight_Tail(p,p_tail,p_full);
p2    = reshape(p_abs,N,Np);

%Normalized thermodynamic p only for diagnostics.
p_th  = p_abs./max(sum(p_abs,3),eps);
bad   = sum(p_abs,3) <= eps;

if any(bad(:))
    [~,idmax] = max(p,[],3);
    for ip = 1:Np
        tmp = p_th(:,:,ip);
        tmp(bad & idmax == ip) = 1;
        p_th(:,:,ip) = tmp;
    end
end

%Hysteresis memory
if isfield(LE_state,'active') && isequal(size(LE_state.active),[N,Np])
    active_old = LE_state.active;
else
    active_old = false(N,Np);
end

%Active decision based on absolute smoothed amount
active = (active_old & (p2 > p_off)) | (~active_old & (p2 > p_on));

%Make sure every grid has at least one active phase
p_phys = reshape(p,N,Np);

for i = 1:N

    if ~any(active(i,:))
        [~,im] = max(p_phys(i,:));
        active(i,im) = true;
    end

    %Limit maximum number of active phases
    if sum(active(i,:)) > Pmax
        score   = p2(i,:) + 0.5*p_on*active_old(i,:);
        score   = score + 1e-12*p_phys(i,:);
        [~,ord] = sort(score,'descend');
        active(i,:) = false;
        active(i,ord(1:Pmax)) = true;
    end
end

%Normalized active p for diagnostics only
p_le = p_abs.*reshape(active,1,N,Np);
p_le = p_le./max(sum(p_le,3),eps);
end


function p_abs = CalcThermoWeight_Tail(p,p_tail,p_full)
x = (p - p_tail)./(p_full - p_tail);
x = min(max(x,0),1);
%Absolute smooth turn-on. Do not normalize here.
p_abs = p.*x.^2.*(3 - 2*x);
end
