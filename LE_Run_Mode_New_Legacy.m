function [STATE] = LE_Run_Mode_New(STATE,PARAM,MODEL)
%This function perform local equilibrium
%Two modes are available
%1) 'LE': Classical KKS type fixed E mode with fixed p and E
%2) 'GP': Grand potential based local equilibrium with fixed mu_e

%Choose the mode
mode       =  'LE';

if isfield(PARAM,'LE_mode')
    mode = PARAM.LE_mode;
end

%p mode on and off
p_tail     =  1e-3;
p_full     =  5e-2;
p_on       =  1e-2;
p_off      =  8e-3;

%Maximal number of active phases for LE/GP
Pmax       =  4;

%Local equilibrium parameters
alpha_LE   = [0.7 0.4 0.3 0.2];
iter_LE    = [100 100 100 100];
iter_GP    = [100 100 100 100];

%Optional active-set controls
if isfield(PARAM,'LE_p_tail')
    p_tail = PARAM.LE_p_tail;
end
if isfield(PARAM,'LE_p_full')
    p_full = PARAM.LE_p_full;
end
if isfield(PARAM,'LE_p_on')
    p_on = PARAM.LE_p_on;
end
if isfield(PARAM,'LE_p_off')
    p_off = PARAM.LE_p_off;
end
if isfield(PARAM,'LE_Pmax')
    Pmax = PARAM.LE_Pmax;
end
if isfield(PARAM,'LE_alpha_LE')
    alpha_LE = PARAM.LE_alpha_LE;
end
if isfield(PARAM,'LE_iter_LE')
    iter_LE = PARAM.LE_iter_LE;
end
if isfield(PARAM,'LE_iter_GP')
    iter_GP = PARAM.LE_iter_GP;
end

%Unpack fields and prepare the calculation
ny         =  size(STATE.p,1);
nx         =  size(STATE.p,2);
N          =  ny*nx;
p          =  UnpackP(STATE.p);
c          =  UnpackC(STATE.c);
E          =  UnpackCell(STATE.E);
mu_e       =  UnpackCell(STATE.mu_e);
chi        =  UnpackChi(STATE.chi);
eta        =  reshape(PARAM.eta,1,[]);
pars       =  MODEL.pars;

%Keep grain-resolved fields before collapse
p_grain    =  p;
c_grain    =  c;

%Grain-resolved active c-update mask.
%This is the only grain-region mask used for writing c back.
%Tiny numerical p tails below p_tail cannot receive updated c.
grain_mask =  CalcThermoWeight_Tail(p_grain,p_tail,p_full) > 0;

%Phase index of the grains
phase_index=  MODEL.phase_index(:).';

%Collapse repeated grains into thermodynamic phases
[pars,c,p,grain_to_phase] = CollapsePhases(pars,c,p,phase_index);

%Phase number after collapse
Np         =  size(p,3);
Ne         =  numel(E);
Pmax       =  min(Pmax,Np);

%Damping of the interface excess energy
pars_orig      = pars;

if isfield(PARAM,'w_scale_phase') && ~isempty(PARAM.w_scale_phase)
    pars_inter = Apply_WScale_Map(pars,PARAM.w_scale_phase,ny,nx);
else
    pars_inter = pars;
end

%Find the active phase with smoothed thermodynamic p and hysteresis
if isfield(STATE,'LE_state') && ~isempty(STATE.LE_state)
    LE_state_old = STATE.LE_state;
else
    LE_state_old = struct();
end

[p_le,p_th,active] = BuildActiveSet_Tail(p,LE_state_old,Pmax,p_tail,p_full,p_on,p_off);
n_active_map       = reshape(sum(active,2),ny,nx);

%Group grids by active phase set, this will be doing for LE or GP
[active_sets,~,set_id] = unique(active,'rows');
set_id_map             = reshape(set_id,ny,nx);

%Elemental composition of all collapsed thermodynamic phases
%This is used only for inactive mass-balance subtraction.
e_all = Calc_e(pars_orig,c);

%MAIN PART OF THE LE
for iset = 1:size(active_sets,1)

    %Active phase
    ph_act       = find(active_sets(iset,:));
    k            = numel(ph_act);
    mask         = (set_id == iset).';

    %Slice active local problem
    if isfield(PARAM,'use_WScale') && PARAM.use_WScale
        pars_loc = SliceParsWScale(pars_inter,ph_act,mask);
    else
        pars_loc = pars(ph_act);
    end

    c_loc        = SliceC(c,ph_act,mask);
    p_loc        = SliceP(p,mask,ph_act);
    eta_loc      = eta(mask);

    %Choose damping based on number of active phases
    kk           = min(k,numel(alpha_LE));
    aa           = alpha_LE(kk);

    %MODE 1: KKS TYPE LOCAL EQUILIBRIUM WITH E AND P FIXED
    if strcmpi(mode,'LE')

        mm                     =  iter_LE(kk);
        [~,~,E_loc]            =  SliceActiveMassBalance_FromE(c,p,E,e_all,ph_act,mask);
        [c_loc,mu_loc,chi_loc] =  LE_Calculator(pars_loc,p_loc,c_loc,E_loc,eta_loc,[aa,mm]);

    elseif strcmpi(mode,'GP')

        [~,~,E_loc]            =  SliceActiveMassBalance_FromE(c,p,E,e_all,ph_act,mask);

        %For one phase, use LE because E is the conserved variable.
        %For multiphase regions, use mass-consistent GP as smooth closure.
        if k == 1
            mm                     =  iter_LE(kk);
            [c_loc,mu_loc,chi_loc] =  LE_Calculator(pars_loc,p_loc,c_loc,E_loc,eta_loc,[aa,mm]);
        else
            mm                     =  iter_GP(kk);
            mu_loc                 =  SliceCell(mu_e,mask);
            [c_loc,chi_loc]        =  GP_Calculator(pars_loc,p_loc,c_loc,mu_loc,E_loc,eta_loc,aa,mm);
        end

    else

        error('LE_Run_Mode_New: unknown LE_mode "%s".',mode)

    end

    %Assign back to collapsed thermodynamic phases
    c    = AssignC(c,ph_act,mask,c_loc);
    mu_e = AssignCell(mu_e,mask,mu_loc);
    chi  = AssignChi(chi,mask,chi_loc);

    %Assign back to grain-resolved c only where:
    % 1) the pixel was actually sent to this LE/GP local solve, and
    % 2) the individual grain is active according to grain_mask.
    c_grain = AssignC_ActiveMask_ToGrains(c_grain,grain_mask,grain_to_phase,ph_act,mask,c_loc);

end

% -------------------------------------------------------------------------
% Pack output
% -------------------------------------------------------------------------
c_out       = PackC(c_grain,ny);
E           = PackCell(E,ny);
mu_e        = PackCell(mu_e,ny);
chi         = PackChi(chi,ny);

%Grain-resolved e and omega
%e_out uses packed grain-resolved c_out.
%omg uses unpacked c_grain so PhaseG receives row-vector fields.
e_out       = Calc_e(MODEL.pars,c_out);
omg_raw     = CalcOmegaLocal(MODEL.pars,c_grain,e_out,mu_e,ny,nx,Ne,numel(c_out));
omg         = omg_raw;

%Save active-set memory
LE_state                 =  struct();
LE_state.active          =  active;
LE_state.p_le            =  p_le;
LE_state.p_th            =  p_th;
LE_state.phase_index     =  phase_index;
LE_state.grain_to_phase  =  grain_to_phase;
LE_state.n_active_map    =  n_active_map;
LE_state.set_id_map      =  set_id_map;
LE_state.p_tail          =  p_tail;
LE_state.p_full          =  p_full;
LE_state.p_on            =  p_on;
LE_state.p_off           =  p_off;
LE_state.Pmax            =  Pmax;

%Regularize chi
if isfield(PARAM,'regularize_chi_all') && PARAM.regularize_chi_all == 1
    if isfield(PARAM,'chi_reg_floor')
        chi_floor = PARAM.chi_reg_floor;
    else
        chi_floor = 1e-8;
    end
    chi = RegularizeChi_OnMask(chi,true(ny,nx),chi_floor);
end

if isfield(PARAM,'smooth_chi_all') && PARAM.smooth_chi_all == 1
    if isfield(PARAM,'smooth_chi_iter')
        nsmooth = PARAM.smooth_chi_iter;
    else
        nsmooth = 1;
    end
    chi = SmoothChi_Local(chi,nsmooth);
    if isfield(PARAM,'chi_reg_floor')
        chi = RegularizeChi_OnMask(chi,true(ny,nx),PARAM.chi_reg_floor);
    end
end

%Update STATE
STATE.c        = c_out;
STATE.e        = e_out;
STATE.E        = E;
STATE.mu_e     = mu_e;
STATE.chi      = chi;
STATE.omg      = omg;
STATE.omg_raw  = omg_raw;
STATE.LE_state = LE_state;

end


% ALL HELPER FUNCTIONS BELOW
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
% WScale slicing
% =========================================================================
function pars_cur = SliceParsWScale(pars_inter,ph_act,mask)
%SLICEPARSWSCALE Extract nodewise w_scale for one active-set block.

pars_cur = pars_inter(ph_act);

for ia = 1:numel(pars_cur)

    if isfield(pars_cur{ia},'w_scale') && ~isempty(pars_cur{ia}.w_scale)

        ws = reshape(pars_cur{ia}.w_scale,1,[]);

        if isscalar(ws)
            pars_cur{ia}.w_scale = ws;
        elseif numel(ws) == numel(mask)
            pars_cur{ia}.w_scale = ws(mask);
        elseif numel(ws) == nnz(mask)
            pars_cur{ia}.w_scale = ws;
        end
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
% Collapse grains
% =========================================================================
function [pars_c,c_c,p_c,grain_to_phase] = CollapsePhases(pars,c,p,phase_index)

phase_id        = unique(phase_index,'stable');
Ngrain          = numel(c);
Nphase          = numel(phase_id);
N               = size(p,2);
grain_to_phase  = zeros(1,Ngrain);

for iph = 1:Nphase
    grain_to_phase(phase_index == phase_id(iph)) = iph;
end

pars_c = cell(1,Nphase);
p_c    = zeros(1,N,Nphase);
c_c    = cell(1,Nphase);

for iph = 1:Nphase

    grains       = find(grain_to_phase == iph);
    ig0          = grains(1);
    pars_c{iph}  = pars{ig0};

    for ig = grains
        p_c(:,:,iph) = p_c(:,:,iph) + p(:,:,ig);
    end

    Nc       = numel(c{ig0});
    den      = reshape(p_c(:,:,iph),1,N);
    good     = den > eps;
    c_c{iph} = cell(1,Nc);

    for ic = 1:Nc

        num = zeros(1,N);

        for ig = grains
            num = num + reshape(p(:,:,ig),1,N).*reshape(c{ig}{ic},1,N);
        end

        tmp = reshape(c{ig0}{ic},1,N);
        tmp(good) = num(good)./den(good);

        c_c{iph}{ic} = tmp;

    end
end

end


% =========================================================================
% Slice / assign
% =========================================================================
function [c_blk,p_blk,E_blk] = SliceActiveMassBalance_FromE(c,p,E,e_all,ph_act,mask)

Np    = numel(c);
Ne    = numel(E);
Nloc  = nnz(mask);
c_blk = SliceC(c,ph_act,mask);
p_blk = SliceP(p,mask,ph_act);
E_blk = SliceCell(E,mask);

inactive = setdiff(1:Np,ph_act);
Efix     = zeros(Ne,Nloc);

for ii = 1:numel(inactive)

    ip   = inactive(ii);
    p_ip = reshape(p(:,mask,ip),1,Nloc);
    e_ip = zeros(Ne,Nloc);

    for ie = 1:Ne
        tmp        = reshape(e_all{ip}{ie},1,[]);
        e_ip(ie,:) = tmp(mask);
    end

    Efix = Efix + e_ip.*p_ip;

end

for ie = 1:Ne
    E_blk{ie} = E_blk{ie} - Efix(ie,:);
end

end


function c2 = SliceC(c,ph,mask)

c2 = cell(1,numel(ph));

for a = 1:numel(ph)

    ip    = ph(a);
    c2{a} = cell(size(c{ip}));

    for ic = 1:numel(c{ip})
        x          = reshape(c{ip}{ic},1,[]);
        c2{a}{ic} = x(mask);
    end
end

end


function p2 = SliceP(p,mask,ph)

p2 = p(:,mask,ph);
p2 = reshape(p2,1,nnz(mask),numel(ph));

end


function f2 = SliceCell(f,mask)

f2 = cell(size(f));

for i = 1:numel(f)
    x      = reshape(f{i},1,[]);
    f2{i} = x(mask);
end

end


function c = AssignC(c,ph,mask,c_sub)

for a = 1:numel(ph)

    ip = ph(a);

    for ic = 1:numel(c{ip})
        x          = reshape(c{ip}{ic},1,[]);
        x(mask)    = reshape(c_sub{a}{ic},1,[]);
        c{ip}{ic}  = x;
    end
end

end


function c_grain = AssignC_ActiveMask_ToGrains(c_grain,grain_mask,grain_to_phase,ph_act,mask,c_loc)

ids = find(mask);

if isempty(ids)
    return
end

for ia = 1:numel(ph_act)

    iph    = ph_act(ia);
    grains = find(grain_to_phase == iph);

    if isempty(grains)
        continue
    end

    for ig = grains

        gm   = reshape(grain_mask(:,:,ig),1,[]);
        take = gm(ids);

        if ~any(take)
            continue
        end

        id_take = ids(take);

        for ic = 1:numel(c_grain{ig})

            A = reshape(c_grain{ig}{ic},1,[]);
            B = reshape(c_loc{ia}{ic},1,[]);

            A(id_take) = B(take);

            c_grain{ig}{ic} = A;

        end
    end
end

end


function f = AssignCell(f,mask,f_sub)

for i = 1:numel(f)
    x       = reshape(f{i},1,[]);
    x(mask) = reshape(f_sub{i},1,[]);
    f{i}    = x;
end

end


function ch = AssignChi(ch,mask,ch_sub)

Ne = size(ch,1);

for i = 1:Ne
    for j = 1:Ne
        x       = reshape(ch{i,j},1,[]);
        x(mask) = reshape(ch_sub{i,j},1,[]);
        ch{i,j} = x;
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
% Chi utilities
% =========================================================================
function chi = SmoothChi_Local(chi,nsmooth)

Ne  = size(chi,1);
ker = [1 2 1; 2 4 2; 1 2 1];
ker = ker/sum(ker(:));

for it = 1:nsmooth
    for i = 1:Ne
        for j = 1:Ne
            chi{i,j} = conv2(chi{i,j},ker,'same');
        end
    end
end

end


function chi = RegularizeChi_OnMask(chi,mask,chi_floor)

Ne  = size(chi,1);
ids = find(mask(:));

for kk = 1:numel(ids)

    id = ids(kk);
    C  = zeros(Ne,Ne);

    for i = 1:Ne
        for j = 1:Ne
            tmp    = chi{i,j};
            C(i,j) = tmp(id);
        end
    end

    C     = 0.5*(C+C.');
    [V,D] = eig(C);
    lam   = diag(D);
    scale = max(1,max(abs(lam)));
    lam   = max(lam,chi_floor*scale);
    Creg  = V*diag(lam)*V.';
    Creg  = 0.5*(Creg+Creg.');

    for i = 1:Ne
        for j = 1:Ne
            tmp      = chi{i,j};
            tmp(id)  = Creg(i,j);
            chi{i,j} = tmp;
        end
    end

end

end


% =========================================================================
% Active set
% =========================================================================
function [p_le,p_th,active] = BuildActiveSet_Tail(p,LE_state,Pmax,p_tail,p_full,p_on,p_off)

N  = size(p,2);
Np = size(p,3);

%Absolute smoothed phase amount used for active-set decision.
%This is not normalized.
p_abs = CalcThermoWeight_Tail(p,p_tail,p_full);
p2    = reshape(p_abs,N,Np);

%Normalized thermodynamic p only for diagnostics.
p_th = p_abs./max(sum(p_abs,3),eps);

bad = sum(p_abs,3) <= eps;

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

        score = p2(i,:) + 0.5*p_on*active_old(i,:);
        score = score + 1e-12*p_phys(i,:);

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
