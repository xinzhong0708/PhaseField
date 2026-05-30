function [STATE,DIAG] = LE_Run_FromMu_AutoClosure_Simple(STATE,PARAM,MODEL)
%LE_RUN_FROMMU_AUTOCLOSURE_SIMPLE
% Clean mu-primary post-projection with automatic E-closure switching.
%
% Purpose
% -------
% This function is used after PF_Coupled_ACCH_LETangent.
%
% It keeps the smooth incoming chemical potential STATE.mu_e, then updates
% phase compositions c and omega.  The only special feature is the automatic
% choice of the output E closure:
%
%   Regular thermodynamics:
%       trust mu, solve fixed-mu GP stationarity, then set
%           E = sum(p*e(c)) + mu/eta.
%
%   Dominant spinodal / nearly singular thermodynamics:
%       trust the conserved CH predictor E, then correct c so that
%           E_in ~= sum(p*e(c)) + mu/eta.
%
% Why LE_Run is still needed before this function
% -----------------------------------------------
% LE_Run is E-primary: it gives the physical reference state and chi tangent
% for PF_Coupled_ACCH_LETangent.  This function is mu-primary: it is only a
% post-projection of the predictor.  The two calls should remain separate in
% the current time scheme.

% -------------------------------------------------------------------------
% 1. Unpack fields and collapse repeated grains
% -------------------------------------------------------------------------
pars        = MODEL.pars;
phase_index = MODEL.phase_index(:).';

p       = STATE.p;
c       = STATE.c;
E_in    = STATE.E;
mu_e    = STATE.mu_e;
chi     = STATE.chi;
eta_vec = PARAM.eta(:);

if isfield(STATE,'LE_state') && ~isempty(STATE.LE_state)
    LE_state = STATE.LE_state;
else
    LE_state = struct();
end

ny = size(p,1);

c    = Unpack_c(c);
p    = Unpack_p(p);
E_in = Unpack_E(E_in);
mu_e = Unpack_E(mu_e);
chi  = Unpack_Chi(chi);

[pars,c,p,grain_to_phase] = Collapse_LE_Phases(pars,c,p,phase_index);

Np = numel(c);
N  = numel(c{1}{1});
Ne = numel(mu_e);

opt = ReadOptions(PARAM);

% -------------------------------------------------------------------------
% 2. Build active thermodynamic phase weights
% -------------------------------------------------------------------------
[p_le,active,LE_state] = BuildActivePhases(p,LE_state,opt);
[active_sets,~,set_id] = unique(active,'rows');

if opt.use_WScale
    pars_inter = Apply_WScale_FromP(pars,p,0.85,0.99,1);
end

% -------------------------------------------------------------------------
% 3. Decide the E closure BEFORE nonlinear correction
% -------------------------------------------------------------------------
% seed_bad(i)        = true if node i has unsafe curvature in an active phase
% seed_bad_weight(i) = largest p_le of an unsafe phase at node i
seed_bad        = false(1,N);
seed_bad_weight = zeros(1,N);

for iset = 1:size(active_sets,1)

    ph_act = find(active_sets(iset,:));
    k      = numel(ph_act);
    mask   = (set_id == iset).';

    if opt.use_WScale && k > 1
        pars_cur = Slice_Pars_WScale(pars_inter,ph_act,mask);
    else
        pars_cur = pars(ph_act);
    end

    c_cur = Slice_c(c,ph_act,mask);
    p_cur = p_le(:,mask,ph_act);

    [bad_cur,bad_w_cur] = ClassifyBadNodes(pars_cur,p_cur,c_cur,opt);

    seed_bad(mask)        = bad_cur;
    seed_bad_weight(mask) = bad_w_cur;
end

old_mode = false;
if isfield(LE_state,'preserve_E_mode') && ~isempty(LE_state.preserve_E_mode)
    old_mode = logical(LE_state.preserve_E_mode);
end

preserve_E = SelectClosure(old_mode,seed_bad_weight,opt);
LE_state.preserve_E_mode = preserve_E;

% -------------------------------------------------------------------------
% 4. Correct compositions for each active phase set
% -------------------------------------------------------------------------
stat_node   = zeros(1,N);
newton_iter = zeros(1,N);
proj_iter   = zeros(1,N);
bad_node    = false(1,N);
bad_weight  = zeros(1,N);

for iset = 1:size(active_sets,1)

    ph_act = find(active_sets(iset,:));
    k      = numel(ph_act);
    mask   = (set_id == iset).';

    if opt.use_WScale && k > 1
        pars_cur = Slice_Pars_WScale(pars_inter,ph_act,mask);
    else
        pars_cur = pars(ph_act);
    end

    c_cur   = Slice_c(c,ph_act,mask);
    p_cur   = p_le(:,mask,ph_act);
    mu_cur  = Slice_E(mu_e,mask);
    E_cur   = Slice_E(E_in,mask);
    eta_cur = eta_vec(mask);

    local_opt                 = opt;
    local_opt.alpha           = opt.alpha_LE(k);
    local_opt.max_newton_iter = opt.max_newton_iter(k);

    [c_tmp,D] = LocalFromMu( ...
        pars_cur,p_cur,c_cur,mu_cur,E_cur,eta_cur,preserve_E,local_opt);

    % chi is intentionally not rebuilt here.  In the current time scheme,
    % the next step starts from LE_Run, which recomputes physical chi.
    chi_tmp = Slice_Chi_Local(chi,mask);

    [c,mu_e,chi] = Assign_LE_Back( ...
        c,mu_e,chi,c_tmp,mu_cur,chi_tmp,ph_act,mask);

    stat_node(mask)   = D.stat;
    newton_iter(mask) = D.newton_iter;
    proj_iter(mask)   = D.proj_iter;
    bad_node(mask)    = D.bad;
    bad_weight(mask)  = D.bad_weight;
end

% -------------------------------------------------------------------------
% 5. Whole-step E closure
% -------------------------------------------------------------------------
% Do not splice E node-by-node.
E_from_mu = Calc_E_FromMu(pars,c,p_le,mu_e,eta_vec);

if preserve_E
    E = E_in;        % spinodal: keep conservative CH predictor
else
    E = E_from_mu;   % regular GP: keep smooth mu closure
end

% -------------------------------------------------------------------------
% 6. Pack output and compute omega
% -------------------------------------------------------------------------
c_col = Pack_c(c,ny);
E     = Pack_E(E,ny);
mu_e  = Pack_E(mu_e,ny);
chi   = Pack_chi(chi,ny);
e_col = Calc_e(pars,c_col);

omg_col = zeros(ny,size(STATE.p,2),Np);
for ip = 1:Np
    omg_col(:,:,ip) = reshape(PhaseG(pars{ip},c{ip}),ny,[]);
    for ie = 1:Ne
        omg_col(:,:,ip) = omg_col(:,:,ip) - e_col{ip}{ie}.*mu_e{ie};
    end
end

% Copy collapsed thermodynamic phases back to grain-resolved variables
Ngrain = numel(grain_to_phase);
c_out  = cell(1,Ngrain);
omg    = zeros(ny,size(STATE.p,2),Ngrain);

for ig = 1:Ngrain
    iph         = grain_to_phase(ig);
    c_out{ig}   = c_col{iph};
    omg(:,:,ig) = omg_col(:,:,iph);
end

LE_state.phase_index    = phase_index;
LE_state.grain_to_phase = grain_to_phase;

STATE.c        = c_out;
STATE.e        = Calc_e(MODEL.pars,c_out);
STATE.E        = E;
STATE.mu_e     = mu_e;
STATE.chi      = chi;
STATE.omg      = omg;
STATE.LE_state = LE_state;

DIAG.preserve_E_mode          = preserve_E;
DIAG.max_stationarity         = max(stat_node);
DIAG.stat_node                = reshape(stat_node,ny,[]);
DIAG.newton_iter              = reshape(newton_iter,ny,[]);
DIAG.project_iter             = reshape(proj_iter,ny,[]);
DIAG.max_newton_iter          = max(newton_iter);
DIAG.max_project_iter         = max(proj_iter);
DIAG.seed_spinodal_fraction   = mean(seed_bad);
DIAG.spinodal_fraction        = mean(bad_node);
DIAG.bad_phase_weight         = reshape(bad_weight,ny,[]);
DIAG.max_E_projection_mismatch = MaxCellDiff(Pack_E(E_from_mu,ny),Pack_E(E_in,ny));
DIAG.active                   = active;
DIAG.p_le                     = p_le;

end


% =========================================================================
% Local correction
% =========================================================================

function [c,D] = LocalFromMu(pars,p,c,mu_e,E_in,eta,preserve_E,opt)
%LOCALFROMMU Apply either regular GP Newton or conserved-E correction.

N  = numel(mu_e{1});

D.stat        = zeros(1,N);
D.newton_iter = zeros(1,N);
D.proj_iter   = zeros(1,N);
D.bad         = false(1,N);
D.bad_weight  = zeros(1,N);

if ~HasCompositionDOF(pars,c)
    return
end

[bad,bad_weight] = ClassifyBadNodes(pars,p,c,opt);
D.bad        = bad;
D.bad_weight = bad_weight;

if preserve_E
    % Spinodal branch:
    % Do not solve fixed-mu stationarity.  Correct c only to reduce
    % E_in - [sum(p*e(c)) + mu/eta].
    live = true(1,N);
    for it = 1:opt.conserved_iter
        merit0 = EMerit(pars,p,c,mu_e,eta,E_in);
        live   = live & merit0 > opt.E_tol;
        if ~any(live)
            break
        end
        dc = EClosureStep(pars,p,c,mu_e,E_in,eta,live,opt);
        dc = CapStep(dc,live,opt.dc_cap_conserved);
        [c,accepted] = AcceptEStep(pars,p,c,mu_e,E_in,eta,dc,live,merit0,opt);
        D.proj_iter(live) = it;
        live = live & ~accepted;
    end
else
    % Regular GP branch:
    % Safe nodes converge by exact fixed-mu Newton.
    live = ~bad;
    for it = 1:opt.max_newton_iter
        if ~any(live)
            break
        end

        dc = GPNewtonStep(pars,c,mu_e,live);
        dc = CapStep(dc,live,opt.dc_cap_good);

        [c,accepted,step] = AcceptGPStep(pars,p,c,mu_e,dc,live,opt);
        D.newton_iter(live) = it;

        stat = Stationarity(pars,c,mu_e);
        done = live & accepted & (step < opt.c_tol | stat < opt.stat_tol);

        % Do not keep iterating nodes where line search failed.
        failed = live & ~accepted;
        live   = live & ~(done | failed);
    end

    % Small unsafe tails in an otherwise regular step get one or two bounded
    % E-guided mixed corrections.  This prevents tiny bad phase tails from
    % creating discontinuities without switching the whole step to conserved-E.
    live = bad;
    for it = 1:opt.max_proj_iter
        if ~any(live)
            break
        end
        merit0 = EMerit(pars,p,c,mu_e,eta,E_in);
        dc = MixedProjectionStep(pars,p,c,mu_e,E_in,eta,live,opt);
        dc = CapStep(dc,live,opt.dc_cap_bad);
        [c,accepted] = AcceptEStep(pars,p,c,mu_e,E_in,eta,dc,live,merit0,opt);
        D.proj_iter(live) = it;
        live = live & ~accepted;
    end
end

D.stat = Stationarity(pars,c,mu_e);

end


function dc = GPNewtonStep(pars,c,mu_e,mask)
% Exact fixed-mu Newton step:
%   H dc = J' mu - mu_c

Np     = numel(c);
Ne     = numel(mu_e);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);
dc     = ZeroStepLike(c);

for ip = 1:Np
    R = PhaseThermo(pars{ip},c{ip});
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    H    = SymPages(R.H_c);
    JT   = permute(R.Jac,[2 1 3]);
    Nc   = size(H,1);
    Rmu  = reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),Nc,N) ...
           - StackFields(R.mu_c,N);

    d = reshape(pagemldivide(H(:,:,mask),reshape(Rmu(:,mask),Nc,1,[])),Nc,[]);

    for ic = 1:Nc
        dc{ip}{ic}(mask) = d(ic,:);
    end
end

end


function dc = EClosureStep(pars,p,c,mu_e,E_in,eta,mask,opt)
% Conserved-E spinodal correction:
%   min ||G dc - R_E||^2 + lam_E ||dc||^2
%
% where
%   R_E = E_in - sum(p*e(c)) - mu/eta
%   G   = [p_1 J_1, p_2 J_2, ...]

Np = numel(c);
Ne = numel(mu_e);
N  = numel(mu_e{1});
Nb = nnz(mask);
dc = ZeroStepLike(c);

if Nb == 0
    return
end

mu_mat = StackFields(mu_e,N);
E_mat  = StackFields(E_in,N);
eta_v  = EtaVector(eta,N);

Rall     = cell(1,Np);
idx      = cell(1,Np);
Nc_total = 0;
Emix     = zeros(Ne,N);

for ip = 1:Np
    Rall{ip} = PhaseThermo(pars{ip},c{ip});
    p_ip     = reshape(p(:,:,ip),1,N);
    Emix     = Emix + StackFields(Rall{ip}.e,N).*p_ip;

    if isempty(Rall{ip}.Jac)
        idx{ip} = [];
        continue
    end

    Nc = size(Rall{ip}.Jac,2);
    idx{ip} = Nc_total + (1:Nc);
    Nc_total = Nc_total + Nc;
end

if Nc_total == 0
    return
end

R_E = E_mat - Emix - mu_mat./eta_v;
G   = zeros(Ne,Nc_total,Nb);

for ip = 1:Np
    id = idx{ip};
    if isempty(id)
        continue
    end

    J     = Rall{ip}.Jac(:,:,mask);
    p_ip  = reshape(p(:,:,ip),1,N);
    p_sel = reshape(p_ip(mask),1,1,Nb);

    G(:,id,:) = J.*p_sel;
end

GT = permute(G,[2 1 3]);
A  = pagemtimes(GT,G) + repmat(eye(Nc_total),1,1,Nb).*opt.lam_E;
b  = pagemtimes(GT,reshape(R_E(:,mask),Ne,1,Nb));
d  = reshape(pagemldivide(A,b),Nc_total,Nb);

for ip = 1:Np
    id = idx{ip};
    if isempty(id)
        continue
    end
    for ic = 1:numel(id)
        dc{ip}{ic}(mask) = d(id(ic),:);
    end
end

end


function dc = MixedProjectionStep(pars,p,c,mu_e,E_in,eta,mask,opt)
% E-guided fallback for small unsafe tails in regular GP mode:
%   min ||G dc - R_E||^2 + w_mu ||H dc - R_mu||^2 + lam ||dc||^2

Np = numel(c);
Ne = numel(mu_e);
N  = numel(mu_e{1});
Nb = nnz(mask);
dc = ZeroStepLike(c);

if Nb == 0
    return
end

mu_mat = StackFields(mu_e,N);
E_mat  = StackFields(E_in,N);
eta_v  = EtaVector(eta,N);

Rall     = cell(1,Np);
idx      = cell(1,Np);
Rmu      = cell(1,Np);
Nc_total = 0;
Emix     = zeros(Ne,N);

for ip = 1:Np
    Rall{ip} = PhaseThermo(pars{ip},c{ip});
    p_ip     = reshape(p(:,:,ip),1,N);
    Emix     = Emix + StackFields(Rall{ip}.e,N).*p_ip;

    if isempty(Rall{ip}.mu_c) || isempty(Rall{ip}.H_c) || isempty(Rall{ip}.Jac)
        idx{ip} = [];
        continue
    end

    Nc = size(Rall{ip}.H_c,1);
    idx{ip} = Nc_total + (1:Nc);
    Nc_total = Nc_total + Nc;

    JT = permute(Rall{ip}.Jac,[2 1 3]);
    Rmu{ip} = reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),Nc,N) ...
              - StackFields(Rall{ip}.mu_c,N);
end

if Nc_total == 0
    return
end

R_E = E_mat - Emix - mu_mat./eta_v;
G   = zeros(Ne,Nc_total,Nb);
T   = zeros(Nc_total,Nc_total,Nb);
r   = zeros(Nc_total,1,Nb);

for ip = 1:Np
    id = idx{ip};
    if isempty(id)
        continue
    end

    J     = Rall{ip}.Jac(:,:,mask);
    H     = SymPages(Rall{ip}.H_c(:,:,mask));
    HT    = permute(H,[2 1 3]);
    p_ip  = reshape(p(:,:,ip),1,N);
    p_sel = reshape(p_ip(mask),1,1,Nb);

    G(:,id,:)  = J.*p_sel;
    T(id,id,:) = pagemtimes(HT,H);
    r(id,1,:)  = pagemtimes(HT,reshape(Rmu{ip}(:,mask),numel(id),1,Nb));
end

GT = permute(G,[2 1 3]);
A  = pagemtimes(GT,G) + opt.wmu_proj*T + repmat(eye(Nc_total),1,1,Nb).*opt.lam_proj;
b  = pagemtimes(GT,reshape(R_E(:,mask),Ne,1,Nb)) + opt.wmu_proj*r;
d  = reshape(pagemldivide(A,b),Nc_total,Nb);

for ip = 1:Np
    id = idx{ip};
    if isempty(id)
        continue
    end
    for ic = 1:numel(id)
        dc{ip}{ic}(mask) = d(id(ic),:);
    end
end

end


% =========================================================================
% Step acceptance
% =========================================================================

function [c,accepted,step] = AcceptGPStep(pars,p,c,mu_e,dc,mask,opt)

N         = numel(mu_e{1});
F0        = GPObjective(pars,p,c,mu_e);
accepted  = false(1,N);
alpha_acc = zeros(1,N);
alpha     = opt.alpha*ones(1,N);
alpha(~mask) = 0;

for ils = 1:opt.max_ls
    test = mask & ~accepted & alpha >= opt.alpha_min;
    if ~any(test)
        break
    end

    c_try = AddStep(c,dc,alpha);
    F_try = GPObjective(pars,p,c_try,mu_e);

    take = test & isfinite(F_try) & ...
           (F_try <= F0 + opt.energy_tol.*max(1,abs(F0)));

    accepted(take)  = true;
    alpha_acc(take) = alpha(take);
    alpha(test & ~take) = 0.2*alpha(test & ~take);
end

c    = AddStep(c,dc,alpha_acc);
step = NodeMaxStep(dc).*alpha_acc;

end


function [c,accepted] = AcceptEStep(pars,p,c,mu_e,E_in,eta,dc,mask,merit0,opt)

N         = numel(mu_e{1});
accepted  = false(1,N);
alpha_acc = zeros(1,N);
alpha     = opt.alpha*ones(1,N);
alpha(~mask) = 0;

for ils = 1:opt.max_proj_ls
    test = mask & ~accepted & alpha >= opt.alpha_min;
    if ~any(test)
        break
    end

    c_try = AddStep(c,dc,alpha);
    merit = EMerit(pars,p,c_try,mu_e,eta,E_in);

    take = test & isfinite(merit) & ...
           (merit < merit0 - opt.proj_accept_abs.*max(1,merit0) | merit <= opt.E_tol);

    accepted(take)  = true;
    alpha_acc(take) = alpha(take);
    alpha(test & ~take) = 0.2*alpha(test & ~take);
end

c = AddStep(c,dc,alpha_acc);

end


% =========================================================================
% Thermodynamic measures
% =========================================================================

function [bad,bad_weight] = ClassifyBadNodes(pars,p,c,opt)

Np = numel(c);
N  = numel(c{1}{1});

bad        = false(1,N);
bad_weight = zeros(1,N);

for ip = 1:Np
    R = PhaseThermo(pars{ip},c{ip});
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    unsafe = ~HessianSafe(SymPages(R.H_c),opt.h_rel_tol,opt.cond_max);
    p_ip   = reshape(p(:,:,ip),1,N);

    bad        = bad | unsafe;
    bad_weight = max(bad_weight,p_ip.*unsafe);
end

end


function stat = Stationarity(pars,c,mu_e)
% max |mu_c - J' mu| over solution phases

Np     = numel(c);
Ne     = numel(mu_e);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);

stat = zeros(1,N);

for ip = 1:Np
    R = PhaseThermo(pars{ip},c{ip});
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    mu_c = StackFields(R.mu_c,N);
    JT   = permute(R.Jac,[2 1 3]);
    JTmu = reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),size(mu_c,1),N);

    stat = max(stat,max(abs(mu_c-JTmu),[],1));
end

end


function merit = EMerit(pars,p,c,mu_e,eta,E_in)
% || E_in - [sum(p*e(c)) + mu/eta] ||

N    = numel(mu_e{1});
E_mu = Calc_E_FromMu(pars,c,p,mu_e,eta);
R    = StackFields(E_in,N) - StackFields(E_mu,N);
merit = sqrt(sum(R.^2,1));

end


function F = GPObjective(pars,p,c,mu_e)
% sum p * [g(c) - mu'*e(c)]

Np     = numel(c);
Ne     = numel(mu_e);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);

F = zeros(1,N);
try
    for ip = 1:Np
        R    = PhaseThermo(pars{ip},c{ip});
        g    = reshape(R.g,1,N);
        e    = StackFields(R.e,N);
        p_ip = reshape(p(:,:,ip),1,N);
        F    = F + p_ip.*(g - sum(e.*mu_mat,1));
    end
catch
    F = inf(1,N);
end

end


function E = Calc_E_FromMu(pars,c,p,mu_e,eta)

Np    = numel(c);
Ne    = numel(mu_e);
N     = numel(mu_e{1});
eta_v = EtaVector(eta,N);

Emix = zeros(Ne,N);

for ip = 1:Np
    R    = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + StackFields(R.e,N).*p_ip;
end

E = cell(1,Ne);
for ie = 1:Ne
    E{ie} = Emix(ie,:) + reshape(mu_e{ie},1,N)./eta_v;
end

end


% =========================================================================
% Active set and options
% =========================================================================

function [p_le,active,LE_state] = BuildActivePhases(p,LE_state,opt)

N  = size(p,2);
Np = size(p,3);

p_th   = Calc_Thermo_p(p,opt.p_tail,opt.p_full);
pThMat = reshape(p_th,N,Np);

if isfield(LE_state,'active') && isequal(size(LE_state.active),[N,Np])
    active_old = LE_state.active;
else
    active_old = false(N,Np);
end

active = (active_old & pThMat > opt.p_off) | (~active_old & pThMat > opt.p_on);

repair = find(~any(active,2) | sum(active,2) > opt.Pmax).';
for i = repair
    if ~any(active(i,:))
        [~,idmax] = max(pThMat(i,:));
        active(i,idmax) = true;
    end
    if sum(active(i,:)) > opt.Pmax
        score = pThMat(i,:) + 0.5*opt.p_on*active_old(i,:);
        [~,ord] = sort(score,'descend');
        active(i,:) = false;
        active(i,ord(1:opt.Pmax)) = true;
    end
end

p_le = p_th.*reshape(active,1,N,Np);
p_le = p_le./max(sum(p_le,3),eps);

LE_state.active = active;
LE_state.p_th   = p_th;
LE_state.p_le   = p_le;

end


function preserve_E = SelectClosure(old_mode,bad_weight,opt)

switch opt.E_mode
    case 'from_mu'
        preserve_E = false;
    case 'conserved'
        preserve_E = true;
    otherwise
        if old_mode
            preserve_E = any(bad_weight >= opt.p_exit);
        else
            preserve_E = any(bad_weight >= opt.p_enter);
        end
end

end


function opt = ReadOptions(PARAM)

opt.Pmax       = GetField(PARAM,'frommu_Pmax',4);
opt.p_tail     = GetField(PARAM,'frommu_p_tail',1e-3);
opt.p_full     = GetField(PARAM,'frommu_p_full',5e-2);
opt.p_on       = GetField(PARAM,'frommu_p_on',5e-3);
opt.p_off      = GetField(PARAM,'frommu_p_off',2e-3);
opt.alpha_LE   = GetField(PARAM,'frommu_alpha_LE',[0.8 0.5 0.4 0.3]);

opt.E_mode     = lower(GetField(PARAM,'frommu_E_mode','auto'));
opt.p_enter    = GetField(PARAM,'frommu_p_enter',0.80);
opt.p_exit     = GetField(PARAM,'frommu_p_exit',0.60);

opt.max_newton_iter = GetField(PARAM,'frommu_max_newton_iter',100);
if numel(opt.max_newton_iter) == 1
    opt.max_newton_iter = repmat(opt.max_newton_iter,1,4);
end
opt.c_tol      = GetField(PARAM,'frommu_c_tol',1e-7);
opt.stat_tol   = GetField(PARAM,'frommu_stat_tol',1e-7);
opt.max_ls     = GetField(PARAM,'frommu_max_ls',8);
opt.energy_tol = GetField(PARAM,'frommu_energy_tol',1e-10);
opt.alpha_min  = GetField(PARAM,'frommu_alpha_min',1e-8);
opt.dc_cap_good = GetField(PARAM,'frommu_dc_cap_good',5e-3);

opt.conserved_iter       = GetField(PARAM,'frommu_conserved_iter',1);
opt.dc_cap_conserved     = GetField(PARAM,'frommu_dc_cap_conserved',5e-4);
opt.lam_E                = GetField(PARAM,'frommu_lam_E',1e-12);

opt.max_proj_iter   = GetField(PARAM,'frommu_max_proj_iter',2);
opt.max_proj_ls     = GetField(PARAM,'frommu_max_proj_ls',3);
opt.dc_cap_bad      = GetField(PARAM,'frommu_dc_cap_bad',5e-4);
opt.E_tol           = GetField(PARAM,'frommu_E_tol',1e-8);
opt.proj_accept_abs = GetField(PARAM,'frommu_proj_accept_abs',1e-12);
opt.lam_proj        = GetField(PARAM,'frommu_lam_proj',1e-8);
opt.wmu_proj        = GetField(PARAM,'frommu_wmu_proj',1);

opt.h_rel_tol = GetField(PARAM,'frommu_h_rel_tol',1e-6);
opt.cond_max  = GetField(PARAM,'frommu_cond_max',1e8);

opt.use_WScale = isfield(PARAM,'use_WScale') && PARAM.use_WScale;

end


% =========================================================================
% Small utilities
% =========================================================================

function value = GetField(S,name,default)
if isfield(S,name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default;
end
end


function tf = HasCompositionDOF(pars,c)
tf = false;
for ip = 1:numel(c)
    R = PhaseThermo(pars{ip},c{ip});
    if ~(isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac))
        tf = true;
        return
    end
end
end


function safe = HessianSafe(H,h_rel_tol,cond_max)

N = size(H,3);
safe = false(1,N);

for i = 1:N
    ev = eig(H(:,:,i));
    sc = max(1,max(abs(ev)));
    safe(i) = min(ev) > h_rel_tol*sc && ...
              max(abs(ev))/max(min(abs(ev)),eps*sc) < cond_max;
end

end


function H = SymPages(H)
H = 0.5*(H + permute(H,[2 1 3]));
end


function A = StackFields(fields,N)
A = zeros(numel(fields),N);
for i = 1:numel(fields)
    A(i,:) = reshape(fields{i},1,N);
end
end


function eta_v = EtaVector(eta,N)
if isscalar(eta)
    eta_v = eta*ones(1,N);
else
    eta_v = reshape(eta,1,N);
end
end


function dc = ZeroStepLike(c)
dc = c;
for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        dc{ip}{ic} = zeros(size(c{ip}{ic}));
    end
end
end


function dc = CapStep(dc,mask,cap)

if ~any(mask)
    return
end

amp = NodeMaxStep(dc);
fac = ones(size(mask));
over = mask & amp > cap;
fac(over) = cap./amp(over);

for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        dc{ip}{ic} = dc{ip}{ic}.*fac;
    end
end

end


function amp = NodeMaxStep(dc)

amp = zeros(size(dc{1}{1}));
for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        amp = max(amp,abs(dc{ip}{ic}));
    end
end

end


function c_new = AddStep(c,dc,alpha_node)

c_new = c;
for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        c_new{ip}{ic} = c{ip}{ic} + alpha_node.*dc{ip}{ic};
    end
end

end


function out = Slice_Chi_Local(chi,mask)

Ne = size(chi,1);
out = cell(Ne,Ne);
for i = 1:Ne
    for j = 1:Ne
        x = reshape(chi{i,j},1,[]);
        out{i,j} = x(mask);
    end
end

end


function d = MaxCellDiff(A,B)

d = 0;
for i = 1:numel(A)
    d = max(d,max(abs(A{i}(:)-B{i}(:))));
end

end


function pars_cur = Slice_Pars_WScale(pars_inter,ph_act,mask)

pars_cur = pars_inter(ph_act);
for ia = 1:numel(pars_cur)
    if isfield(pars_cur{ia},'w_scale') && ~isempty(pars_cur{ia}.w_scale)
        pars_cur{ia}.w_scale = pars_cur{ia}.w_scale(mask);
    end
end

end


function p_th = Calc_Thermo_p(p,p_tail,p_full)

x     = (p-p_tail)./(p_full-p_tail);
x     = min(max(x,0),1);
w     = p.*x.^2.*(3-2*x);
wsum  = sum(w,3);
p_th  = w./max(wsum,eps);

bad = wsum <= eps;
if any(bad(:))
    [~,idmax] = max(p,[],3);
    for ip = 1:size(p,3)
        tmp = p_th(:,:,ip);
        tmp(bad & idmax == ip) = 1;
        p_th(:,:,ip) = tmp;
    end
end

end


function [pars_c,c_c,p_c,grain_to_phase] = Collapse_LE_Phases(pars,c,p,phase_index)

phase_id       = unique(phase_index,'stable');
Ngrain         = numel(c);
Nphase         = numel(phase_id);
N              = size(p,2);
grain_to_phase = zeros(1,Ngrain);

for iph = 1:Nphase
    grain_to_phase(phase_index == phase_id(iph)) = iph;
end

pars_c = cell(1,Nphase);
for iph = 1:Nphase
    ig = find(grain_to_phase == iph,1,'first');
    pars_c{iph} = pars{ig};
end

p_c = zeros(1,N,Nphase);
for ig = 1:Ngrain
    iph = grain_to_phase(ig);
    p_c(:,:,iph) = p_c(:,:,iph) + p(:,:,ig);
end

c_c = cell(1,Nphase);
for iph = 1:Nphase

    grains = find(grain_to_phase == iph);
    ig0    = grains(1);
    Nc     = numel(c{ig0});
    den    = reshape(sum(p(:,:,grains),3),size(c{ig0}{1}));
    good   = den > eps;

    c_c{iph} = cell(1,Nc);

    for ic = 1:Nc
        num = zeros(size(c{ig0}{ic}));
        for ig = grains
            num = num + reshape(p(:,:,ig),size(num)).*c{ig}{ic};
        end

        tmp = c{ig0}{ic};
        tmp(good) = num(good)./den(good);
        c_c{iph}{ic} = tmp;
    end
end

end
