function [STATE,DIAG] = LE_Run_FromMu_AutoClosure_Lean(STATE,PARAM,MODEL)
%LE_RUN_FROMMU_AUTOCLOSURE_LEAN  Mu-primary post-projection with automatic E closure.
%
% Intended use:
%   STATE_OLD   = LE_Run(STATE_OLD,PARAM,MODEL);
%   STATE_RAW   = PF_Coupled_ACCH_LETangent(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);
%   [STATE_TRIAL,DIAG] = LE_Run_FromMu_AutoClosure_Fast(STATE_RAW,PARAM,MODEL);
%
% Design:
%   1) Convex, well-conditioned nodes:
%        solve the fixed-mu local equilibrium by Newton iteration.
%   2) Non-convex / ill-conditioned nodes:
%        apply only a few bounded E-guided projection sweeps.
%        These nodes are not asked to converge to fixed-mu stationarity.
%   3) E closure:
%        normal GP regime     -> E = E_from_mu for smooth mu/omega;
%        spinodal/singular    -> E = incoming conserved CH field.
%
% Controls can optionally be stored in PARAM:
%   frommu_E_mode          = 'auto' (default), 'from_mu', or 'conserved'
%   frommu_p_enter         = 0.80     % enter conserved-E closure
%   frommu_p_exit          = 0.60     % leave conserved-E closure
%   frommu_max_proj_iter   = 3        % fallback projection cap outside conserved-E mode
%   frommu_conserved_iter  = 1        % dominant spinodal: E-only correction sweeps
%   frommu_lam_E           = 1e-12    % regularization of E-only correction
%   frommu_h_rel_tol       = 1e-6     % near-singular Hessian cutoff
%   frommu_compute_chi     = false    % keep input chi; next LE_Run overwrites it
%   frommu_wmu_proj        = 1
%   frommu_lam_proj        = 1e-8
%   frommu_dc_cap_bad      = 5e-4
%
% Note:
%   This function is a post-projection. In the spinodal regime the incoming
%   E field is preserved; the returned mu/c are therefore a bounded smooth
%   representation, not an exact solution of an invertible E(mu) map.

opt = Read_FromMu_Options(PARAM);

% Unpack node fields
pars        = MODEL.pars;
phase_index = MODEL.phase_index(:).';
ny          = size(STATE.p,1);

p           = Unpack_p(STATE.p);
c           = Unpack_c(STATE.c);
E_in        = Unpack_E(STATE.E);
mu_e        = Unpack_E(STATE.mu_e);
chi         = Unpack_Chi(STATE.chi);
eta_vec     = reshape(PARAM.eta,1,[]);

% Collapse grains to thermodynamic phases
[pars,c,p,grain_to_phase] = Collapse_LE_Phases(pars,c,p,phase_index);

% Number of phase and grid
Np          = numel(c);
N           = numel(c{1}{1});

if isfield(STATE,'LE_state') && ~isempty(STATE.LE_state)
    LE_state = STATE.LE_state;
else
    LE_state = struct();
end

% ---------- Active phases and smooth thermodynamic weights ----------
[p_le,active,LE_state] = Build_Active_Phases(p,LE_state,opt);
[active_sets,~,set_id] = unique(active,'rows');

stat_node       = zeros(1,N);
newton_iter     = zeros(1,N);
proj_iter       = zeros(1,N);
bad_node        = false(1,N);
bad_weight      = zeros(1,N);
proj_stalled    = false(1,N);
newton_stalled  = false(1,N);

if opt.use_WScale
    pars_inter = Apply_WScale_FromP(pars,p,0.85,0.99,1);
end

% ---------- Decide the closure before nonlinear correction ----------
% This is deliberately based on the incoming predictor state.  Inner
% projection iterations must not change the global closure mode.
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
    [seed_bad(mask),seed_bad_weight(mask)] = Classify_Bad_Nodes(pars_cur,p_cur,c_cur,opt);
end

old_mode = false;
if isfield(LE_state,'preserve_E_mode') && ~isempty(LE_state.preserve_E_mode)
    old_mode = logical(LE_state.preserve_E_mode);
end
preserve_E = Select_E_Closure(old_mode,seed_bad_weight,opt);
opt.preserve_E = preserve_E;

% ---------- Solve each active-phase combination ----------
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
    E_cur   = Slice_E(E_in,mask);
    mu_cur  = Slice_E(mu_e,mask);
    eta_cur = eta_vec(mask);

    local_opt = opt;
    local_opt.alpha           = opt.alpha_LE(k);
    local_opt.max_newton_iter = opt.max_newton_iter(k);

    [c_tmp,chi_tmp,D] = Solve_Local_FromMu(pars_cur,p_cur,c_cur,mu_cur,eta_cur,E_cur,seed_bad(mask),seed_bad_weight(mask),local_opt);

    if ~opt.compute_chi
        chi_tmp = Slice_Chi_Local(chi,mask);
    end

    [c,mu_e,chi] = Assign_LE_Back(c,mu_e,chi,c_tmp,mu_cur,chi_tmp,ph_act,mask);

    stat_node(mask)      = D.stat;
    newton_iter(mask)    = D.newton_iter;
    proj_iter(mask)      = D.proj_iter;
    bad_node(mask)       = D.bad;
    bad_weight(mask)     = D.bad_weight;
    proj_stalled(mask)   = D.proj_stalled;
    newton_stalled(mask) = D.newton_stalled;
end

% ---------- Whole-step E closure; never splice E point-by-point ----------
E_from_mu = Calc_E_FromMu(pars,c,p_le,mu_e,eta_vec);
E_mu_pack = Pack_E(E_from_mu,ny);
E_in_pack = Pack_E(E_in,ny);

DIAG.max_E_projection_mismatch = Max_Cell_Diff_Local(E_mu_pack,E_in_pack);

if preserve_E
    E = E_in;       % Singular/spinodal map: retain conserved PDE result.
else
    E = E_from_mu;  % Regular GP map: retain smooth constitutive closure.
end

LE_state.preserve_E_mode = preserve_E;

% ---------- Pack and compute outputs ----------
c_col = Pack_c(c,ny);
E     = Pack_E(E,ny);
mu_e  = Pack_E(mu_e,ny);
chi   = Pack_chi(chi,ny);
e_col = Calc_e(pars,c_col);

Ne      = numel(mu_e);
omg_col = zeros(ny,size(STATE.p,2),Np);
for ip = 1:Np
    omg_col(:,:,ip) = reshape(PhaseG(pars{ip},c{ip}),ny,[]);
    for ie = 1:Ne
        omg_col(:,:,ip) = omg_col(:,:,ip) - e_col{ip}{ie} .* mu_e{ie};
    end
end

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

% ---------- Diagnostics ----------
total_iter              = newton_iter + proj_iter;
DIAG.max_stationarity   = max(stat_node);
DIAG.stat_node          = reshape(stat_node,ny,[]);
DIAG.iter_max           = max(total_iter);          % backwards-compatible total work
DIAG.iter_node          = reshape(total_iter,ny,[]);
DIAG.max_newton_iter    = max(newton_iter);
DIAG.max_project_iter   = max(proj_iter);
DIAG.newton_iter        = reshape(newton_iter,ny,[]);
DIAG.project_iter       = reshape(proj_iter,ny,[]);
DIAG.active             = active;
DIAG.p_le               = p_le;
DIAG.spinodal_node      = reshape(bad_node,ny,[]);
DIAG.spinodal_fraction  = mean(bad_node);
DIAG.bad_phase_weight   = reshape(bad_weight,ny,[]);
DIAG.project_stalled    = reshape(proj_stalled,ny,[]);
DIAG.newton_stalled     = reshape(newton_stalled,ny,[]);
DIAG.preserve_E_mode    = preserve_E;
DIAG.compute_output_chi  = opt.compute_chi;
DIAG.seed_spinodal_fraction = mean(seed_bad);

end


% ========================================================================
% Local solver
% ========================================================================

function [c,chi,D] = Solve_Local_FromMu( ...
    pars,p,c,mu_e,eta,E_in,bad_seed,bad_weight_seed,opt)
%SOLVE_LOCAL_FROMMU Mode-aware local correction.
%
% Regular GP closure:
%   Safe nodes converge with exact fixed-mu Newton.
%   Unsafe tails receive a small bounded mixed projection.
%
% Conserved-E spinodal closure:
%   Unsafe nodes are not asked to satisfy fixed-mu stationarity.
%   Their composition is corrected only toward E_in, normally in one sweep.

Ne = numel(mu_e);
N  = numel(mu_e{1});

D.stat            = zeros(1,N);
D.newton_iter     = zeros(1,N);
D.proj_iter       = zeros(1,N);
D.bad             = bad_seed;
D.bad_weight      = bad_weight_seed;
D.proj_stalled    = false(1,N);
D.newton_stalled  = false(1,N);

if ~HasCompositionDOF(pars,c)
    if opt.compute_chi
        chi = EyeOverEta(eta,Ne,N);
    else
        chi = [];
    end
    return
end

bad_locked = bad_seed;

% ---------- A. Conserved-E mode: composition follows E everywhere ----------
% Once the whole-step closure is conserved-E, a nodewise mixture of
% fixed-mu Newton and E-projection is inconsistent.  Use one inexpensive
% E-only correction for all nodes and return.
if opt.preserve_E
    live = true(1,N);
    for it = 1:opt.conserved_iter
        merit_old = E_Closure_Merit(pars,p,c,mu_e,eta,E_in);
        live = live & merit_old > opt.E_tol;
        if ~any(live)
            break
        end

        dc = Block_EClosure_Step(pars,p,c,mu_e,eta,E_in,live,opt);
        dc = CapStepByNode(dc,live,opt.dc_cap_conserved);

        trial_opt = opt;
        trial_opt.max_proj_ls = opt.conserved_ls;
        [c,accepted,merit_new] = Accept_Projected_Step( ...
            pars,p,c,mu_e,eta,E_in,dc,live,merit_old,trial_opt);

        D.proj_iter(live) = it;
        failed = live & ~accepted;
        D.proj_stalled(failed) = true;
        live(failed) = false;

        if any(accepted)
            rel_drop = (merit_old-merit_new) ./ max(merit_old,eps);
            done = live & accepted & ...
                (merit_new <= opt.E_tol | rel_drop <= opt.proj_min_rel_drop);
            live(done) = false;
        end
    end
    D.proj_stalled(live) = true;

    if opt.compute_chi
        chi = Calc_Chi_Hybrid(pars,p,c,eta,Ne,opt);
    else
        chi = [];
    end
    D.stat = Mu_Stationarity_Residual(pars,c,mu_e);
    return
end

% ---------- B. Regular GP mode: exact fixed-mu Newton for safe nodes ----------
live = ~bad_locked;

for it = 1:opt.max_newton_iter
    if ~any(live)
        break
    end

    dc = Exact_Newton_Step(pars,c,mu_e,live);
    dc = CapStepByNode(dc,live,opt.dc_cap_good);

    [c,accepted,step_taken] = Accept_Newton_Step( ...
        pars,p,c,mu_e,dc,live,opt);

    D.newton_iter(live) = it;

    failed = live & ~accepted;
    D.newton_stalled(failed) = true;
    live(failed) = false;

    if any(accepted)
        stat = Mu_Stationarity_Residual(pars,c,mu_e);
        done = live & accepted & ...
            (step_taken < opt.c_tol | stat < opt.stat_tol);
        live(done) = false;
    end

    % Reclassify only the nodes that are still undergoing Newton.
    % A crossed node is locked into projection for the remainder of this call.
    if any(live)
        [bad_now,bad_w_now] = Classify_Bad_Nodes(pars,p,c,opt);
        crossed = live & bad_now;
        bad_locked(crossed) = true;
        live(crossed) = false;
        D.bad_weight = max(D.bad_weight,bad_w_now);
    end
end

D.newton_stalled(live) = true;
D.bad = bad_locked;

% ---------- C. Regular-GP fallback for small unsafe tails ----------
live = bad_locked;

for it = 1:opt.max_proj_iter
    if ~any(live)
        break
    end

    merit_old = E_Closure_Merit(pars,p,c,mu_e,eta,E_in);
    already_done = live & merit_old <= opt.E_tol;
    live(already_done) = false;
    if ~any(live)
        break
    end

    % Regular-GP regime: only small unsafe tails reach this branch.
    dc = Block_Projected_Step(pars,p,c,mu_e,eta,E_in,live,opt);
    trial_opt = opt;

    dc = CapStepByNode(dc,live,opt.dc_cap_bad);

    [c,accepted,merit_new] = Accept_Projected_Step( ...
        pars,p,c,mu_e,eta,E_in,dc,live,merit_old,trial_opt);

    D.proj_iter(live) = it;

    failed = live & ~accepted;
    D.proj_stalled(failed) = true;
    live(failed) = false;

    if any(accepted)
        rel_drop = (merit_old-merit_new) ./ max(merit_old,eps);
        done = live & accepted & ...
            (merit_new <= opt.E_tol | rel_drop <= opt.proj_min_rel_drop);
        live(done) = false;
    end
end

D.proj_stalled(live) = true;

if opt.compute_chi
    chi = Calc_Chi_Hybrid(pars,p,c,eta,Ne,opt);
else
    chi = [];
end
D.stat = Mu_Stationarity_Residual(pars,c,mu_e);

end


function [bad,bad_weight] = Classify_Bad_Nodes(pars,p,c,opt)
Np = numel(c);
N  = numel(c{1}{1});

bad        = false(1,N);
bad_weight = zeros(1,N);

for ip = 1:Np
    R = PhaseThermo(pars{ip},c{ip});
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end
    unsafe = ~HessianSafeMask(SymmetricPages(R.H_c),opt.h_rel_tol,opt.cond_max);
    p_ip   = reshape(p(:,:,ip),1,N);
    bad        = bad | unsafe;
    bad_weight = max(bad_weight,p_ip .* unsafe);
end
end


function dc = Exact_Newton_Step(pars,c,mu_e,mask)
Np     = numel(c);
Ne     = numel(mu_e);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);

dc = ZeroStepLike(c);

for ip = 1:Np
    R = PhaseThermo(pars{ip},c{ip});
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    H    = SymmetricPages(R.H_c);
    JT   = permute(R.Jac,[2 1 3]);
    Nc   = size(H,1);
    Rmu  = reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),Nc,N) ...
           - StackFields(R.mu_c,N);
    d    = reshape(pagemldivide(H(:,:,mask), ...
           reshape(Rmu(:,mask),Nc,1,[])),Nc,[]);

    for ic = 1:Nc
        dc{ip}{ic}(mask) = d(ic,:);
    end
end
end


function dc = Block_EClosure_Step(pars,p,c,mu_e,eta,E_in,mask,opt)
%BLOCK_ECLOSURE_STEP Composition correction in conserved-E spinodal mode.
%
% Solve only
%       min ||G*dc - (E_in - E_mu)||^2 + lam_E ||dc||^2
% because E is the conserved variable in this mode.  The fixed-mu
% stationarity residual is intentionally excluded.

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

Rall = cell(1,Np);
idx  = cell(1,Np);
Nc_total = 0;
Emix = zeros(Ne,N);

for ip = 1:Np
    Rall{ip} = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + bsxfun(@times,StackFields(Rall{ip}.e,N),p_ip);

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

R_E = E_mat - Emix - bsxfun(@rdivide,mu_mat,eta_v);
G   = zeros(Ne,Nc_total,Nb);

for ip = 1:Np
    id = idx{ip};
    if isempty(id)
        continue
    end
    J     = Rall{ip}.Jac(:,:,mask);
    p_ip  = reshape(p(:,:,ip),1,N);
    p_sel = reshape(p_ip(mask),1,1,Nb);
    G(:,id,:) = bsxfun(@times,J,p_sel);
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


function dc = Block_Projected_Step(pars,p,c,mu_e,eta,E_in,mask,opt)
% Solve the E-guided least-squares step only on selected non-convex nodes:
%
% min ||G dc - R_E||^2 + w_mu ||B dc - R_mu||^2 + lam ||dc||^2
%
% G = [p_1 J_1, p_2 J_2, ...] is the derivative of mixture E.

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

Rall = cell(1,Np);
idx  = cell(1,Np);
Rmu  = cell(1,Np);
Nc_total = 0;
Emix = zeros(Ne,N);

for ip = 1:Np
    Rall{ip} = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + bsxfun(@times,StackFields(Rall{ip}.e,N),p_ip);

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

R_E = E_mat - Emix - bsxfun(@rdivide,mu_mat,eta_v);
G   = zeros(Ne,Nc_total,Nb);
T   = zeros(Nc_total,Nc_total,Nb);
r   = zeros(Nc_total,1,Nb);

for ip = 1:Np
    id = idx{ip};
    if isempty(id)
        continue
    end

    J     = Rall{ip}.Jac(:,:,mask);
    H     = SymmetricPages(Rall{ip}.H_c(:,:,mask));
    HT    = permute(H,[2 1 3]);
    p_ip  = reshape(p(:,:,ip),1,N);
    p_sel = reshape(p_ip(mask),1,1,Nb);

    G(:,id,:)  = bsxfun(@times,J,p_sel);
    T(id,id,:) = pagemtimes(HT,H);
    r(id,1,:)  = pagemtimes(HT,reshape(Rmu{ip}(:,mask),numel(id),1,Nb));
end

GT = permute(G,[2 1 3]);
A  = pagemtimes(GT,G) + opt.wmu_proj*T + ...
     repmat(eye(Nc_total),1,1,Nb).*opt.lam_proj;
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


function [c,accepted,step_taken] = Accept_Newton_Step(pars,p,c,mu_e,dc,mask,opt)
N          = numel(mu_e{1});
F0         = LE_Mu_Objective(pars,p,c,mu_e);
accepted   = false(1,N);
alpha_acc  = zeros(1,N);
alpha_try  = opt.alpha*ones(1,N);
alpha_try(~mask) = 0;

for ils = 1:opt.max_ls
    test = mask & ~accepted & alpha_try >= opt.alpha_min;
    if ~any(test)
        break
    end

    c_try = AddStep(c,dc,alpha_try);
    F_try = LE_Mu_Objective(pars,p,c_try,mu_e);

    take = test & isfinite(F_try) & ...
        (F_try <= F0 + opt.energy_tol.*max(1,abs(F0)));

    accepted(take)  = true;
    alpha_acc(take) = alpha_try(take);
    alpha_try(test & ~take) = 0.2*alpha_try(test & ~take);
end

c = AddStep(c,dc,alpha_acc);
step_taken = NodeMaxStep(dc).*alpha_acc;
end


function [c,accepted,merit_new] = Accept_Projected_Step( ...
    pars,p,c,mu_e,eta,E_in,dc,mask,merit_old,opt)

N          = numel(mu_e{1});
accepted   = false(1,N);
alpha_acc  = zeros(1,N);
alpha_try  = opt.alpha*ones(1,N);
alpha_try(~mask) = 0;
merit_new  = merit_old;

for ils = 1:opt.max_proj_ls
    test = mask & ~accepted & alpha_try >= opt.alpha_min;
    if ~any(test)
        break
    end

    c_try = AddStep(c,dc,alpha_try);
    m_try = E_Closure_Merit(pars,p,c_try,mu_e,eta,E_in);

    decrease = m_try < merit_old - opt.proj_accept_abs.*max(1,merit_old);
    take = test & isfinite(m_try) & (decrease | m_try <= opt.E_tol);

    accepted(take)  = true;
    alpha_acc(take) = alpha_try(take);
    merit_new(take) = m_try(take);
    alpha_try(test & ~take) = 0.2*alpha_try(test & ~take);
end

c = AddStep(c,dc,alpha_acc);
end


function merit = E_Closure_Merit(pars,p,c,mu_e,eta,E_in)
E_mu = Calc_E_FromMu(pars,c,p,mu_e,eta);
N    = numel(mu_e{1});
R    = StackFields(E_in,N) - StackFields(E_mu,N);
merit = sqrt(sum(R.^2,1));
end


% ========================================================================
% Constitutive tangent and E closure
% ========================================================================

function chi = Calc_Chi_Hybrid(pars,p,c,eta,Ne,opt)
% Exact thermodynamic tangent at convex nodes; bounded projection tangent at
% non-convex nodes. When conserved-E closure is selected the next main-loop
% LE_Run normally overwrites this chi before it enters the PDE solve.

Np    = numel(c);
N     = numel(c{1}{1});
eta_v = EtaVector(eta,N);

Rall = cell(1,Np);
idx  = cell(1,Np);
Nc_total = 0;
bad = false(1,N);

for ip = 1:Np
    Rall{ip} = PhaseThermo(pars{ip},c{ip});
    if isempty(Rall{ip}.mu_c) || isempty(Rall{ip}.H_c) || isempty(Rall{ip}.Jac)
        idx{ip} = [];
        continue
    end
    Nc = size(Rall{ip}.H_c,1);
    idx{ip} = Nc_total + (1:Nc);
    Nc_total = Nc_total + Nc;
    bad = bad | ~HessianSafeMask(SymmetricPages(Rall{ip}.H_c), ...
                                  opt.h_rel_tol,opt.cond_max);
end

chi_page = repmat(eye(Ne),1,1,N).*reshape(1./eta_v,1,1,N);
if Nc_total == 0
    chi = PagesToCell(chi_page,Ne);
    return
end

good = ~bad;

if any(good)
    C = zeros(Ne,Ne,nnz(good));
    for ip = 1:Np
        if isempty(idx{ip}), continue; end
        J      = Rall{ip}.Jac(:,:,good);
        H      = SymmetricPages(Rall{ip}.H_c(:,:,good));
        p_ip   = reshape(p(:,:,ip),1,N);
        p_good = reshape(p_ip(good),1,1,[]);
        C = C + bsxfun(@times, ...
            pagemtimes(J,pagemldivide(H,permute(J,[2 1 3]))),p_good);
    end
    chi_page(:,:,good) = chi_page(:,:,good) + C;
end

if any(bad)
    Nb = nnz(bad);
    G  = zeros(Ne,Nc_total,Nb);
    T  = zeros(Nc_total,Nc_total,Nb);
    Q  = zeros(Nc_total,Ne,Nb);

    for ip = 1:Np
        id = idx{ip};
        if isempty(id), continue; end
        J     = Rall{ip}.Jac(:,:,bad);
        H     = SymmetricPages(Rall{ip}.H_c(:,:,bad));
        p_ip  = reshape(p(:,:,ip),1,N);
        p_bad = reshape(p_ip(bad),1,1,Nb);
        G(:,id,:)  = bsxfun(@times,J,p_bad);
        T(id,id,:) = pagemtimes(permute(H,[2 1 3]),H);
        Q(id,:,:)  = pagemtimes(permute(H,[2 1 3]),permute(J,[2 1 3]));
    end

    GT = permute(G,[2 1 3]);
    A  = pagemtimes(GT,G) + opt.wmu_proj*T + ...
         repmat(eye(Nc_total),1,1,Nb).*opt.lam_proj;
    dB = -bsxfun(@rdivide,GT,reshape(eta_v(bad),1,1,Nb)) + opt.wmu_proj*Q;

    chi_page(:,:,bad) = chi_page(:,:,bad) + ...
        pagemtimes(G,pagemldivide(A,dB));
end

chi = PagesToCell(chi_page,Ne);
end


function preserve_E = Select_E_Closure(old_mode,bad_weight,opt)
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


function E = Calc_E_FromMu(pars,c,p,mu_e,eta)
Np    = numel(c);
Ne    = numel(mu_e);
N     = numel(mu_e{1});
eta_v = EtaVector(eta,N);
Emix  = zeros(Ne,N);

for ip = 1:Np
    R    = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + bsxfun(@times,StackFields(R.e,N),p_ip);
end

E = cell(1,Ne);
for ie = 1:Ne
    E{ie} = Emix(ie,:) + reshape(mu_e{ie},1,N)./eta_v;
end
end


% ========================================================================
% Active phases and controls
% ========================================================================

function opt = Read_FromMu_Options(PARAM)
opt.Pmax             = GetField(PARAM,'frommu_Pmax',4);
opt.p_tail           = GetField(PARAM,'frommu_p_tail',1e-3);
opt.p_full           = GetField(PARAM,'frommu_p_full',5e-2);
opt.p_on             = GetField(PARAM,'frommu_p_on',5e-3);
opt.p_off            = GetField(PARAM,'frommu_p_off',2e-3);
opt.alpha_LE         = GetField(PARAM,'frommu_alpha_LE',[0.8 0.5 0.4 0.3]);
opt.max_newton_iter  = GetField(PARAM,'frommu_max_newton_iter',[100 100 100 100]);

opt.c_tol             = GetField(PARAM,'frommu_c_tol',1e-7);
opt.stat_tol          = GetField(PARAM,'frommu_stat_tol',1e-7);
opt.max_ls            = GetField(PARAM,'frommu_max_ls',8);
opt.energy_tol        = GetField(PARAM,'frommu_energy_tol',1e-10);
opt.alpha_min         = GetField(PARAM,'frommu_alpha_min',1e-8);
opt.dc_cap_good       = GetField(PARAM,'frommu_dc_cap_good',5e-3);

% Key speed controls for nonconvex nodes:
opt.max_proj_iter     = GetField(PARAM,'frommu_max_proj_iter',2);
opt.max_proj_ls       = GetField(PARAM,'frommu_max_proj_ls',3);
opt.dc_cap_bad        = GetField(PARAM,'frommu_dc_cap_bad',5e-4);
opt.E_tol             = GetField(PARAM,'frommu_E_tol',1e-8);
opt.proj_min_rel_drop = GetField(PARAM,'frommu_proj_min_rel_drop',1e-3);
opt.proj_accept_abs   = GetField(PARAM,'frommu_proj_accept_abs',1e-12);

% In conserved-E spinodal mode, one E-only correction is normally enough
% for a single solution phase because e(c) is locally affine.
opt.conserved_iter    = GetField(PARAM,'frommu_conserved_iter',1);
opt.conserved_ls      = GetField(PARAM,'frommu_conserved_ls',2);
opt.dc_cap_conserved  = GetField(PARAM,'frommu_dc_cap_conserved',5e-4);
opt.lam_E             = GetField(PARAM,'frommu_lam_E',1e-12);

% Switch away from exact H-inversion before a positive curvature becomes
opt.h_rel_tol         = GetField(PARAM,'frommu_h_rel_tol',1e-6);
opt.cond_max          = GetField(PARAM,'frommu_cond_max',1e8);
opt.lam_proj          = GetField(PARAM,'frommu_lam_proj',1e-8);
opt.wmu_proj          = GetField(PARAM,'frommu_wmu_proj',1);

% With the stated loop, output chi is overwritten by LE_Run at the next
% step before PF_Coupled_ACCH_LETangent uses it. Skip its expensive rebuild.
opt.compute_chi       = logical(GetField(PARAM,'frommu_compute_chi',false));

opt.p_enter           = GetField(PARAM,'frommu_p_enter',0.80);
opt.p_exit            = GetField(PARAM,'frommu_p_exit',0.60);
opt.E_mode            = lower(GetField(PARAM,'frommu_E_mode','auto'));
opt.use_WScale        = isfield(PARAM,'use_WScale') && PARAM.use_WScale;

assert(opt.p_off < opt.p_on, 'frommu_p_off must be smaller than frommu_p_on.');
assert(opt.p_tail < opt.p_full, 'frommu_p_tail must be smaller than frommu_p_full.');
assert(numel(opt.alpha_LE) >= opt.Pmax && numel(opt.max_newton_iter) >= opt.Pmax, ...
    'frommu alpha/iteration controls do not cover Pmax.');
end


function [p_le,active,LE_state] = Build_Active_Phases(p,LE_state,opt)
Np = size(p,3);
N  = size(p,2);

p_th   = Calc_Thermo_p(p,opt.p_tail,opt.p_full);
pThMat = reshape(p_th,N,Np);

if ~isfield(LE_state,'active') || ~isequal(size(LE_state.active),[N,Np])
    active_old = false(N,Np);
else
    active_old = LE_state.active;
end

active = (active_old & pThMat > opt.p_off) | ...
         (~active_old & pThMat > opt.p_on);

repair = find(~any(active,2) | sum(active,2) > opt.Pmax).';
for i = repair
    if ~any(active(i,:))
        [~,imax] = max(pThMat(i,:));
        active(i,imax) = true;
    end
    if sum(active(i,:)) > opt.Pmax
        score = pThMat(i,:) + 0.5*opt.p_on*active_old(i,:);
        [~,ord] = sort(score,'descend');
        active(i,:) = false;
        active(i,ord(1:opt.Pmax)) = true;
    end
end

p_le = p_th .* reshape(active,1,N,Np);
p_le = p_le ./ max(sum(p_le,3),eps);

LE_state.active = active;
LE_state.p_th   = p_th;
LE_state.p_le   = p_le;
end


function value = GetField(S,name,default)
if isfield(S,name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default;
end
end


% ========================================================================
% Generic local algebra and utility routines
% ========================================================================

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


function H = SymmetricPages(H)
H = 0.5*(H + permute(H,[2 1 3]));
end


function safe = HessianSafeMask(H,h_rel_tol,cond_max)
N = size(H,3);
safe = false(1,N);
for i = 1:N
    ev = eig(H(:,:,i));
    scale = max(1,max(abs(ev)));
    safe(i) = min(ev) > h_rel_tol*scale && max(abs(ev))/max(min(abs(ev)),eps*scale) < cond_max;
end
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


function out = Slice_Chi_Local(chi,mask)
% Return the incoming chi subset without recomputing a post-projection tangent.
Ne = size(chi,1);
out = cell(Ne,Ne);
for i = 1:Ne
    for j = 1:Ne
        x = reshape(chi{i,j},1,[]);
        out{i,j} = x(mask);
    end
end
end


function chi = EyeOverEta(eta,Ne,N)
eta_v = EtaVector(eta,N);
chi = PagesToCell(repmat(eye(Ne),1,1,N).*reshape(1./eta_v,1,1,N),Ne);
end


function chi = PagesToCell(A,Ne)
chi = cell(Ne,Ne);
for i = 1:Ne
    for j = 1:Ne
        chi{i,j} = reshape(A(i,j,:),1,[]);
    end
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


function dc = CapStepByNode(dc,mask,cap)
if ~any(mask), return; end
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
N = numel(dc{1}{1});
amp = zeros(1,N);
for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        amp = max(amp,abs(reshape(dc{ip}{ic},1,N)));
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


function stat = Mu_Stationarity_Residual(pars,c,mu_e)
Np     = numel(c);
Ne     = numel(mu_e);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);
stat   = zeros(1,N);

for ip = 1:Np
    R = PhaseThermo(pars{ip},c{ip});
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end
    mu_c = StackFields(R.mu_c,N);
    JTmu = reshape(pagemtimes(permute(R.Jac,[2 1 3]), ...
        reshape(mu_mat,Ne,1,N)),size(mu_c,1),N);
    stat = max(stat,max(abs(mu_c-JTmu),[],1));
end
end


function F = LE_Mu_Objective(pars,p,c,mu_e)
Np     = numel(c);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);
F      = zeros(1,N);

try
    for ip = 1:Np
        R    = PhaseThermo(pars{ip},c{ip});
        p_ip = reshape(p(:,:,ip),1,N);
        omg  = R.g(:).' - sum(StackFields(R.e,N).*mu_mat,1);
        F    = F + p_ip.*omg;
    end
catch
    F = inf(1,N);
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
x    = min(max((p-p_tail)./(p_full-p_tail),0),1);
w    = p .* x.^2 .* (3-2*x);
wsum = sum(w,3);
p_th = w ./ max(wsum,eps);

bad = wsum <= eps;
if any(bad(:))
    [~,imax] = max(p,[],3);
    for ip = 1:size(p,3)
        tmp = p_th(:,:,ip);
        tmp(bad & imax == ip) = 1;
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
p_c    = zeros(1,N,Nphase);
c_c    = cell(1,Nphase);

for iph = 1:Nphase
    grains      = find(grain_to_phase == iph);
    ig0         = grains(1);
    pars_c{iph} = pars{ig0};

    for ig = grains
        p_c(:,:,iph) = p_c(:,:,iph) + p(:,:,ig);
    end

    Nc = numel(c{ig0});
    c_c{iph} = cell(1,Nc);
    den = reshape(sum(p(:,:,grains),3),size(c{ig0}{1}));
    valid = den > eps;

    for ic = 1:Nc
        num = zeros(size(c{ig0}{ic}));
        for ig = grains
            num = num + reshape(p(:,:,ig),size(num)).*c{ig}{ic};
        end
        tmp = c{ig0}{ic};
        tmp(valid) = num(valid)./den(valid);
        c_c{iph}{ic} = tmp;
    end
end
end


function d = Max_Cell_Diff_Local(A,B)
d = 0;
for i = 1:numel(A)
    d = max(d,max(abs(A{i}(:)-B{i}(:))));
end
end
