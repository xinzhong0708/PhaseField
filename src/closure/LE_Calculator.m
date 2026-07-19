function [c,mu_e,chi,DIAG] = LE_Calculator(pars,p,c,E,eta,level)
%LE_CALCULATOR Local equilibrium with branch-tracking damped update.
%
% This version is stabilized for spinodal/interface cases by using
% residual decrease rather than energy decrease in the line search.
%
% Fast v3 change:
%   The residual line search uses masked local blocks and reuses the
%   already-computed current correction as R_old. This keeps the v2
%   numerics but avoids repeated full-domain PhaseThermo calls.
%
% Main stabilization ideas:
%   1) The line search accepts a trial only when the LE stationarity
%      residual does not increase. It does not force energy decrease.
%   2) The Newton/Picard step is capped component-wise to avoid branch jumps.
%   3) Optional E-continuation can move from the current mixture E to the
%      requested E in several substeps.
%   4) Failed stages preserve the input composition instead of writing a
%      partially converged branch back to STATE.
%
% Optional level entries:
%   level = [alpha, Miter]
%   level = [alpha, Miter, c_tol, MaxLS, lam_c]
%   level = [alpha, Miter, c_tol, MaxLS, lam_c, c_bound_margin]
%   level = [alpha, Miter, c_tol, MaxLS, lam_c, c_bound_margin, dc_cap, E_substep, fast_accept_dc]

%Prepare
c_init       =  c;
Np           =  length(c);
Ne           =  length(E);
N            =  numel(E{1});

if nargin < 6 || isempty(level)
    level = [0.6,100];
end

alpha        =  level(1);
Miter        =  level(2);
c_tol        =  1e-6;

%Line-search & damping parameters
MaxLS        =  10;
amin         =  1e-7;
res_tol      =  1e-10;
lam_c        =  1e-7;

%Composition bounds for LE trial and accepted states.
%Trials that hit/exceed these artificial bounds are rejected in the
%line search instead of being clipped to +/-1 and accepted.
c_min        = -1.0;
c_max        =  1.0;
c_bound_margin = 1e-12;

%Branch tracking parameters
dc_cap       = 1e-2;
E_substep    = 1;
fast_accept_dc = NaN;

DIAG                 = struct();
DIAG.converged       = true;
DIAG.failed          = false;
DIAG.has_dof         = true;
DIAG.n_nodes         = N;
DIAG.alpha_initial   = alpha;
DIAG.alpha_stage     = NaN;
DIAG.n_iter          = 0;
DIAG.max_dc          = 0;
DIAG.max_cchg        = 0;
DIAG.max_res         = NaN;
DIAG.inexact         = false;
DIAG.bound_hits      = 0;
DIAG.c_tol           = c_tol;
DIAG.lam_c           = lam_c;
DIAG.dc_cap          = dc_cap;
DIAG.E_substep       = E_substep;
DIAG.fast_accept_dc  = fast_accept_dc;
DIAG.message         = '';

if numel(level) >= 3 && ~isempty(level(3))
    c_tol = level(3);
end
if numel(level) >= 4 && ~isempty(level(4))
    MaxLS = level(4);
end
if numel(level) >= 5 && ~isempty(level(5))
    lam_c = level(5);
end
if numel(level) >= 6 && ~isempty(level(6))
    c_bound_margin = level(6);
end
if numel(level) >= 7 && ~isempty(level(7))
    dc_cap = level(7);
end
if numel(level) >= 8 && ~isempty(level(8))
    E_substep = max(1,round(level(8)));
end
if numel(level) >= 9 && ~isempty(level(9))
    fast_accept_dc = level(9);
end
if ~isfinite(fast_accept_dc)
    fast_accept_dc = 0.25*dc_cap;
end

DIAG.c_tol     = c_tol;
DIAG.lam_c     = lam_c;
DIAG.dc_cap    = dc_cap;
DIAG.E_substep = E_substep;
DIAG.fast_accept_dc = fast_accept_dc;

%Make sure the input state is already bounded
c_init = BoundC_Local(c_init,c_min,c_max);
c      = c_init;

%Check whether any phase has internal composition degrees of freedom
has_dof = false(1,Np);
for ip = 1:Np
    has_dof(ip) = HasDofFromPars(pars{ip});
end

%If all phases are pure/no-DOF phases, no c-iteration is needed
if ~any(has_dof)
    [~,mu_mat,chi_page] = LE_Quadratic_Step_Fast(pars,p,c,E,eta,0);
    mu_e = CellFromMat_Local(mu_mat,Ne);
    chi  = CellFromPage_Local(chi_page,Ne);
    DIAG.has_dof = false;
    DIAG.alpha_stage = 0;
    return
end

% -------------------------------------------------------------------------
% E-continuation start.
% Since LE_Calculator receives only the new target E, the safest available
% start is the element mixture represented by c_init. If E_substep = 1, this
% is not used and the code behaves as direct fixed-E LE.
% -------------------------------------------------------------------------
if E_substep > 1
    E_start = LE_MixtureElement_Local(pars,p,c_init,E);
else
    E_start = E;
end

% -------------------------------------------------------------------------
% Non-recursive damping schedule.
% If the original alpha fails, retry with smaller alpha, but never call
% LE_Calculator recursively. If all stages fail, preserve c_init.
% -------------------------------------------------------------------------
alpha_stage0 = [alpha, min(alpha,0.3), min(alpha,0.1), min(alpha,0.03)];
alpha_stage  = [];

for ia = 1:numel(alpha_stage0)
    if alpha_stage0(ia) > 0 && ~any(abs(alpha_stage-alpha_stage0(ia)) < eps)
        alpha_stage(end+1) = alpha_stage0(ia); %#ok<AGROW>
    end
end

converged_all = false;
last_alpha    = NaN;
last_iter     = 0;
last_dcmax    = 0;
last_cchg     = 0;
last_res      = NaN;
last_bound_hits = 0;
last_inexact = false;

for istage = 1:numel(alpha_stage)

    alpha_now = alpha_stage(istage);
    last_alpha = alpha_now;

    if istage == numel(alpha_stage)
        Miter_now = max(Miter,300);
    else
        Miter_now = Miter;
    end

    % Always restart from the input c for a smaller-alpha retry.
    % This avoids carrying a bad partially converged branch.
    c = c_init;

    failed_stage = false;

    for isub = 1:E_substep

        theta = isub / E_substep;
        E_now = BlendE_Local(E_start,E,theta);

        [c,converged_sub,ST] = LE_Solve_Target(pars,p,c,E_now,eta,alpha_now, ...
            Miter_now,c_tol,MaxLS,amin,res_tol,lam_c,c_min,c_max, ...
            c_bound_margin,dc_cap,fast_accept_dc,has_dof);

        last_iter       = ST.n_iter;
        last_dcmax      = ST.max_dc;
        last_cchg       = ST.max_cchg;
        last_res        = ST.max_res;
        last_bound_hits = last_bound_hits + ST.bound_hits;
        last_inexact    = last_inexact | ST.inexact;

        if ~converged_sub
            failed_stage = true;
            break
        end
    end

    if ~failed_stage
        converged_all = true;
        break
    end

    if failed_stage && istage < numel(alpha_stage)
        fprintf('LE_Calculator: retry with smaller alpha = %.3e\n',alpha_stage(istage+1));
    end
end

%If all stages fail, preserve previous input composition.
%This prevents infinite recursion and avoids accepting a bad partial branch.
if ~converged_all
    disp('LE_Calculator: not converged after branch-tracking retries, preserving input c and returning DIAG.failed = true.')
    c = c_init;
end

DIAG.converged   = converged_all;
DIAG.failed      = ~converged_all;
DIAG.alpha_stage = last_alpha;
DIAG.n_iter      = last_iter;
DIAG.max_dc      = last_dcmax;
DIAG.max_cchg    = last_cchg;
DIAG.max_res     = last_res;
DIAG.bound_hits  = last_bound_hits;
DIAG.inexact     = last_inexact;
DIAG.c_tol       = c_tol;
DIAG.lam_c       = lam_c;
DIAG.dc_cap      = dc_cap;
DIAG.E_substep   = E_substep;
DIAG.fast_accept_dc = fast_accept_dc;

if ~converged_all
    if last_bound_hits > 0
        DIAG.message = sprintf('not converged after branch-tracking retries; rejected %d bound-hit trial nodes and preserved input c',last_bound_hits);
    else
        DIAG.message = 'not converged after branch-tracking retries, preserving input c';
    end
elseif last_inexact
    DIAG.message = 'bounded inexact branch-tracking LE update accepted';
elseif last_bound_hits > 0
    DIAG.message = sprintf('converged after rejecting %d bound-hit trial nodes',last_bound_hits);
end

%Safety bound before final thermodynamic response
c = BoundC_Local(c,c_min,c_max);

%Recalculate final mu_e and chi at accepted c using the same regularized
%Hessian family used by the LE iteration. Returning raw-Hessian chi here can
%hand AC-CH/CH-LE a different local tangent than the one LE converged with.
try
    [~,mu_mat,chi_page] = LE_Quadratic_Step_Fast(pars,p,c,E,eta,lam_c);
catch ME
    error('LE_Calculator: regularized final mu_e/chi response failed: %s',ME.message)
end

%If final response produced non-finite values, fail fast instead of
%returning a misleading local tangent.
if any(~isfinite(mu_mat(:))) || any(~isfinite(chi_page(:)))
    error('LE_Calculator: non-finite final mu_e/chi response.')
end

%Diffusion potential and susceptibility
mu_e = CellFromMat_Local(mu_mat,Ne);
chi  = CellFromPage_Local(chi_page,Ne);

end

function [c,converged,ST] = LE_Solve_Target(pars,p,c,E,eta,alpha_now, ...
    Miter,c_tol,MaxLS,amin,res_tol,lam_c,c_min,c_max,c_bound_margin,dc_cap,fast_accept_dc,has_dof)

N = numel(E{1});

ST = struct();
ST.n_iter = 0;
ST.max_dc = 0;
ST.max_cchg = 0;
ST.max_res = NaN;
ST.inexact = false;
ST.bound_hits = 0;

converged = false;
failed    = false;

for it = 1:Miter

    ST.n_iter = it;

    %Save old c
    c_old = c;

    %Quadratic branch-tracking step.
    [dc_all,~,~] = LE_Quadratic_Step_Fast(pars,p,c_old,E,eta,lam_c);

    dcmax_raw = MaxAbsStep(dc_all);
    ST.max_dc = dcmax_raw;

    if dcmax_raw < c_tol
        converged = true;
        break
    end

    %Residual at the current state. This is the globalization target.
    %It is the next regularized LE correction size. The current correction
    %has already been computed above, so do not call PhaseThermo again here.
    R_old = MaxAbsStepNode(dc_all,N);
    R_old(~isfinite(R_old)) = inf;
    tmp_res = R_old(isfinite(R_old));
    if isempty(tmp_res)
        ST.max_res = inf;
    else
        ST.max_res = max(tmp_res);
    end

    %Trust-region cap on the nonlinear update.
    dc_all = CapStep_Local(dc_all,dc_cap);
    dcnode = MaxAbsStepNode(dc_all,N);

    %Backtracking line search, independently for each grid point
    good_node = false(1,N);
    alpha_try = alpha_now*ones(1,N);
    alpha_acc = zeros(1,N);

    for ils = 1:MaxLS

        bad = ~good_node & alpha_try >= amin;

        if ~any(bad)
            break
        end

        %Add trial for all nodes.
        %Do not clip an out-of-bound trial to +/-1 and accept it.
        c_try_raw = AddStep(c_old,dc_all,alpha_try);
        bound_bad = BoundHitNode_Local(pars,c_try_raw,c_min,c_max,c_bound_margin);
        ST.bound_hits = ST.bound_hits + nnz(bad & bound_bad);

        %Small bounded updates are accepted directly. These nodes are
        %not where branch jumps happen, and skipping their residual check
        %avoids many PhaseThermo calls.
        step_node = alpha_try .* dcnode;
        fast_ok   = bad & ~bound_bad & fast_accept_dc > 0 & step_node <= fast_accept_dc;

        if any(fast_ok)
            alpha_acc(fast_ok) = alpha_try(fast_ok);
            good_node(fast_ok) = true;
        end

        %Bound only for safe residual evaluation. Bound-hit nodes are
        %assigned residual = inf below and therefore cannot be accepted.
        c_try = BoundC_Local(c_try_raw,c_min,c_max);

        R_try_bad = inf(1,nnz(bad));
        eval_mask = bad & ~bound_bad & ~fast_ok;

        if any(eval_mask)
            R_eval = LE_ModelResidual_MaskFast(pars,p,c_try,E,eta,lam_c,eval_mask);
            bad_ids_all = find(bad);
            pos_eval = find(eval_mask(bad_ids_all));
            R_try_bad(pos_eval) = R_eval;
        end

        R_old_bad = R_old(bad);
        fast_bad = fast_ok(bad);
        good_bad  = fast_bad | (isfinite(R_try_bad) & ...
                    (R_try_bad <= R_old_bad.*(1 + 1e-4) + res_tol.*max(1,abs(R_old_bad))));

        bad_ids  = find(bad);
        good_ids = bad_ids(good_bad);
        fail_ids = bad_ids(~good_bad);

        alpha_acc(good_ids) = alpha_try(good_ids);
        good_node(good_ids) = true;

        alpha_try(fail_ids) = 0.2 * alpha_try(fail_ids);

        if all(good_node | alpha_try < amin)
            break
        end
    end

    %If no node accepts the residual-decrease step, do not immediately fail.
    %In spinodal branch tracking, the regularized local tangent may not give
    %a strictly monotone residual every iteration. Accept a very small finite
    %bounded step as an inexact branch-tracking update.
    if ~any(good_node)

        alpha_small = min(alpha_now,0.03)*ones(1,N);
        accepted_small = false(1,N);

        for ksmall = 1:8

            c_small_raw = AddStep(c_old,dc_all,alpha_small);
            bound_small = BoundHitNode_Local(pars,c_small_raw,c_min,c_max,c_bound_margin);
            c_small     = BoundC_Local(c_small_raw,c_min,c_max);

            R_small = LE_ModelResidual_MaskFast(pars,p,c_small,E,eta,lam_c,true(1,N));
            accepted_small = isfinite(R_small) & ~bound_small & alpha_small >= amin;

            if any(accepted_small)
                break
            end

            alpha_small = 0.2*alpha_small;
        end

        if any(accepted_small)
            alpha_acc(accepted_small) = alpha_small(accepted_small);
            good_node(accepted_small) = true;
            ST.inexact = true;
        else
            failed = true;
            break
        end
    end

    %Update only accepted nodes. Unaccepted nodes get alpha = 0.
    c_new_raw = AddStep(c_old,dc_all,alpha_acc);
    bound_bad_acc = BoundHitNode_Local(pars,c_new_raw,c_min,c_max,c_bound_margin) & (alpha_acc > 0);

    if any(bound_bad_acc)
        alpha_acc(bound_bad_acc) = 0;
        good_node(bound_bad_acc) = false;
        ST.bound_hits = ST.bound_hits + nnz(bound_bad_acc);
    end

    c = AddStep(c_old,dc_all,alpha_acc);
    c = BoundC_Local(c,c_min,c_max);

    %Check whether large-step nodes failed the line search
    bad_large = any((~good_node) & (dcnode > c_tol));

    %Check convergence by accepted composition change
    cchg = 0;
    for ip = 1:length(c)
        if ~has_dof(ip)
            continue
        end
        for ic = 1:length(c{ip})
            cchg = max(cchg, max(abs(c{ip}{ic} - c_old{ip}{ic})));
        end
    end
    ST.max_cchg = cchg;

    if cchg < c_tol && ~bad_large
        converged = true;
        break
    end

    if it == Miter
        %Do not treat a bounded finite inexact branch update as a hard LE failure.
        %The caller can see DIAG.inexact/max_cchg, but STATE is not rolled back.
        ST.inexact = true;
        converged = true;
    end
end

if failed
    converged = false;
end

end

%This is the core quadratic solver for the LE
function [dc_all,mu_mat,chi_page] = LE_Quadratic_Step_Fast(pars,p,c,E,eta,lam_c)
%Prepare
Np         = length(c);
Ne         = length(E);
N          = numel(E{1});
E_mat      = cell2mat(E(:));
if isscalar(eta)
    eta_vec = eta * ones(1,N);
else
    eta_vec = eta(:).';
end

%Initialize step
for ip = 1:Np
    for ic = 1:length(c{ip})
        dc_all{ip}{ic} = zeros(size(c{ip}{ic})); %#ok<AGROW>
    end
end

%Prepare thermodynamic for each phase
R     = cell(1,Np);
Hreg  = cell(1,Np);
Bmix  = zeros(Ne, N);
Cmix  = zeros(Ne, Ne, N);

for ip = 1:Np

    R{ip} = PhaseThermo(pars{ip}, c{ip}, 'LE');
    e_ref = cell2mat(R{ip}.e(:));
    p_ip  = reshape(p(:,:,ip), 1, N);

    %Pure phase with no internal DOF
    if IsFrozenPhase(pars{ip}) || isempty(R{ip}.mu_c) || isempty(R{ip}.H_c) || isempty(R{ip}.Jac)
        B_phase  = e_ref;
        C_phase  = zeros(Ne, Ne, N);
    else
        %Normal phases
        mu_c     = cell2mat(R{ip}.mu_c(:));
        Hc       = R{ip}.H_c;
        J        = R{ip}.Jac;
        Nc       = size(mu_c,1);

        %Regularized positive definite Hessian, cached for the dc update
        Hreg{ip} = RegularizeHessian(Hc,lam_c);

        %H^{-1} * mu_c
        Hinv_mu3 = pagemldivide(Hreg{ip}, reshape(mu_c, Nc, 1, N));
        Hinv_mu  = reshape(Hinv_mu3, Nc, N);

        %H^{-1} * J^T
        JT       = permute(J, [2 1 3]);
        Hinv_JT  = pagemldivide(Hreg{ip}, JT);

        %B = e - J H^{-1} mu_c
        JHinvmu3 = pagemtimes(J, reshape(Hinv_mu, Nc, 1, N));
        JHinvmu  = reshape(JHinvmu3, Ne, N);
        B_phase  = e_ref - JHinvmu;

        %C = J H^{-1} J^T
        C_phase  = pagemtimes(J, Hinv_JT);
    end

    %Add Bmix and Cmix
    Bmix = Bmix + B_phase .* p_ip;
    Cmix = Cmix + C_phase .* reshape(p_ip, 1, 1, N);
end

%Calculate mu_e = (I/eta + Cmix)^(-1) * (E - Bmix)
IpagesE  = repmat(eye(Ne), 1, 1, N);
chi_page = Cmix + IpagesE .* reshape(1 ./ eta_vec, 1, 1, N);
rhs      = E_mat - Bmix;
mu3      = pagemldivide(chi_page, reshape(rhs, Ne, 1, N));
mu_mat   = reshape(mu3, Ne, N);

%Calculate dc = H^{-1}(J^T mu_e - mu_c)
for ip = 1:Np

    %If pure phase, no need to update c
    if IsFrozenPhase(pars{ip}) || isempty(R{ip}.mu_c) || isempty(R{ip}.H_c) || isempty(R{ip}.Jac)
        continue
    end

    mu_c    = cell2mat(R{ip}.mu_c(:));
    J       = R{ip}.Jac;
    Nc      = size(mu_c,1);
    JT      = permute(J, [2 1 3]);
    JTmu3   = pagemtimes(JT, reshape(mu_mat, Ne, 1, N));
    JTmu    = reshape(JTmu3, Nc, N);
    rhs_dc  = JTmu - mu_c;
    dc3     = pagemldivide(Hreg{ip}, reshape(rhs_dc, Nc, 1, N));
    dc      = reshape(dc3, Nc, N);

    for ic = 1:length(c{ip})
        dc_all{ip}{ic} = dc(ic,:);
    end
end
end


function Rn = LE_ModelResidual_MaskFast(pars,p,c,E,eta,lam_c,mask)
%LE_MODELRESIDUAL_MASKFAST Residual measured by next regularized LE correction.
%
%The residual is exactly the same type used by the branch-tracking line
%search in v2: the norm of the next regularized LE correction.  The speedup
%is that only the requested mask is sliced and evaluated, instead of
%recomputing PhaseThermo on all nodes for every line-search trial.

Nloc = nnz(mask);

if Nloc == 0
    Rn = zeros(1,0);
    return
end

try
    [pars_m,p_m,c_m,E_m,eta_m] = SliceLEBlock_Local(pars,p,c,E,eta,mask);
    [dc_next,~,~] = LE_Quadratic_Step_Fast(pars_m,p_m,c_m,E_m,eta_m,lam_c);
    Rn = MaxAbsStepNode(dc_next,Nloc);
    Rn(~isfinite(Rn)) = inf;
catch
    Rn = inf(1,Nloc);
end

end

function [pars_m,p_m,c_m,E_m,eta_m] = SliceLEBlock_Local(pars,p,c,E,eta,mask)

Np   = length(c);
Ne   = length(E);
Nloc = nnz(mask);

pars_m = cell(1,Np);
p_m    = zeros(1,Nloc,Np);
c_m    = cell(1,Np);
E_m    = cell(1,Ne);

for ip = 1:Np

    pars_m{ip} = SliceParsWScale_Local(pars{ip},mask);

    pp = reshape(p(:,:,ip),1,[]);
    p_m(1,:,ip) = pp(mask);

    c_m{ip} = cell(size(c{ip}));

    for ic = 1:numel(c{ip})
        A = reshape(c{ip}{ic},1,[]);
        c_m{ip}{ic} = A(mask);
    end
end

for ie = 1:Ne
    A = reshape(E{ie},1,[]);
    E_m{ie} = A(mask);
end

if isscalar(eta)
    eta_m = eta;
else
    eta_vec = reshape(eta,1,[]);
    eta_m = eta_vec(mask);
end

end

%This function regularize Hessian by adding a ridge and making it symmetric
function Hreg = RegularizeHessian(Hc,lam_c)

[Nc,~,N] = size(Hc);
Hs       = 0.5*(Hc + permute(Hc,[2 1 3]));

%Fast page-wise path. If pageeig is unavailable, use the loop fallback.
try
    [V,D] = pageeig(Hs);

    lam = zeros(Nc,N);
    for i = 1:Nc
        lam(i,:) = reshape(D(i,i,:),1,N);
    end

    scale = sqrt(reshape(sum(sum(Hs.^2,1),2),1,N))/max(1,Nc);
    scale = max(1,scale);

    evmin = min(lam,[],1);
    shift = max(0,1e-10*scale - evmin);

    lam_use = lam + shift + lam_c*scale;

    Duse = zeros(Nc,Nc,N);
    for i = 1:Nc
        Duse(i,i,:) = reshape(lam_use(i,:),1,1,N);
    end

    Hreg = pagemtimes(pagemtimes(V,Duse),permute(V,[2 1 3]));
    Hreg = 0.5*(Hreg + permute(Hreg,[2 1 3]));

catch

    Hreg = zeros(Nc,Nc,N);

    for i = 1:N
        H           = Hs(:,:,i);
        scale       = max(1,norm(H,'fro')/max(1,Nc));
        evmin       = min(eig(H));
        shift       = max(0,1e-10*scale - evmin);
        Hreg(:,:,i) = H + (shift + lam_c*scale) * eye(Nc);
    end
end
end

%This function bounds all independent composition variables.
%It is used only for safe evaluation after bound-hit nodes have been rejected.
function c = BoundC_Local(c,c_min,c_max)

for ip = 1:length(c)
    for ic = 1:length(c{ip})
        A = c{ip}{ic};
        A(~isfinite(A)) = 0;
        A = min(max(A,c_min),c_max);
        c{ip}{ic} = A;
    end
end

end

function hit = BoundHitNode_Local(pars,c,c_min,c_max,margin)

N   = numel(c{1}{1});
hit = false(1,N);

for ip = 1:length(c)

    % Pure/no-DOF phases may have c = 1 by construction.
    % Only phases with real composition DOF are tested.
    if ~HasDofFromPars(pars{ip})
        continue
    end

    for ic = 1:length(c{ip})

        A = reshape(c{ip}{ic},1,[]);

        hit = hit | ~isfinite(A);
        hit = hit | (A <= c_min + margin);
        hit = hit | (A >= c_max - margin);

    end
end

end

%This function simply update c with a step fraction alpha
function c_new = AddStep(c_old,dc_all,alpha_node)
c_new = c_old;
for ip = 1:length(c_old)
    for ic = 1:length(c_old{ip})
        c_new{ip}{ic} = c_old{ip}{ic} + alpha_node.*dc_all{ip}{ic};
    end
end
end

function dc_all = CapStep_Local(dc_all,dc_cap)

if isempty(dc_cap) || ~isfinite(dc_cap) || dc_cap <= 0
    return
end

for ip = 1:length(dc_all)
    for ic = 1:length(dc_all{ip})
        A = dc_all{ip}{ic};
        A = min(max(A,-dc_cap),dc_cap);
        dc_all{ip}{ic} = A;
    end
end

end

function dcmax = MaxAbsStep(dc_all)
dcmax = 0;
for ip = 1:length(dc_all)
    for ic = 1:length(dc_all{ip})
        dcmax = max(dcmax,max(abs(dc_all{ip}{ic})));
    end
end
end

function dcnode = MaxAbsStepNode(dc_all,N)

dcnode = zeros(1,N);

for ip = 1:length(dc_all)
    for ic = 1:length(dc_all{ip})
        dcnode = max(dcnode,abs(dc_all{ip}{ic}));
    end
end

end

function mu_e = CellFromMat_Local(mu_mat,Ne)
mu_e = cell(1,Ne);
for ie = 1:Ne
    mu_e{ie} = mu_mat(ie,:);
end
end

function chi = CellFromPage_Local(chi_page,Ne)
chi = cell(Ne,Ne);
for i = 1:Ne
    for j = 1:Ne
        chi{i,j} = reshape(chi_page(i,j,:),1,[]);
    end
end
end

function E0 = LE_MixtureElement_Local(pars,p,c,E)

Np = length(c);
Ne = length(E);
N  = numel(E{1});

Emix = zeros(Ne,N);

for ip = 1:Np

    try
        R = PhaseThermo(pars{ip},c{ip},'LE');
        e_ref = cell2mat(R.e(:));
    catch
        E0 = E;
        return
    end

    if any(~isfinite(e_ref(:)))
        E0 = E;
        return
    end

    p_ip  = reshape(p(:,:,ip),1,N);
    Emix  = Emix + e_ref .* p_ip;

end

E0 = CellFromMat_Local(Emix,Ne);

end

function Eout = BlendE_Local(E0,E1,theta)

Ne = length(E1);
Eout = cell(1,Ne);

for ie = 1:Ne
    Eout{ie} = (1-theta)*E0{ie} + theta*E1{ie};
end

end

function tf = HasDofFromPars(pars)

nReal = numel(pars.g0);

% A phase with only one real endmember is stoichiometric.
% Penalty/dependent endmembers should not create a physical DOF.
if nReal <= 1
    tf = false;
    return
end

tf = true;

end

function tf = IsFrozenPhase(pars)

nReal = numel(pars.g0);

% Fixed-composition stoichiometric phase.
% This includes one-real-endmember phases with penalty endmembers.
tf = nReal <= 1;

end

function pars_m = SliceParsWScale_Local(pars_m,mask)

if ~isfield(pars_m,'w_scale') || isempty(pars_m.w_scale)
    return
end

ws = reshape(pars_m.w_scale,1,[]);

if isscalar(ws)

    pars_m.w_scale = ws;

elseif numel(ws) == numel(mask)

    pars_m.w_scale = ws(mask);

elseif numel(ws) == nnz(mask)

    pars_m.w_scale = ws;

else

    error('SliceParsWScale_Local: w_scale length %d does not match mask length %d or nnz(mask) %d.', ...
          numel(ws),numel(mask),nnz(mask))

end

end
