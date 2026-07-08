function [c,mu_e,chi,DIAG] = LE_Calculator_Stable(pars,p,c,E,eta,level,anchor)
%LE_CALCULATOR_STABLE Stable wrapper for LE_Calculator.
%
% Purpose:
%   Keep the existing fast LE_Calculator as the normal path, but when it
%   fails, retry through a small eta/prox continuation path from several
%   seeds.  This is a drop-in wrapper and does not assume phase names or
%   known composition domains.
%
% Usage:
%   [c,mu_e,chi,DIAG] = LE_Calculator_Stable(pars,p,c,E,eta,level)
%   [c,mu_e,chi,DIAG] = LE_Calculator_Stable(pars,p,c,E,eta,level,anchor)
%
% Optional anchor fields:
%   anchor.valid = true/false
%   anchor.c     = previous successful local LE composition, same format as c
%   anchor.mu_e  = previous successful local LE mu_e, same format as mu_e
%   anchor.chi   = previous successful local LE chi, same format as chi
%
% Notes:
%   1) This file calls your existing LE_Calculator.m internally.
%   2) LE_ContinuationTrust is included below as a local subfunction.
%   3) This is not a bound/domain hot fix.  It tries to avoid bad-basin
%      locking by changing the nonlinear path: weak eta/high proximal ->
%      full eta/normal proximal.

if nargin < 7 || isempty(anchor)
    anchor = struct();
end
if ~isfield(anchor,'valid') || isempty(anchor.valid)
    anchor.valid = false;
end

c_input = c;

DIAG = struct();
DIAG.failed       = false;
DIAG.soft_failed  = false;
DIAG.converged    = false;
DIAG.solver       = '';
DIAG.fast_diag    = [];
DIAG.cont_diag    = [];
DIAG.n_seeds      = 0;
DIAG.best_seed    = 0;
DIAG.message      = '';

% -------------------------------------------------------------------------
% 1. Fast path: original LE_Calculator
% -------------------------------------------------------------------------
fast_ok = false;
try
    [c_fast,mu_fast,chi_fast,diag_fast] = LE_Calculator(pars,p,c,E,eta,level);
    DIAG.fast_diag = diag_fast;

    if ~IsFailed_Local(diag_fast)
        c     = c_fast;
        mu_e  = mu_fast;
        chi   = chi_fast;
        DIAG  = MergeDiag_Local(DIAG,diag_fast);
        DIAG.failed    = false;
        DIAG.converged = true;
        DIAG.solver    = 'fast_LE';
        DIAG.message   = 'fast LE succeeded';
        return
    end

    % Even when failed, the current LE_Calculator normally returns a safe
    % response at its preserved input c.  Keep it as the no-anchor fallback.
    c_fallback    = c_fast;
    mu_fallback   = mu_fast;
    chi_fallback  = chi_fast;
    fast_ok       = true;

catch ME
    diag_fast = struct();
    diag_fast.failed = true;
    diag_fast.message = ME.message;
    DIAG.fast_diag = diag_fast;
    c_fallback = c_input;
    mu_fallback = [];
    chi_fallback = [];
end

% -------------------------------------------------------------------------
% 2. Build seed bank
% -------------------------------------------------------------------------
seeds = {};
seed_names = {};

seeds{end+1} = c_input;
seed_names{end+1} = 'current';

if IsValidAnchorC_Local(anchor,c_input)
    seeds{end+1} = anchor.c;
    seed_names{end+1} = 'anchor';
end

DIAG.n_seeds = numel(seeds);

% -------------------------------------------------------------------------
% 3. Continuation/trust fallback from each seed
% -------------------------------------------------------------------------
best = struct();
best.ok = false;
best.F  = inf;
best.c  = [];
best.mu = [];
best.chi = [];
best.diag = [];
best.seed = 0;

for iseed = 1:numel(seeds)

    c_seed = seeds{iseed};

    try
        [ct,mut,chit,dt] = LE_ContinuationTrust(pars,p,c_seed,E,eta,level);
    catch ME
        dt = struct();
        dt.failed = true;
        dt.message = ME.message;
        ct = c_seed;
        mut = [];
        chit = [];
    end

    if ~isfield(DIAG,'seed_diag') || isempty(DIAG.seed_diag)
        DIAG.seed_diag = cell(1,numel(seeds));
    end
    DIAG.seed_diag{iseed} = dt;

    if ~IsFailed_Local(dt)
        F = LE_Objective_All_Local(pars,p,ct,E,eta);
        Fmax = MaxFiniteOrInf_Local(F);

        if isfinite(Fmax) && Fmax < best.F
            best.ok   = true;
            best.F    = Fmax;
            best.c    = ct;
            best.mu   = mut;
            best.chi  = chit;
            best.diag = dt;
            best.seed = iseed;
        end
    end
end

if best.ok
    c     = best.c;
    mu_e  = best.mu;
    chi   = best.chi;

    DIAG = MergeDiag_Local(DIAG,best.diag);
    DIAG.failed    = false;
    DIAG.converged = true;
    DIAG.solver    = 'continuation_trust_LE';
    DIAG.best_seed = best.seed;
    DIAG.message   = ['continuation trust LE succeeded from seed ',seed_names{best.seed}];
    return
end

% -------------------------------------------------------------------------
% 4. If all retry paths fail, use last good anchor if it contains a full
% closure.  This is a soft fallback.  If no anchor exists, return original
% LE fallback and keep failed=true.
% -------------------------------------------------------------------------
if IsValidAnchorClosure_Local(anchor,c_input,E)
    c     = anchor.c;
    mu_e  = anchor.mu_e;
    chi   = anchor.chi;

    DIAG.failed      = false;
    DIAG.soft_failed = true;
    DIAG.converged   = false;
    DIAG.solver      = 'last_good_anchor';
    DIAG.message     = 'fast LE and continuation failed; used last good LE anchor';
    return
end

if fast_ok
    c     = c_fallback;
    mu_e  = mu_fallback;
    chi   = chi_fallback;
else
    error('LE_Calculator_Stable: fast LE failed before producing mu_e/chi and no valid anchor was supplied.')
end

DIAG.failed      = true;
DIAG.soft_failed = false;
DIAG.converged   = false;
DIAG.solver      = 'failed_return_fast_fallback';
DIAG.message     = 'fast LE and continuation failed; no valid anchor was available';

end


function [c,mu_e,chi,DIAG] = LE_ContinuationTrust(pars,p,c,E,eta,level)
%LE_CONTINUATIONTRUST Eta/prox continuation wrapper around LE_Calculator.
%
% This function deliberately keeps the implementation simple:
%   - weak eta first, so the phase is not forced to fit a difficult E all
%     at once;
%   - strong lam_c first, so the local quadratic update is more proximal;
%   - gradually return to full eta and the requested lam_c.
%
% It calls the existing LE_Calculator at each stage.

alpha0 = level(1);
Miter0 = level(2);

c_tol = 1e-5;
MaxLS = 10;
lam0  = 1e-7;

if numel(level) >= 3 && ~isempty(level(3)), c_tol = level(3); end
if numel(level) >= 4 && ~isempty(level(4)), MaxLS = level(4); end
if numel(level) >= 5 && ~isempty(level(5)), lam0  = level(5); end

eta_fac = [0.01 0.03 0.1 0.3 1.0];
lam_fac = [100 30 10 3 1];
alpha_fac = [0.25 0.35 0.5 0.75 1.0];

DIAG = struct();
DIAG.failed = false;
DIAG.converged = false;
DIAG.n_stage = numel(eta_fac);
DIAG.stage_failed = false(1,numel(eta_fac));
DIAG.stage_message = cell(1,numel(eta_fac));
DIAG.stage_diag = cell(1,numel(eta_fac));
DIAG.message = '';
DIAG.max_dc = 0;
DIAG.max_cchg = 0;
DIAG.n_iter = 0;
DIAG.alpha_stage = NaN;

for is = 1:numel(eta_fac)

    eta_s = ScaleEta_Local(eta,eta_fac(is));
    alpha_s = min(alpha0,alpha0*alpha_fac(is));
    if alpha_s <= 0
        alpha_s = alpha0;
    end

    Miter_s = max(Miter0,150);
    if is == numel(eta_fac)
        Miter_s = max(Miter_s,300);
    end

    lam_s = max(lam0*lam_fac(is),lam0);

    level_s = [alpha_s,Miter_s,c_tol,MaxLS,lam_s];

    %Try this stage. If it fails, retry the same continuation stage with
    %smaller alpha before giving up.
    ok_stage = false;
    retry_fac = [1 0.3 0.1 0.03];
    last_dt = [];
    last_msg = '';

    for ir = 1:numel(retry_fac)

        level_try = level_s;
        level_try(1) = max(alpha_s*retry_fac(ir),1e-4);

        try
            [c_try,mu_try,chi_try,dt] = LE_Calculator(pars,p,c,E,eta_s,level_try);
        catch ME
            dt = struct();
            dt.failed = true;
            dt.message = ME.message;
        end

        last_dt = dt;
        if isfield(dt,'message'), last_msg = dt.message; end

        if ~IsFailed_Local(dt)
            c = c_try;
            mu_e = mu_try;
            chi = chi_try;
            ok_stage = true;
            break
        end
    end

    DIAG.stage_diag{is} = last_dt;
    DIAG.stage_message{is} = last_msg;

    if isfield(last_dt,'max_dc') && isfinite(last_dt.max_dc)
        DIAG.max_dc = max(DIAG.max_dc,last_dt.max_dc);
    end
    if isfield(last_dt,'max_cchg') && isfinite(last_dt.max_cchg)
        DIAG.max_cchg = max(DIAG.max_cchg,last_dt.max_cchg);
    end
    if isfield(last_dt,'n_iter') && isfinite(last_dt.n_iter)
        DIAG.n_iter = DIAG.n_iter + last_dt.n_iter;
    end
    if isfield(last_dt,'alpha_stage') && isfinite(last_dt.alpha_stage)
        DIAG.alpha_stage = last_dt.alpha_stage;
    end

    if ~ok_stage
        DIAG.failed = true;
        DIAG.stage_failed(is) = true;
        DIAG.message = sprintf('continuation failed at stage %d, eta_fac = %.3e',is,eta_fac(is));
        return
    end
end

DIAG.failed = false;
DIAG.converged = true;
DIAG.message = 'eta/prox continuation succeeded';

end


function eta_s = ScaleEta_Local(eta,fac)
if isscalar(eta)
    eta_s = eta*fac;
else
    eta_s = eta.*fac;
end
end


function tf = IsFailed_Local(diag)
tf = false;
if isstruct(diag) && isfield(diag,'failed') && diag.failed
    tf = true;
end
end


function tf = IsValidAnchorC_Local(anchor,c_ref)
tf = false;
if ~isstruct(anchor) || ~isfield(anchor,'valid') || ~anchor.valid
    return
end
if ~isfield(anchor,'c') || isempty(anchor.c)
    return
end
if ~SameCSize_Local(anchor.c,c_ref)
    return
end
tf = true;
end


function tf = IsValidAnchorClosure_Local(anchor,c_ref,E_ref)
tf = false;
if ~IsValidAnchorC_Local(anchor,c_ref)
    return
end
if ~isfield(anchor,'mu_e') || isempty(anchor.mu_e)
    return
end
if ~isfield(anchor,'chi') || isempty(anchor.chi)
    return
end
if numel(anchor.mu_e) ~= numel(E_ref)
    return
end
tf = true;
end


function tf = SameCSize_Local(a,b)
tf = false;
if numel(a) ~= numel(b)
    return
end
for ip = 1:numel(a)
    if numel(a{ip}) ~= numel(b{ip})
        return
    end
    for ic = 1:numel(a{ip})
        if ~isequal(size(a{ip}{ic}),size(b{ip}{ic}))
            return
        end
    end
end
tf = true;
end


function F = LE_Objective_All_Local(pars,p,c,E,eta)
Np = numel(c);
Ne = numel(E);
N  = numel(E{1});

E_mat = cell2mat(E(:));
if isscalar(eta)
    eta_vec = eta*ones(1,N);
else
    eta_vec = reshape(eta,1,[]);
end

Gmix = zeros(1,N);
Emix = zeros(Ne,N);

try
    for ip = 1:Np
        R = PhaseThermo(pars{ip},c{ip},'LE');
        p_ip = reshape(p(:,:,ip),1,N);
        Gmix = Gmix + p_ip.*reshape(R.g,1,N);
        Emix = Emix + cell2mat(R.e(:)).*p_ip;
    end
    res = E_mat - Emix;
    F = Gmix + 0.5*eta_vec.*sum(res.^2,1);
    F(~isfinite(F)) = inf;
catch
    F = inf(1,N);
end
end


function m = MaxFiniteOrInf_Local(F)
F = F(:).';
F = F(isfinite(F));
if isempty(F)
    m = inf;
else
    m = max(F);
end
end


function D = MergeDiag_Local(D,src)
if ~isstruct(src)
    return
end
fn = fieldnames(src);
for i = 1:numel(fn)
    D.(fn{i}) = src.(fn{i});
end
end
