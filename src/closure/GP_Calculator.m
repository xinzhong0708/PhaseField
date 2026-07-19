% =========================================================================
% GP calculator
% =========================================================================
function [c,chi,DIAG] = GP_Calculator(pars,p,c,mu_e,E_in,eta,alpha,Miter,PARAM)
%GP_CALCULATOR Fixed-mu GP projection with convex-split Hessian option.
%
% This function is intended for mixed-interface nodes in LE_Run_Mode when
% PARAM.LE_mode = 'GP'.
%
% Design:
%   - Raw PhaseThermo is used for e, mu_c, and objective evaluation.
%   - A convexified H_c can be used only as the numerical tangent in the GP
%     update, consistent with the convex-split CH/AC solver.
%   - Returned chi is convex-projected by default for stability.  Set
%     PARAM.GP_return_CS_chi = 0 to return the raw-projected chi.
%
% Backward compatibility:
%   The old call without PARAM still works, using conservative defaults.

if nargin < 9 || isempty(PARAM)
    PARAM = struct();
end

% ------------------------------------------------------------
% Parameters
% ------------------------------------------------------------
c_tol       = 1e-6;
wmu         = 1;
lam         = 1e-9;

MaxLS       = 10;
amin        = 1e-8;
obj_tol     = 1e-10;
obj_stop    = 1e-8;

dc_cap      = inf;
verbose     = 0;

%Artificial bounds for independent composition variables.
%Trials that hit/exceed these bounds are rejected instead of accepted.
c_min       = -1.0;
c_max       =  1.0;
c_bound_margin = 1e-12;

% Stable defaults. These are only used when the caller did not set them.
% Do not overwrite PARAM values coming from metadata or Run_2D.
if ~isfield(PARAM,'use_CS_chi'),       PARAM.use_CS_chi       = 1;       end
if ~isfield(PARAM,'GP_use_CS_H'),      PARAM.GP_use_CS_H      = 1;       end
if ~isfield(PARAM,'GP_return_CS_chi'), PARAM.GP_return_CS_chi = 1;       end
if ~isfield(PARAM,'GP_H_floor'),       PARAM.GP_H_floor       = 1e-10;   end
if ~isfield(PARAM,'GP_lam'),           PARAM.GP_lam           = 1e-7;    end
if ~isfield(PARAM,'GP_wmu'),           PARAM.GP_wmu           = 1;       end
if ~isfield(PARAM,'GP_c_tol'),         PARAM.GP_c_tol         = 1e-5;    end
if ~isfield(PARAM,'GP_dc_cap'),        PARAM.GP_dc_cap        = 0.005;   end
if ~isfield(PARAM,'GP_MaxLS'),         PARAM.GP_MaxLS         = 8;       end
if ~isfield(PARAM,'GP_obj_stop'),      PARAM.GP_obj_stop      = 1e-6;    end
if ~isfield(PARAM,'GP_verbose'),       PARAM.GP_verbose       = 0;       end

if isfield(PARAM,'GP_c_tol'),        c_tol    = PARAM.GP_c_tol;        end
if isfield(PARAM,'GP_wmu'),          wmu      = PARAM.GP_wmu;          end
if isfield(PARAM,'GP_lam'),          lam      = PARAM.GP_lam;          end
if isfield(PARAM,'GP_MaxLS'),        MaxLS    = PARAM.GP_MaxLS;        end
if isfield(PARAM,'GP_amin'),         amin     = PARAM.GP_amin;         end
if isfield(PARAM,'GP_obj_tol'),      obj_tol  = PARAM.GP_obj_tol;      end
if isfield(PARAM,'GP_obj_stop'),     obj_stop = PARAM.GP_obj_stop;     end
if isfield(PARAM,'GP_dc_cap'),       dc_cap   = PARAM.GP_dc_cap;       end
if isfield(PARAM,'GP_verbose'),      verbose  = PARAM.GP_verbose;      end
if isfield(PARAM,'GP_c_bound_margin')
    c_bound_margin = PARAM.GP_c_bound_margin;
elseif isfield(PARAM,'LE_c_bound_margin')
    c_bound_margin = PARAM.LE_c_bound_margin;
end

c_tol   = max(c_tol,0);
wmu     = max(wmu,0);
lam     = max(lam,0);
MaxLS   = max(round(MaxLS),1);
amin    = max(amin,0);
dc_cap  = max(dc_cap,0);

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
DIAG = struct();
DIAG.converged     = false;
DIAG.failed        = false;
DIAG.line_failed   = false;
DIAG.n_nodes       = numel(mu_e{1});
DIAG.alpha_stage   = alpha;
DIAG.n_iter        = 0;
DIAG.max_dc        = 0;
DIAG.max_cchg      = 0;
DIAG.bound_hits    = 0;
DIAG.iter          = 0;
DIAG.max_step      = inf;
DIAG.max_F_old     = inf;
DIAG.max_F_new     = inf;
DIAG.message       = '';
DIAG.use_CS_H      = Use_CS_H_Local(PARAM);
DIAG.return_CS_chi = isfield(PARAM,'GP_return_CS_chi') && PARAM.GP_return_CS_chi == 1;

converged   = false;
line_failed = false;
bound_hits  = 0;

% ------------------------------------------------------------
% Iteration
% ------------------------------------------------------------
for it = 1:Miter

    DIAG.iter   = it;
    DIAG.n_iter = it;

    % Current fixed-mu residual objective using raw thermodynamics.
    F_old = GP_Objective(pars,p,c,mu_e,E_in,eta,wmu);

    % Convex-split GP step.
    dc = GP_Step_MassConsistent(pars,p,c,mu_e,E_in,eta,lam,wmu,PARAM);

    % Optional hard safety cap on composition increment.
    if isfinite(dc_cap) && dc_cap > 0
        dc = CapStep(dc,dc_cap);
    end

    step_max = MaxStep(dc);
    DIAG.max_step = step_max;
    DIAG.max_dc   = max(DIAG.max_dc,step_max);

    finite_F = F_old(isfinite(F_old));

    if isempty(finite_F)
        DIAG.max_F_old = inf;
    else
        DIAG.max_F_old = max(finite_F);
    end

    % Convergence check before line search.
    if alpha*step_max < c_tol
        converged = true;
        break
    end

    % Backtracking line search, independently for each grid point.
    N          = numel(mu_e{1});
    good_node  = false(1,N);
    alpha_try  = alpha*ones(1,N);
    alpha_acc  = zeros(1,N);
    F_acc      = F_old;

    for ils = 1:MaxLS

        c_try_raw = AddStep(c,dc,alpha_try);
        bound_bad = BoundHitNode_Local(pars,c_try_raw,c_min,c_max,c_bound_margin);
        bound_hits = bound_hits + nnz((~good_node) & bound_bad);

        %Bound only for safe objective evaluation. Bound-hit nodes are
        %assigned F = inf below and therefore cannot be accepted.
        c_try = BoundC_Local(c_try_raw,c_min,c_max);
        F_try = GP_Objective(pars,p,c_try,mu_e,E_in,eta,wmu);
        F_try(bound_bad) = inf;

        good = isfinite(F_try) & (F_try <= F_old + obj_tol.*max(1,abs(F_old)));

        good_new            = good & ~good_node;
        alpha_acc(good_new) = alpha_try(good_new);
        F_acc(good_new)     = F_try(good_new);
        good_node(good_new) = true;

        bad = ~good_node;
        alpha_try(bad) = 0.2*alpha_try(bad);

        if all(good_node | alpha_try < amin)
            break
        end
    end

    % Accept only nodes with decreasing raw fixed-mu residual objective.
    if any(good_node)

        %Safety: reject accepted nodes if the final accepted trial still
        %hits/exceeds the artificial composition bounds.
        c_new_raw = AddStep(c,dc,alpha_acc);
        bound_bad_acc = BoundHitNode_Local(pars,c_new_raw,c_min,c_max,c_bound_margin) & (alpha_acc > 0);

        if any(bound_bad_acc)
            alpha_acc(bound_bad_acc) = 0;
            good_node(bound_bad_acc) = false;
            bound_hits = bound_hits + nnz(bound_bad_acc);
        end

        if any(alpha_acc > 0)
            c = AddStep(c,dc,alpha_acc);
            c = BoundC_Local(c,c_min,c_max);
        else
            line_failed = true;

            if verbose == 1
                disp('GP line search failed because all accepted nodes hit composition bounds...')
            end

            break
        end

    else
        line_failed = true;

        if verbose == 1
            disp('GP line search failed, keeping last accepted GP iterate...')
        end

        break
    end

    % Stop if accepted step is tiny.
    step_now = MaxStepScaled(dc,alpha_acc);
    DIAG.max_step = step_now;

    if step_now < c_tol
        converged = true;
        break
    end

    % Stop if objective no longer improves appreciably.
    good = good_node & isfinite(F_old) & isfinite(F_acc);

    if any(good)
        obj_change = max(abs(F_old(good)-F_acc(good))./max(1,abs(F_old(good))));
        DIAG.max_F_new = max(F_acc(good));

        if obj_change < obj_stop
            converged = true;
            break
        end
    end
end

if ~converged && ~line_failed && verbose == 1
    disp('GP warning: reached Miter, keeping last accepted GP iterate...')
end

DIAG.converged   = converged;
DIAG.line_failed = line_failed;
DIAG.failed      = line_failed;
DIAG.max_cchg    = DIAG.max_step;
DIAG.bound_hits  = bound_hits;

if line_failed
    if bound_hits > 0
        DIAG.message = sprintf('GP line search failed after rejecting %d bound-hit trial nodes',bound_hits);
    else
        DIAG.message = 'GP line search failed';
    end
elseif ~converged
    if bound_hits > 0
        DIAG.message = sprintf('GP reached Miter with safeguarded accepted iterates; rejected %d bound-hit trial nodes',bound_hits);
    else
        DIAG.message = 'GP reached Miter with safeguarded accepted iterates';
    end
elseif bound_hits > 0
    DIAG.message = sprintf('GP converged after rejecting %d bound-hit trial nodes',bound_hits);
end

% ------------------------------------------------------------
% Final susceptibility
% ------------------------------------------------------------
Ne  = numel(mu_e);
chi = Chi_FromMu_Projected(pars,p,c,eta,Ne,lam,wmu,PARAM);

end


% =========================================================================
% GP fixed-mu residual objective
% =========================================================================
function F = GP_Objective(pars,p,c,mu_e,E_in,eta,wmu)
%GP_OBJECTIVE Nodewise objective used only for line search.
%
% Raw thermodynamics is deliberately used here.
%
% Residuals:
%   mass residual:      E = Emix + mu/eta
%   chemical residual:  mu_c = J' * mu_e
%
% The line search accepts a step only when this raw residual objective
% decreases.

Np     = numel(c);
Ne     = numel(mu_e);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);
E_mat  = StackFields(E_in,N);
eta_v  = EtaVector(eta,N);

Emix   = zeros(Ne,N);
F      = zeros(1,N);

try
    Rall = cell(1,Np);

    for ip = 1:Np
        Rall{ip} = PhaseThermo(pars{ip},c{ip},'LE');
        p_ip     = reshape(p(:,:,ip),1,N);
        Emix     = Emix + StackFields(Rall{ip}.e,N).*p_ip;
    end

    R_E = E_mat - Emix - mu_mat./eta_v;
    F   = F + sum(R_E.^2,1);

    for ip = 1:Np

        R = Rall{ip};

        if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
            continue
        end

        J     = R.Jac;
        JT    = permute(J,[2 1 3]);
        Nc    = size(JT,1);
        mu_c  = StackFields(R.mu_c,N);
        p_ip  = reshape(p(:,:,ip),1,N);

        R_mu  = reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),Nc,N) - mu_c;

        F     = F + wmu*p_ip.*sum(R_mu.^2,1);
    end

    bad = ~isfinite(F);
    F(bad) = inf;

catch

    F = inf(1,N);

end

end


% =========================================================================
% GP step: fixed mu_e, raw residual, convex-split Hessian tangent
% =========================================================================
function dc = GP_Step_MassConsistent(pars,p,c,mu_e,E_in,eta,lam,wmu,PARAM)

Np     = numel(c);
Ne     = numel(mu_e);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);
E_mat  = StackFields(E_in,N);
eta_v  = EtaVector(eta,N);
dc     = ZeroStepLike(c);

Emix   = zeros(Ne,N);
Rall   = cell(1,Np);

% Current mixture composition.
for ip = 1:Np

    try
        Rall{ip} = PhaseThermo(pars{ip},c{ip},'LE');
    catch
        Rall{ip} = [];
        continue
    end

    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + StackFields(Rall{ip}.e,N).*p_ip;
end

% Mass residual: E = Emix + mu/eta.
R_E_global = E_mat - Emix - mu_mat./eta_v;

% Residual distribution weight.
psum2 = zeros(1,N);

for jp = 1:Np
    p_j   = reshape(p(:,:,jp),1,N);
    psum2 = psum2 + p_j.^2;
end

psum2 = max(psum2,eps);

% Phase-fraction weight for chemical projection.
p_freeze = 1e-4;
p_solve  = 1e-2;

if isfield(PARAM,'GP_p_freeze'), p_freeze = PARAM.GP_p_freeze; end
if isfield(PARAM,'GP_p_solve'),  p_solve  = PARAM.GP_p_solve;  end

if p_solve <= p_freeze
    p_solve = p_freeze + eps;
end

for ip = 1:Np

    R = Rall{ip};

    if isempty(R)
        continue
    end

    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    J     = R.Jac;
    Hraw  = SymPages(R.H_c);
    Huse  = ConvexifyH_ForGP(Hraw,PARAM);

    JT    = permute(J,[2 1 3]);
    HT    = permute(Huse,[2 1 3]);
    Nc    = size(Huse,1);
    mu_c  = StackFields(R.mu_c,N);

    % Chemical-potential residual uses raw mu_c.
    R_mu = reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),Nc,N) - mu_c;

    % This phase contributes to the mass residual in proportion to p_ip.
    p_ip = reshape(p(:,:,ip),1,N);
    R_E  = R_E_global.*p_ip./psum2;

    w  = (p_ip - p_freeze)./(p_solve - p_freeze);
    w  = min(max(w,0),1);
    w  = w.^2.*(3 - 2*w);
    w3 = reshape(w,1,1,N);

    A  = pagemtimes(JT,J) + wmu.*w3.*pagemtimes(HT,Huse) + ...
         repmat(eye(Nc),1,1,N).*lam;

    b  = pagemtimes(JT,reshape(R_E,Ne,1,N)) + ...
         wmu.*reshape(w,1,1,N).*pagemtimes(HT,reshape(R_mu,Nc,1,N));

    try
        d = reshape(pagemldivide(A,b),Nc,N);
    catch
        d = zeros(Nc,N);
    end

    d(~isfinite(d)) = 0;

    for ic = 1:Nc
        dc{ip}{ic} = d(ic,:);
    end
end

end


% =========================================================================
% GP susceptibility: chi = dE/dmu_e = sum_ip p_ip * de_ip/dmu_e + I/eta
% =========================================================================
function chi = Chi_FromMu_Projected(pars,p,c,eta,Ne,lam,wmu,PARAM)
% By default this returns the raw projected GP chi.
% Set PARAM.GP_return_CS_chi = 1 to return convex-H projected chi.

Np    = numel(c);
N     = numel(c{1}{1});
eta_v = EtaVector(eta,N);
Cmix  = zeros(Ne,Ne,N);

return_cs = isfield(PARAM,'GP_return_CS_chi') && PARAM.GP_return_CS_chi == 1;

for ip = 1:Np

    try
        R = PhaseThermo(pars{ip},c{ip},'LE');
    catch
        R = [];
    end

    p_ip = reshape(p(:,:,ip),1,N);

    if isempty(R) || isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)

        Cphase = zeros(Ne,Ne,N);

    else

        J      = R.Jac;
        Hraw   = SymPages(R.H_c);

        if return_cs == 1
            Huse = ConvexifyH_ForGP(Hraw,PARAM);
        else
            Huse = Hraw;
        end

        JT     = permute(J,[2 1 3]);
        HT     = permute(Huse,[2 1 3]);
        Nc     = size(Huse,1);

        A      = pagemtimes(JT,J) + wmu*pagemtimes(HT,Huse) + ...
                 repmat(eye(Nc),1,1,N).*lam;

        B      = wmu*pagemtimes(HT,JT);

        try
            Cphase = pagemtimes(J,pagemldivide(A,B));
        catch
            Cphase = zeros(Ne,Ne,N);
        end

        Cphase(~isfinite(Cphase)) = 0;

    end

    Cmix = Cmix + Cphase.*reshape(p_ip,1,1,N);
end

I   = repmat(eye(Ne),1,1,N);
chi = PagesToChi(Cmix + I.*reshape(1./eta_v,1,1,N),Ne);

end


% =========================================================================
% Convex-split Hessian for GP numerical tangent
% =========================================================================
function Huse = ConvexifyH_ForGP(H,PARAM)
%CONVEXIFYH_FORGP Positive Hessian used only for GP update/tangent.
%
% This mirrors the convex split idea used in the ACCH solver:
% raw thermodynamics is kept for the residual/objective, but the implicit
% tangent is made positive and bounded.

use_cs = Use_CS_H_Local(PARAM);

if use_cs ~= 1
    Huse = H;
    return
end

[Nc,~,N] = size(H);

h_floor = 1e-8;
h_cap   = inf;

if isfield(PARAM,'GP_H_floor')
    h_floor = PARAM.GP_H_floor;
elseif isfield(PARAM,'CS_chi_floor')
    h_floor = PARAM.CS_chi_floor;
end

if isfield(PARAM,'GP_H_cap')
    h_cap = PARAM.GP_H_cap;
end

h_floor = max(h_floor,0);

H = 0.5*(H + permute(H,[2 1 3]));

bad = squeeze(any(any(~isfinite(H),1),2));

if any(bad)
    H(:,:,bad) = repmat(eye(Nc),1,1,nnz(bad));
end

[V,D] = pageeig(H);

lam = zeros(Nc,N);

for i = 1:Nc
    lam(i,:) = reshape(D(i,i,:),1,N);
end

scale     = max(1,max(abs(lam),[],1));
floor_val = h_floor.*scale;

% Negative curvature is treated as positive implicit stiffness.
lam_use = max(abs(lam),floor_val);

if isfinite(h_cap)
    lam_use = min(lam_use,h_cap.*scale);
end

Duse = zeros(Nc,Nc,N);

for i = 1:Nc
    Duse(i,i,:) = reshape(lam_use(i,:),1,1,N);
end

Huse = pagemtimes(pagemtimes(V,Duse),permute(V,[2 1 3]));
Huse = 0.5*(Huse + permute(Huse,[2 1 3]));

end


function use_cs = Use_CS_H_Local(PARAM)

use_cs = 0;

if isfield(PARAM,'GP_use_CS_H')
    use_cs = PARAM.GP_use_CS_H;
elseif isfield(PARAM,'use_CS_chi')
    use_cs = PARAM.use_CS_chi;
end

end


% =========================================================================
% GP helpers
% =========================================================================
function H = SymPages(H)

H = 0.5*(H + permute(H,[2 1 3]));

end


function dc = ZeroStepLike(c)

dc = c;

for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        dc{ip}{ic} = zeros(size(dc{ip}{ic}));
    end
end

end


function dc = CapStep(dc,dc_cap)

for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        dc{ip}{ic} = max(min(dc{ip}{ic},dc_cap),-dc_cap);
    end
end

end


function m = MaxStep(dc)

m = 0;

for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        m = max(m,max(abs(dc{ip}{ic})));
    end
end

end


function m = MaxStepScaled(dc,alpha)

m = 0;

for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        m = max(m,max(abs(alpha.*dc{ip}{ic})));
    end
end

end


function c_new = AddStep(c,dc,alpha)

c_new = c;

for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        c_new{ip}{ic} = c{ip}{ic} + alpha.*dc{ip}{ic};
    end
end

end


function eta_v = EtaVector(eta,N)

if isscalar(eta)
    eta_v = eta*ones(1,N);
else
    eta_v = reshape(eta,1,N);
end

end


function chi = PagesToChi(A,Ne)

chi = cell(Ne,Ne);

for i = 1:Ne
    for j = 1:Ne
        chi{i,j} = reshape(A(i,j,:),1,[]);
    end
end

end


function A = StackFields(fields,N)

A = zeros(numel(fields),N);

for i = 1:numel(fields)
    A(i,:) = reshape(fields{i},1,N);
end

end


function c = BoundC_Local(c,c_min,c_max)

for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
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

for ip = 1:numel(c)

    % Pure/no-DOF phases may have c = 1 by construction.
    % Only phases with real composition DOF are tested.
    if ~HasCompositionDof_Local(pars{ip})
        continue
    end

    for ic = 1:numel(c{ip})

        A = reshape(c{ip}{ic},1,[]);

        hit = hit | ~isfinite(A);
        hit = hit | (A <= c_min + margin);
        hit = hit | (A >= c_max - margin);

    end
end

end


function tf = HasCompositionDof_Local(pars)

tf = true;

if isfield(pars,'g0')
    if numel(pars.g0) <= 1
        tf = false;
    end
end

end

