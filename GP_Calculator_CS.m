% =========================================================================
% GP calculator with residual-based convex-split fixed-mu mapping
% =========================================================================
function [c,chi] = GP_Calculator_CS(pars,p,c,mu_e,E_in,eta,alpha,Miter,PARAM)
%GP_CALCULATOR_CS Proximal fixed-mu local GP closure.
%
% This maps fixed mu_e to phase compositions by solving:
%
%     mu_c(c_new) + S*(c_new - c_ref) = J'*mu_e
%
% where c_ref is the composition at the beginning of this GP call.
%
% The line search minimizes the residual of this same equation:
%
%     R = J'*mu_e - mu_c(c) - S*(c-c_ref)
%
% This avoids the earlier problem where the Newton step and the merit
% function were not consistent.
%
% E_in is kept in the input list for compatibility but is not used.

if nargin < 9
    PARAM = struct();
end

% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------
c_tol        = 1e-6;
prox_base    = 1e-2;
prox_target  = 5e-2;
prox_posdef  = 1;
dc_cap       = 0.01;
move_cap     = inf;
c_floor      = 1e-9;

MaxLS        = 8;
amin         = 1e-7;
obj_tol      = 1e-8;

if isfield(PARAM,'LE_GP_c_tol')
    c_tol = PARAM.LE_GP_c_tol;
end
if isfield(PARAM,'LE_GP_prox_base')
    prox_base = PARAM.LE_GP_prox_base;
end
if isfield(PARAM,'LE_GP_prox_target')
    prox_target = PARAM.LE_GP_prox_target;
end
if isfield(PARAM,'LE_GP_prox_posdef')
    prox_posdef = PARAM.LE_GP_prox_posdef;
end
if isfield(PARAM,'LE_GP_dc_cap')
    dc_cap = PARAM.LE_GP_dc_cap;
end
if isfield(PARAM,'LE_GP_move_cap')
    move_cap = PARAM.LE_GP_move_cap;
end
if isfield(PARAM,'LE_GP_c_floor')
    c_floor = PARAM.LE_GP_c_floor;
end
if isfield(PARAM,'LE_GP_MaxLS')
    MaxLS = PARAM.LE_GP_MaxLS;
end
if isfield(PARAM,'LE_GP_obj_tol')
    obj_tol = PARAM.LE_GP_obj_tol;
end

% Reference branch at the beginning of this local GP call
c_ref = c;

% Fixed proximal shift for this GP call.
% This keeps the branch selector consistent during the line search.
S_ref = BuildProxShiftReference(pars,c_ref,prox_base,prox_target,prox_posdef);

% -------------------------------------------------------------------------
% Iteration
% -------------------------------------------------------------------------
for it = 1:Miter

    F_old = GP_Objective_Residual(pars,p,c,c_ref,mu_e,S_ref);

    dc = GP_ProxStep(pars,c,c_ref,mu_e,S_ref);

    dc = CapCellStep(dc,dc_cap);

    if alpha*MaxStep(dc) < c_tol
        break
    end

    % Nodewise residual line search.
    % Nodes that cannot find a decreasing residual simply keep alpha=0.
    N          = numel(mu_e{1});
    good_node  = false(1,N);
    alpha_try  = alpha*ones(1,N);
    alpha_acc  = zeros(1,N);

    for ils = 1:MaxLS

        active_ls = ~good_node & alpha_try >= amin;

        if ~any(active_ls)
            break
        end

        c_try = AddStepBounded(c,dc,alpha_try,c_floor,c_ref,move_cap);

        F_try = GP_Objective_Residual(pars,p,c_try,c_ref,mu_e,S_ref);

        good = active_ls & isfinite(F_try) & ...
               F_try <= F_old + obj_tol.*max(1,abs(F_old));

        alpha_acc(good) = alpha_try(good);
        good_node(good) = true;

        bad = active_ls & ~good;
        alpha_try(bad) = 0.2*alpha_try(bad);
    end

    % If no node can reduce the stabilized residual, keep last accepted c.
    if ~any(alpha_acc > 0)
        break
    end

    c_old = c;
    c     = AddStepBounded(c,dc,alpha_acc,c_floor,c_ref,move_cap);

    if MaxChange(c,c_old) < c_tol
        break
    end
end

% -------------------------------------------------------------------------
% Susceptibility from stabilized local map
% -------------------------------------------------------------------------
Ne  = numel(mu_e);
chi = Chi_FromMu_Prox(pars,p,c,eta,Ne,S_ref);

end


% =========================================================================
% Residual merit function for convex-split fixed-mu equation
% =========================================================================
function F = GP_Objective_Residual(pars,p,c,c_ref,mu_e,S_ref)

Np     = numel(c);
Ne     = numel(mu_e);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);

F = zeros(1,N);

try

    for ip = 1:Np

        R = PhaseThermo(pars{ip},c{ip});

        if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
            continue
        end

        J     = R.Jac;
        JT    = permute(J,[2 1 3]);
        Nc    = size(JT,1);

        mu_c  = StackFields(R.mu_c,N);
        c_now = StackFields(c{ip},N);
        c_old = StackFields(c_ref{ip},N);
        S     = S_ref{ip};

        % Residual of:
        %
        %     mu_c(c) + S*(c-c_ref) = J'*mu_e
        %
        Rmu = reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),Nc,N) - mu_c;
        Rmu = Rmu - c_now.*S + c_old.*S;

        p_ip = reshape(p(:,:,ip),1,N);

        F = F + p_ip.*sum(Rmu.^2,1);
    end

    bad = ~isfinite(F);
    F(bad) = inf;

catch

    F = inf(1,N);

end

end


% =========================================================================
% One proximal fixed-mu Newton step
% =========================================================================
function dc = GP_ProxStep(pars,c,c_ref,mu_e,S_ref)

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

    J     = R.Jac;
    H     = SymPages(R.H_c);
    JT    = permute(J,[2 1 3]);
    Nc    = size(H,1);

    mu_c  = StackFields(R.mu_c,N);
    c_now = StackFields(c{ip},N);
    c_old = StackFields(c_ref{ip},N);
    S     = S_ref{ip};

    % Raw fixed-mu residual:
    %
    %     J' mu_e - mu_c(c)
    %
    rhs = reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),Nc,N) - mu_c;

    % Convex-split/proximal residual:
    %
    %     J' mu_e - mu_c(c) - S*(c-c_ref)
    %
    rhs = rhs - c_now.*S + c_old.*S;

    % Linearized equation:
    %
    %     (H + S I) dc = rhs
    %
    Hs = AddDiagPages(H,S);

    d = reshape(pagemldivide(Hs,reshape(rhs,Nc,1,N)),Nc,N);

    for ic = 1:Nc
        dc{ip}{ic} = d(ic,:);
    end
end

end


% =========================================================================
% Susceptibility from proximal map
% =========================================================================
function chi = Chi_FromMu_Prox(pars,p,c,eta,Ne,S_ref)

Np    = numel(c);
N     = size(p,2);

if isscalar(eta)
    eta_v = eta*ones(1,N);
else
    eta_v = reshape(eta,1,N);
end

Cmix = zeros(Ne,Ne,N);

for ip = 1:Np

    R    = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);

    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)

        Cphase = zeros(Ne,Ne,N);

    else

        J  = R.Jac;
        H  = SymPages(R.H_c);
        JT = permute(J,[2 1 3]);
        S  = S_ref{ip};

        Hs = AddDiagPages(H,S);

        % From:
        %
        %     (H + S I) dc = J' dmu
        %
        % so:
        %
        %     de/dmu = J * inv(H + S I) * J'
        %
        Cphase = pagemtimes(J,pagemldivide(Hs,JT));
    end

    Cmix = Cmix + Cphase.*reshape(p_ip,1,1,N);
end

I   = repmat(eye(Ne),1,1,N);
chi = PagesToChi(Cmix + I.*reshape(1./eta_v,1,1,N),Ne);

end


% =========================================================================
% Build fixed proximal shift from reference branch
% =========================================================================
function S_ref = BuildProxShiftReference(pars,c_ref,prox_base,prox_target,prox_posdef)

Np    = numel(c_ref);
S_ref = cell(1,Np);

for ip = 1:Np

    R = PhaseThermo(pars{ip},c_ref{ip});

    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        S_ref{ip} = [];
        continue
    end

    H = SymPages(R.H_c);
    S_ref{ip} = ProxShiftPages(H,prox_base,prox_target,prox_posdef);

end

end


function S = ProxShiftPages(H,prox_base,prox_target,prox_posdef)

[Nc,~,N] = size(H);
S = zeros(1,N);

for i = 1:N

    A  = 0.5*(H(:,:,i) + H(:,:,i).');
    sc = max(1,norm(A,'fro')/max(1,Nc));

    % Always add a base proximal curvature
    s0 = prox_base*sc;

    if prox_posdef == 1

        evmin = min(eig(A));

        % Target a meaningful positive curvature, not barely positive.
        % H + S I should have min eigenvalue about prox_target*scale.
        s1 = prox_target*sc - evmin;

        S(i) = max(s0,s1);

    else

        S(i) = s0;

    end
end

S = max(S,0);

end


function Hs = AddDiagPages(H,S)

[Nc,~,N] = size(H);
Hs = H;

for i = 1:N
    Hs(:,:,i) = Hs(:,:,i) + S(i)*eye(Nc);
end

end


% =========================================================================
% Step control and bounds
% =========================================================================
function dc = CapCellStep(dc,dc_cap)

if dc_cap <= 0 || isinf(dc_cap)
    return
end

amp = zeros(size(dc{1}{1}));

for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        amp = max(amp,abs(dc{ip}{ic}));
    end
end

fac = ones(size(amp));
bad = amp > dc_cap;
fac(bad) = dc_cap./amp(bad);

for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        dc{ip}{ic} = dc{ip}{ic}.*fac;
    end
end

end


function c_new = AddStepBounded(c,dc,alpha,c_floor,c_ref,move_cap)

c_new = c;

for ip = 1:numel(c)

    Nc = numel(c{ip});
    N  = numel(c{ip}{1});
    X  = zeros(Nc,N);
    X0 = zeros(Nc,N);

    for ic = 1:Nc
        X(ic,:)  = c{ip}{ic} + alpha.*dc{ip}{ic};
        X0(ic,:) = c_ref{ip}{ic};
    end

    X(~isfinite(X)) = c_floor;
    X = max(X,c_floor);
    X = min(X,1-c_floor);

    % Optional general trust radius relative to the branch reference.
    % This is not phase-specific. Default is inf, so it is off.
    if isfinite(move_cap)
        dX = X - X0;
        dX = max(dX,-move_cap);
        dX = min(dX, move_cap);
        X  = X0 + dX;
    end

    % Keep endmember-like variables inside the simplex if needed.
    s = sum(X,1);
    bad = s > 1-c_floor;

    if any(bad)
        X(:,bad) = X(:,bad)./s(bad).*(1-c_floor);
    end

    for ic = 1:Nc
        c_new{ip}{ic} = reshape(X(ic,:),size(c{ip}{ic}));
    end
end

end


% =========================================================================
% Small helpers
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


function m = MaxStep(dc)

m = 0;

for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        m = max(m,max(abs(dc{ip}{ic})));
    end
end

end


function m = MaxChange(c,c_old)

m = 0;

for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        m = max(m,max(abs(c{ip}{ic} - c_old{ip}{ic})));
    end
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