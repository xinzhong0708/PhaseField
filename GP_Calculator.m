% =========================================================================
% GP calculator
% =========================================================================
function [c,chi] = GP_Calculator(pars,p,c,mu_e,E_in,eta,alpha,Miter)
c_tol       = 1e-8;
wmu         = 1;
lam         = 1e-9;

%Line-search safeguard
MaxLS       = 8;
amin        = 1e-7;
obj_tol     = 1e-8;
obj_stop    = 1e-8;

converged   = false;
line_failed = false;

for it = 1:Miter

    %Current fixed-mu residual objective
    F_old = GP_Objective(pars,p,c,mu_e,E_in,eta,wmu);

    %Original GP step
    dc = GP_Step_MassConsistent(pars,p,c,mu_e,E_in,eta,lam,wmu);

    %Convergence check before line search
    if alpha*MaxStep(dc) < c_tol
        converged = true;
        break
    end

    %Backtracking line search, independently for each grid point
    N          = numel(mu_e{1});
    good_node  = false(1,N);
    alpha_try  = alpha*ones(1,N);
    alpha_acc  = zeros(1,N);
    F_acc      = F_old;

    for ils = 1:MaxLS

        c_try = AddStep(c,dc,alpha_try);
        F_try = GP_Objective(pars,p,c_try,mu_e,E_in,eta,wmu);

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

    %Accept only nodes with decreasing fixed-mu residual objective
    if any(good_node)
        c = AddStep(c,dc,alpha_acc);
    else
        line_failed = true;
        disp('GP line search failed, keeping last accepted GP iterate...')
        break
    end

    %Stop if accepted step is tiny
    step_now = MaxStepScaled(dc,alpha_acc);
    if step_now < c_tol
        converged = true;
        break
    end

    %Stop if objective no longer improves appreciably.
    good = good_node & isfinite(F_old) & isfinite(F_acc);
    if any(good)
        obj_change = max(abs(F_old(good)-F_acc(good))./max(1,abs(F_old(good))));
        if obj_change < obj_stop
            converged = true;
            break
        end
    end
end

%This is now a warning, not a failed returned state.
%With line search, reaching Miter means no step-size convergence was reached
%within the iteration budget, but all accepted steps were safeguarded.
if ~converged && ~line_failed
    disp('GP warning: reached Miter, keeping last accepted GP iterate...')
end

Ne  = numel(mu_e);
chi = Chi_FromMu_Projected(pars,p,c,eta,Ne,lam,wmu);
end


% =========================================================================
% GP fixed-mu residual objective
% =========================================================================
function F = GP_Objective(pars,p,c,mu_e,E_in,eta,wmu)
%GP_OBJECTIVE Nodewise objective used only for line search.
%
% This is not a new thermodynamic model. It measures whether the current
% fixed-mu GP residual is reduced:
%
%   mass residual:      E = Emix + mu/eta
%   chemical residual:  mu_c = J' * mu_e
%
% The line search accepts a step only when this residual objective decreases.

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
        Rall{ip} = PhaseThermo(pars{ip},c{ip});
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
% GP step: fixed mu_e
% =========================================================================
function dc = GP_Step_MassConsistent(pars,p,c,mu_e,E_in,eta,lam,wmu)
Np     = numel(c);
Ne     = numel(mu_e);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);
E_mat  = StackFields(E_in,N);
eta_v  = EtaVector(eta,N);
dc     = ZeroStepLike(c);
Emix   = zeros(Ne,N);
Rall   = cell(1,Np);
%Current mixture composition
for ip = 1:Np
    Rall{ip} = PhaseThermo(pars{ip},c{ip});
    p_ip     = reshape(p(:,:,ip),1,N);
    Emix     = Emix + StackFields(Rall{ip}.e,N).*p_ip;
end
%Mass residual: E = Emix + mu/eta
R_E_global = E_mat - Emix - mu_mat./eta_v;

%Residual distribution weight
psum2 = zeros(1,N);
for jp = 1:Np
    p_j   = reshape(p(:,:,jp),1,N);
    psum2 = psum2 + p_j.^2;
end
psum2 = max(psum2,eps);

for ip = 1:Np
    R = Rall{ip};
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end
    J    = R.Jac;
    H    = SymPages(R.H_c);
    JT   = permute(J,[2 1 3]);
    HT   = permute(H,[2 1 3]);
    Nc   = size(H,1);
    mu_c = StackFields(R.mu_c,N);
    %Chemical-potential residual: %   mu_c = J' * mu_e
    R_mu = reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),Nc,N) - mu_c;
    %This phase contributes to the mass residual in proportion to p_ip
    p_ip = reshape(p(:,:,ip),1,N);
    R_E  = R_E_global.*p_ip./psum2;
    % %Solve mixed local projection: J dc  ~= R_E  and   H dc  ~= R_mu
    % A = pagemtimes(JT,J) + wmu*pagemtimes(HT,H) + repmat(eye(Nc),1,1,N).*lam;
    % b = pagemtimes(JT,reshape(R_E,Ne,1,N)) + wmu*pagemtimes(HT,reshape(R_mu,Nc,1,N));

    % Phase-fraction weight for chemical projection
    p_freeze = 1e-4;
    p_solve  = 1e-2;
    w  = (p_ip - p_freeze)./(p_solve - p_freeze);
    w  = min(max(w,0),1);
    w  = w.^2.*(3 - 2*w);
    w3 = reshape(w,1,1,N);
    A  = pagemtimes(JT,J) + wmu.*w3.*pagemtimes(HT,H) + repmat(eye(Nc),1,1,N).*lam;
    b  = pagemtimes(JT,reshape(R_E,Ne,1,N)) + wmu.*reshape(w,1,1,N).*pagemtimes(HT,reshape(R_mu,Nc,1,N));
    d  = reshape(pagemldivide(A,b),Nc,N);
    for ic = 1:Nc
        dc{ip}{ic} = d(ic,:);
    end
end
end


% =========================================================================
% GP susceptibility: chi = dE/dmu_e = sum_ip p_ip * de_ip/dmu_e + I/eta
% =========================================================================
function chi = Chi_FromMu_Projected(pars,p,c,eta,Ne,lam,wmu)
Np    = numel(c);
N     = numel(c{1}{1});
eta_v = EtaVector(eta,N);
Cmix  = zeros(Ne,Ne,N);
for ip = 1:Np
    R    = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        Cphase = zeros(Ne,Ne,N);
    else
        J      = R.Jac;
        H      = SymPages(R.H_c);
        JT     = permute(J,[2 1 3]);
        HT     = permute(H,[2 1 3]);
        Nc     = size(H,1);
        A      = pagemtimes(JT,J) + wmu*pagemtimes(HT,H) + repmat(eye(Nc),1,1,N).*lam;
        B      = wmu*pagemtimes(HT,JT);
        Cphase = pagemtimes(J,pagemldivide(A,B));
    end
    Cmix = Cmix + Cphase.*reshape(p_ip,1,1,N);
end
I   = repmat(eye(Ne),1,1,N);
chi = PagesToChi(Cmix + I.*reshape(1./eta_v,1,1,N),Ne);
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
