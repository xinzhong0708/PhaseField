function [c,mu_e,chi,DIAG] = LE_Calculator(pars,p,c,E,eta,level)
%LE_CALCULATOR Vectorized documented quadratic local-equilibrium solver.
%
% This version follows the equations in the LE document directly:
%
%   chi = I/eta + sum_ip p_ip * J_ip * H_ip^{-1} * J_ip'
%
%   chi * mu = E - sum_ip p_ip * (e_ip - J_ip*H_ip^{-1}*mu_c_ip)
%
%   dc_ip = H_ip^{-1} * (J_ip' * mu - mu_c_ip)
%
% Main design choices:
%   1. No positive Hessian regularization.  Negative H eigenvalues are kept,
%      so spinodal negative chi can be returned.
%   2. No nodewise line search.  Only one outer iteration loop is used.
%   3. All H^{-1} operations use pagemldivide/pagemtimes.
%   4. The update is simply c <- c + alpha*dc.
%
% level:
%   level(1) = alpha, damping factor, e.g. 0.1 to 0.5
%   level(2) = max iterations
%   level(3) = optional tiny diagonal shift added to H before solve.
%              Default 0. Use only if H is exactly singular.
%
% Recommended first spinodal test:
%   [0.2, 100, 0]
%
% If the update oscillates:
%   [0.05, 200, 0]
%
% If exact singularity gives NaN:
%   [0.05, 200, 1e-12]
%
% The optional shift is not a positive-definite convexification. It is only
% H <- H + h_shift*I.  If h_shift is too large, it changes the physics.

% -------------------------------------------------------------------------
% Sizes and controls
% -------------------------------------------------------------------------
Np    = numel(c);
Ne    = numel(E);
N     = numel(E{1});

alpha = level(1);
Miter = level(2);

if numel(level) >= 3 && ~isempty(level(3))
    h_shift = level(3);
else
    h_shift = 0;
end

tol_dc  = 1e-10;
tol_res = 1e-8;

if isscalar(eta)
    eta_vec = eta*ones(1,N);
else
    eta_vec = eta(:).';
end

E_mat = cell2mat(E(:));

% -------------------------------------------------------------------------
% Main fixed-point loop
% -------------------------------------------------------------------------
DIAG.res_hist = nan(1,Miter);
DIAG.dc_hist  = nan(1,Miter);
DIAG.mu_hist  = nan(1,Miter);

mu_old = zeros(Ne,N);

for it = 1:Miter

    [dc,mu_mat,~,STEP] = QuadraticStep_Documented(pars,p,c,E_mat,eta_vec,h_shift);

    dcmax              = MaxAbsCell(dc);
    res                = ResidualNorm_Documented(pars,p,c,E_mat,eta_vec,mu_mat);

    DIAG.dc_hist(it)   = dcmax;
    DIAG.res_hist(it)  = max(res);
    DIAG.mu_hist(it)   = max(abs(mu_mat(:)-mu_old(:)));

    if dcmax < tol_dc || max(res) < tol_res
        break
    end

    % One simple damped update, no inner line search.
    c = AddCellStep(c,dc,alpha);

    mu_old = mu_mat;

    if STEP.has_nonfinite
        break
    end
end

if it == Miter
    disp('LE fails')
end

% -------------------------------------------------------------------------
% Final output at accepted c
% -------------------------------------------------------------------------
[~,mu_mat,chi_page,STEP] = QuadraticStep_Documented( ...
    pars,p,c,E_mat,eta_vec,h_shift);

mu_e = cell(1,Ne);
for ie = 1:Ne
    mu_e{ie} = reshape(mu_mat(ie,:),size(E{ie}));
end

chi = cell(Ne,Ne);
for ie = 1:Ne
    for je = 1:Ne
        chi{ie,je} = reshape(chi_page(ie,je,:),size(E{1}));
    end
end

DIAG.iter              = it;
DIAG.max_residual      = max(ResidualNorm_Documented(pars,p,c,E_mat,eta_vec,mu_mat));
DIAG.max_dc            = MaxAbsCell(dc);
DIAG.has_nonfinite     = STEP.has_nonfinite;
DIAG.h_shift           = h_shift;

Cdiag = ChiDiagnostics(chi_page);
DIAG.min_chi_eig        = Cdiag.min_chi_eig;
DIAG.negative_chi_nodes = Cdiag.negative_chi_nodes;
DIAG.max_abs_chi        = Cdiag.max_abs_chi;
DIAG.has_nonfinite_chi  = Cdiag.has_nonfinite_chi;

DIAG.res_hist = DIAG.res_hist(1:it);
DIAG.dc_hist  = DIAG.dc_hist(1:it);
DIAG.mu_hist  = DIAG.mu_hist(1:it);

end


% =========================================================================
% Vectorized documented quadratic step
% =========================================================================
function [dc_all,mu_mat,chi_page,DIAG] = QuadraticStep_Documented(pars,p,c,E_mat,eta_vec,h_shift)

Np = numel(c);
Ne = size(E_mat,1);
N  = size(E_mat,2);

Bmix = zeros(Ne,N);
Cmix = zeros(Ne,Ne,N);

dc_all = c;
for ip = 1:Np
    for ic = 1:numel(c{ip})
        dc_all{ip}{ic} = zeros(size(c{ip}{ic}));
    end
end

R = cell(1,Np);

DIAG.has_nonfinite = false;

% -------------------------------------------------------------------------
% First pass: build Bmix and Cmix
% -------------------------------------------------------------------------
for ip = 1:Np

    R{ip} = PhaseThermo(pars{ip},c{ip});

    e_ref = cell2mat(R{ip}.e(:));
    p_ip  = reshape(p(:,:,ip),1,N);

    if isempty(R{ip}.mu_c) || isempty(R{ip}.H_c) || isempty(R{ip}.Jac)

        B_phase = e_ref;
        C_phase = zeros(Ne,Ne,N);

    else

        mu_c = cell2mat(R{ip}.mu_c(:));
        H    = SymPages(R{ip}.H_c);
        J    = R{ip}.Jac;
        Nc   = size(mu_c,1);

        if h_shift ~= 0
            H = H + repmat(eye(Nc),1,1,N).*h_shift;
        end

        % H^{-1} * mu_c
        Hinv_mu3 = pagemldivide(H,reshape(mu_c,Nc,1,N));
        Hinv_mu  = reshape(Hinv_mu3,Nc,N);

        % H^{-1} * J'
        JT       = permute(J,[2 1 3]);
        Hinv_JT  = pagemldivide(H,JT);

        % J*H^{-1}*mu_c
        JHinvmu3 = pagemtimes(J,reshape(Hinv_mu,Nc,1,N));
        JHinvmu  = reshape(JHinvmu3,Ne,N);

        % Documented shorthand:
        % B_i = e_i - J H^{-1} mu_c
        % C_i = J H^{-1} J'
        B_phase = e_ref - JHinvmu;
        C_phase = pagemtimes(J,Hinv_JT);

        if any(~isfinite(B_phase(:))) || any(~isfinite(C_phase(:)))
            DIAG.has_nonfinite = true;
        end
    end

    Bmix = Bmix + B_phase.*p_ip;
    Cmix = Cmix + C_phase.*reshape(p_ip,1,1,N);
end

% -------------------------------------------------------------------------
% Solve common mu:
%   (I/eta + sum p*C_i) mu = E - sum p*B_i
% -------------------------------------------------------------------------
Ipages   = repmat(eye(Ne),1,1,N);
chi_page = Cmix + Ipages.*reshape(1./eta_vec,1,1,N);
rhs      = E_mat - Bmix;

mu3 = pagemldivide(SymPages(chi_page),reshape(rhs,Ne,1,N));
mu_mat = reshape(mu3,Ne,N);

if any(~isfinite(mu_mat(:)))
    DIAG.has_nonfinite = true;
end

% -------------------------------------------------------------------------
% Second pass: dc_i = H^{-1}(J' mu - mu_c)
% -------------------------------------------------------------------------
for ip = 1:Np

    if isempty(R{ip}.mu_c) || isempty(R{ip}.H_c) || isempty(R{ip}.Jac)
        continue
    end

    mu_c = cell2mat(R{ip}.mu_c(:));
    H    = SymPages(R{ip}.H_c);
    J    = R{ip}.Jac;
    Nc   = size(mu_c,1);

    if h_shift ~= 0
        H = H + repmat(eye(Nc),1,1,N).*h_shift;
    end

    JTmu3 = pagemtimes(permute(J,[2 1 3]),reshape(mu_mat,Ne,1,N));
    JTmu  = reshape(JTmu3,Nc,N);

    rhs_dc = JTmu - mu_c;

    dc3 = pagemldivide(H,reshape(rhs_dc,Nc,1,N));
    dc  = reshape(dc3,Nc,N);

    if any(~isfinite(dc(:)))
        DIAG.has_nonfinite = true;
    end

    for ic = 1:numel(c{ip})
        dc_all{ip}{ic} = reshape(dc(ic,:),size(c{ip}{ic}));
    end
end

end


% =========================================================================
% Residual diagnostics
% =========================================================================
function res = ResidualNorm_Documented(pars,p,c,E_mat,eta_vec,mu_mat)

Np = numel(c);
Ne = size(E_mat,1);
N  = size(E_mat,2);

Emix = zeros(Ne,N);
res  = zeros(1,N);

for ip = 1:Np

    R = PhaseThermo(pars{ip},c{ip});
    e = cell2mat(R.e(:));
    p_ip = reshape(p(:,:,ip),1,N);

    Emix = Emix + e.*p_ip;

    if ~isempty(R.mu_c) && ~isempty(R.Jac)
        mu_c = cell2mat(R.mu_c(:));
        JTmu3 = pagemtimes(permute(R.Jac,[2 1 3]),reshape(mu_mat,Ne,1,N));
        JTmu  = reshape(JTmu3,size(mu_c,1),N);
        res   = max(res,max(abs(mu_c-JTmu),[],1));
    end
end

mu_penalty = (E_mat-Emix).*reshape(eta_vec,1,N);
res = max(res,max(abs(mu_mat-mu_penalty),[],1));

bad = any(~isfinite(mu_mat),1) | any(~isfinite(Emix),1);
res(bad) = inf;

end


% =========================================================================
% Small utilities
% =========================================================================
function Hs = SymPages(H)
Hs = 0.5*(H+permute(H,[2 1 3]));
end


function c = AddCellStep(c,dc,alpha)

for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        step = dc{ip}{ic};
        good = isfinite(step);
        tmp = c{ip}{ic};
        tmp(good) = tmp(good) + alpha.*step(good);
        c{ip}{ic} = tmp;
    end
end

end


function m = MaxAbsCell(C)

m = 0;

for ip = 1:numel(C)
    for ic = 1:numel(C{ip})
        v = abs(C{ip}{ic});
        v = v(isfinite(v));
        if ~isempty(v)
            m = max(m,max(v));
        end
    end
end

end


function D = ChiDiagnostics(chi_page)

[~,~,N] = size(chi_page);
mineig = nan(1,N);
maxabs = 0;
bad = false;

for inode = 1:N
    A = 0.5*(chi_page(:,:,inode)+chi_page(:,:,inode).');

    if any(~isfinite(A(:)))
        bad = true;
        continue
    end

    ev = eig(A);
    mineig(inode) = min(ev);
    maxabs = max(maxabs,max(abs(A(:))));
end

valid = mineig(isfinite(mineig));
if isempty(valid)
    D.min_chi_eig = nan;
else
    D.min_chi_eig = min(valid);
end

D.negative_chi_nodes = nnz(mineig < 0);
D.max_abs_chi        = maxabs;
D.has_nonfinite_chi  = bad;

end
