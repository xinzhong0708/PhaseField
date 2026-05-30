function [c,mu_e,chi,DIAG] = LE_Calculator(pars,p,c,E,eta,level)
%LE_CALCULATOR Spinodal-preserving projected local-equilibrium calculator.
%
% Drop-in replacement for the LE_Calculator used by LE_Run.
% External input/output structure is unchanged:
%   [c,mu_e,chi] = LE_Calculator(pars,p,c,E,eta,level)
%
% Main change relative to PFM_Update3:
%   - The thermodynamic response used in the common-mu projection is the
%     physical R.chi returned by PhaseThermo. It is NOT made positive
%     definite. Therefore negative susceptibility in a spinodal region is
%     retained and returned to the CH solver.
%   - Composition updates are performed by remapping target elemental
%     composition through R.Jac (de/dc), not by using a convexified H_c.
%   - Globalization is based on reduction of the LE residual, not decrease
%     of local Gibbs energy. A spinodal state is not a local Gibbs minimum,
%     so an energy-decrease line search would suppress the instability.
%
% Equations used at each projection iteration:
%   S_ip       = d e_ip / d mu = R_ip.chi          (may be indefinite)
%   chi_eff    = I/eta + sum_ip p_ip * S_ip        (may be indefinite)
%   e_ip(mu)   = e_ip^0 + S_ip*(mu - mu_ip^0)
%   E          = mu/eta + sum_ip p_ip*e_ip(mu)
%
% The last equation is solved for the common mu, and each target e_ip is
% remapped back to c_ip through the kinematic Jacobian J_ip = de_ip/dc_ip.
%
% IMPORTANT:
%   This allows spinodal chi to reach PF_Coupled_ACCH_LETangent. The CH
%   step must include nonzero kappa and may need a restricted timestep.

% -------------------------------------------------------------------------
% Settings and sizes
% -------------------------------------------------------------------------
Np        = numel(c);
Ne        = numel(E);
N         = numel(E{1});
alpha     = level(1);
Miter     = level(2);
c_tol     = 1e-7;
res_tol   = 1e-7;
MaxLS     = 10;
amin      = 1e-8;
accept_eps= 1e-12;

DIAG.converged     = false;
DIAG.iter          = 0;
DIAG.max_res       = inf;
DIAG.min_chi_eig   = nan;
DIAG.failed_nodes  = [];

% -------------------------------------------------------------------------
% Projection / Newton iteration
% -------------------------------------------------------------------------
for it = 1:Miter

    c_old = c;

    % Physical signed local tangent and projected common chemical potential
    [e_tar,mu_old,chi_old,R_old] = LE_Physical_Projection(pars,p,c_old,E,eta);
    res_old = LE_Residual(R_old,p,E,eta,mu_old);

    % Map target elemental compositions to composition increments.
    % This uses only de/dc and does not alter spinodal curvature.
    dc_all = Map_TargetE_To_c_Step(c_old,R_old,e_tar);
    dcmax  = MaxAbsStep(dc_all);

    if max(res_old) < res_tol || dcmax < c_tol
        DIAG.converged = max(res_old) < res_tol;
        DIAG.iter      = it - 1;
        break
    end

    % Residual-decreasing, nodewise line search. Do NOT minimize G here:
    % spinodal compositions are allowed to be locally unstable.
    good_node = false(1,N);
    alpha_try = alpha * ones(1,N);
    alpha_acc = zeros(1,N);

    for ils = 1:MaxLS
        c_try = AddStep(c_old,dc_all,alpha_try);

        try
            [~,mu_try,~,R_try] = LE_Physical_Projection(pars,p,c_try,E,eta);
            res_try = LE_Residual(R_try,p,E,eta,mu_try);
            good = isfinite(res_try) & ...
                   (res_try <= res_old .* (1 + 1e-10) + accept_eps);
        catch
            good = false(1,N);
        end

        good_new            = good & ~good_node;
        alpha_acc(good_new) = alpha_try(good_new);
        good_node(good_new) = true;

        bad = ~good_node;
        alpha_try(bad) = 0.2 * alpha_try(bad);

        if all(good_node | alpha_try < amin)
            break
        end
    end

    if ~any(good_node)
        % Return the last finite iterate, with diagnostics instead of
        % silently reporting false convergence.
        DIAG.iter         = it;
        DIAG.max_res      = max(res_old);
        DIAG.failed_nodes = find(~good_node);
        break
    end

    c = AddStep(c_old,dc_all,alpha_acc);
    DIAG.iter = it;

    if it == Miter
        % Retain previous fallback style, but remain in signed-chi mode.
        if alpha > 0.05
            level(1) = 0.05;
            [c,mu_e,chi,DIAG] = LE_Calculator(pars,p,c,E,eta,level);
            return
        end
    end
end

% -------------------------------------------------------------------------
% Final physical signed response returned to CH
% -------------------------------------------------------------------------
[~,mu_mat,chi_page,R_final] = LE_Physical_Projection(pars,p,c,E,eta);
res_final = LE_Residual(R_final,p,E,eta,mu_mat);

for ie = 1:Ne
    mu_e{ie} = mu_mat(ie,:);
end
for i = 1:Ne
    for j = 1:Ne
        chi{i,j} = reshape(chi_page(i,j,:),1,[]);
    end
end

DIAG.max_res      = max(res_final);
DIAG.converged    = DIAG.max_res < res_tol;
DIAG.failed_nodes = find(res_final >= res_tol);
DIAG.min_chi_eig  = MinPageEigenvalue(chi_page);

end


% =========================================================================
% Physical signed LE projection: this is the key spinodal-preserving part.
% =========================================================================
function [e_tar,mu_mat,chi_page,R] = LE_Physical_Projection(pars,p,c,E,eta)

Np    = numel(c);
Ne    = numel(E);
N     = numel(E{1});
E_mat = cell2mat(E(:));

if isscalar(eta)
    eta_vec = eta * ones(1,N);
else
    eta_vec = eta(:).';
end

Bmix  = zeros(Ne,N);
Cmix  = zeros(Ne,Ne,N);
e_tar = cell(Np,1);
R     = cell(Np,1);
S     = cell(Np,1);
mu_ref= cell(Np,1);
e_ref = cell(Np,1);

for ip = 1:Np
    R{ip}   = PhaseThermo(pars{ip},c{ip});
    e_ref{ip} = cell2mat(R{ip}.e(:));
    p_ip    = reshape(p(:,:,ip),1,N);

    if isempty(R{ip}.mu_e) || isempty(R{ip}.chi) || isempty(R{ip}.Jac)
        % Pure/no-DOF phase: fixed elemental composition, zero response.
        S{ip}      = zeros(Ne,Ne,N);
        mu_ref{ip} = zeros(Ne,N);
        B_phase    = e_ref{ip};
    else
        % PHYSICAL thermodynamic response. Do not regularize it positive.
        S{ip}      = 0.5 * (R{ip}.chi + permute(R{ip}.chi,[2 1 3]));
        mu_ref{ip} = cell2mat(R{ip}.mu_e(:));

        Smu3       = pagemtimes(S{ip},reshape(mu_ref{ip},Ne,1,N));
        Smu        = reshape(Smu3,Ne,N);
        B_phase    = e_ref{ip} - Smu;
    end

    Bmix = Bmix + B_phase .* p_ip;
    Cmix = Cmix + S{ip} .* reshape(p_ip,1,1,N);
end

Ipages    = repmat(eye(Ne),1,1,N);
chi_page  = Cmix + Ipages .* reshape(1 ./ eta_vec,1,1,N);
rhs       = E_mat - Bmix;
mu3       = pagemldivide(chi_page,reshape(rhs,Ne,1,N));
mu_mat    = reshape(mu3,Ne,N);

for ip = 1:Np
    if isempty(R{ip}.mu_e) || isempty(R{ip}.chi) || isempty(R{ip}.Jac)
        e_tar{ip} = e_ref{ip};
    else
        dmu3      = reshape(mu_mat - mu_ref{ip},Ne,1,N);
        de3       = pagemtimes(S{ip},dmu3);
        e_tar{ip} = e_ref{ip} + reshape(de3,Ne,N);
    end
end

end


% =========================================================================
% Convert projected target elemental compositions to dc using de/dc only.
% =========================================================================
function dc_all = Map_TargetE_To_c_Step(c,R,e_tar)

Np = numel(c);
N  = numel(c{1}{1});
dc_all = c;

for ip = 1:Np
    for ic = 1:numel(c{ip})
        dc_all{ip}{ic} = zeros(size(c{ip}{ic}));
    end

    if isempty(R{ip}.Jac)
        continue
    end

    e_ref = cell2mat(R{ip}.e(:));
    de    = e_tar{ip} - e_ref;
    J     = R{ip}.Jac;                 % Ne x Nc x N = de/dc
    Nc    = size(J,2);

    dc3   = pagemldivide(J,reshape(de,size(de,1),1,N));
    dc    = reshape(dc3,Nc,N);

    for ic = 1:numel(c{ip})
        dc_all{ip}{ic} = dc(ic,:);
    end
end

end


% =========================================================================
% Residual of the finite-eta local-equilibrium equations.
% =========================================================================
function res = LE_Residual(R,p,E,eta,mu_mat)

Np    = numel(R);
Ne    = numel(E);
N     = numel(E{1});
E_mat = cell2mat(E(:));

if isscalar(eta)
    eta_vec = eta * ones(1,N);
else
    eta_vec = eta(:).';
end

Emix = zeros(Ne,N);
res  = zeros(1,N);

for ip = 1:Np
    p_ip = reshape(p(:,:,ip),1,N);
    e_ip = cell2mat(R{ip}.e(:));
    Emix = Emix + e_ip .* p_ip;

    if ~isempty(R{ip}.mu_e) && ~isempty(R{ip}.Jac)
        mu_ip = cell2mat(R{ip}.mu_e(:));
        r_ip  = max(abs(mu_ip - mu_mat),[],1);
        active = p_ip > 1e-14;
        r_ip(~active) = 0;
        res = max(res,r_ip);
    end
end

% finite-eta relation: mu = eta * (E - sum p e)
mu_pen = (E_mat - Emix) .* eta_vec;
res    = max(res,max(abs(mu_pen - mu_mat),[],1));

end


function c_new = AddStep(c_old,dc_all,alpha_node)
c_new = c_old;
for ip = 1:numel(c_old)
    for ic = 1:numel(c_old{ip})
        c_new{ip}{ic} = c_old{ip}{ic} + alpha_node .* dc_all{ip}{ic};
    end
end
end


function dcmax = MaxAbsStep(dc_all)
dcmax = 0;
for ip = 1:numel(dc_all)
    for ic = 1:numel(dc_all{ip})
        dcmax = max(dcmax,max(abs(dc_all{ip}{ic})));
    end
end
end


function val = MinPageEigenvalue(A)
[~,~,N] = size(A);
val = inf;
for inode = 1:N
    As = 0.5 * (A(:,:,inode) + A(:,:,inode).');
    val = min(val,min(eig(As)));
end
end
