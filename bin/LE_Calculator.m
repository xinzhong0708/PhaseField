function [c,mu_e,chi] = LE_Calculator(pars,p,c,E,eta,level)
%LE_CALCULATOR Local equilibrium with damped Newton/Picard update.
%
% Drop-in accelerated version.
%
% Main speedups compared with the original:
%   1) PhaseThermo values from the quadratic step are reused to build F_old.
%   2) Regularized Hessians are built once per phase per iteration and reused.
%   3) Line-search energy is evaluated only on nodes that are still undecided.
%   4) Hessian regularization uses pageeig when available, with loop fallback.
%
% Optional level entries:
%   level = [alpha, Miter]
%   level = [alpha, Miter, c_tol, MaxLS, lam_c]

%Prepare
c_init       =  c;
Np           =  length(c);
Ne           =  length(E);
N            =  numel(E{1});
alpha        =  level(1);
Miter        =  level(2);
c_tol        =  1e-6;

%Line-search & damping parameters
MaxLS        =  10;
amin         =  1e-7;
energy_tol   =  1e-9;
lam_c        =  1e-5;

if numel(level) >= 3 && ~isempty(level(3))
    c_tol = level(3);
end
if numel(level) >= 4 && ~isempty(level(4))
    MaxLS = level(4);
end
if numel(level) >= 5 && ~isempty(level(5))
    lam_c = level(5);
end

%Check whether any phase has internal composition degrees of freedom
has_dof      =  false(1,Np);
for ip = 1:Np
    R0          =   PhaseThermo(pars{ip}, c{ip});
    has_dof(ip) = ~(isempty(R0.mu_c) || isempty(R0.H_c) || isempty(R0.Jac));
end

%If all phases are pure/no-DOF phases, no c-iteration is needed
if ~any(has_dof)
    [~,mu_mat,chi_page] = LE_Quadratic_Step_Fast(pars,p,c,E,eta,0,0);
    mu_e = CellFromMat_Local(mu_mat,Ne);
    chi  = CellFromPage_Local(chi_page,Ne);
    return
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

for istage = 1:numel(alpha_stage)

    alpha_now = alpha_stage(istage);

    if istage == numel(alpha_stage)
        Miter_now = max(Miter,300);
    else
        Miter_now = Miter;
    end

    % Always restart from the input c for a smaller-alpha retry.
    % This avoids carrying a bad partially converged branch.
    c = c_init;

    converged_stage = false;
    failed_stage    = false;

    %Picard/Newton iteration with globalization
    for it = 1:Miter_now

        %Save old c
        c_old      = c;

        %Analytical solution of local quadratic model.
        %F_old is built from the same PhaseThermo calls, avoiding one
        %extra full LE_Objective call per iteration.
        [dc_all,~,~,F_old] = LE_Quadratic_Step_Fast(pars,p,c_old,E,eta,lam_c,1);

        %If proposed step is essentially zero, stop
        dcmax      = MaxAbsStep(dc_all);
        if dcmax < c_tol
            converged_stage = true;
            break
        end

        %Nodewise step size, used to detect unresolved bad nodes
        dcnode     = MaxAbsStepNode(dc_all,N);

        %Backtracking line search, independently for each grid point
        good_node  = false(1,N);
        alpha_try  = alpha_now*ones(1,N);
        alpha_acc  = zeros(1,N);

        %Line search
        for ils = 1:MaxLS

            bad = ~good_node & alpha_try >= amin;

            if ~any(bad)
                break
            end

            %Add c trial for all nodes. This is cheap relative to PhaseThermo.
            c_try = AddStep(c_old,dc_all,alpha_try);

            %Evaluate objective only on nodes still in the line search.
            F_try_bad = LE_Objective_Mask(pars,p,c_try,E,eta,bad);

            good_bad = isfinite(F_try_bad) & ...
                (F_try_bad <= F_old(bad) + energy_tol.*max(1,abs(F_old(bad))));

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

        %If no node accepts the step, this alpha stage failed
        if ~any(good_node)
            failed_stage = true;
            break
        end

        %Update only accepted nodes. Unaccepted nodes get alpha = 0.
        c = AddStep(c_old,dc_all,alpha_acc);

        %Check whether large-step nodes failed the line search
        bad_large = any((~good_node) & (dcnode > c_tol));

        %Check convergence
        cchg = 0;
        for ip = 1:Np
            if ~has_dof(ip)
                continue
            end
            for ic = 1:length(c{ip})
                cchg = max(cchg, max(abs(c{ip}{ic} - c_old{ip}{ic})));
            end
        end

        %Jump out if tolerance is satisfied and no large failed nodes remain
        if cchg < c_tol && ~bad_large
            converged_stage = true;
            break
        end

        %If this stage reaches the end, mark it as failed
        if it == Miter_now
            failed_stage = true;
        end
    end

    if converged_stage
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
    disp('LE_Calculator: not converged after damping retries, preserving input c.')
    c = c_init;
end

%Recalculate final mu_e and original chi at accepted c.
%The iteration above still uses the regularized Hessian.
%Only the returned chi is computed from the raw thermodynamic Hessian.
try
    [~,mu_mat,chi_page] = LE_Quadratic_Step_OriginalChi(pars,p,c,E,eta);
catch
    disp('LE_Calculator: raw chi failed, returning regularized chi.')
    [~,mu_mat,chi_page] = LE_Quadratic_Step_Fast(pars,p,c,E,eta,lam_c,0);
end

%If raw chi produced non-finite values, fall back to regularized chi.
if any(~isfinite(mu_mat(:))) || any(~isfinite(chi_page(:)))
    disp('LE_Calculator: non-finite raw response, returning regularized chi.')
    [~,mu_mat,chi_page] = LE_Quadratic_Step_Fast(pars,p,c,E,eta,lam_c,0);
end

%Diffusion potential and susceptibility
mu_e = CellFromMat_Local(mu_mat,Ne);
chi  = CellFromPage_Local(chi_page,Ne);

end

%This is the core quadratic solver for the LE
function [dc_all,mu_mat,chi_page,F] = LE_Quadratic_Step_Fast(pars,p,c,E,eta,lam_c,need_F)
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

if need_F
    Gmix = zeros(1,N);
    Emix = zeros(Ne,N);
end

for ip = 1:Np

    R{ip} = PhaseThermo(pars{ip}, c{ip});

    e_ref = cell2mat(R{ip}.e(:));
    p_ip  = reshape(p(:,:,ip), 1, N);

    if need_F
        Gmix = Gmix + p_ip .* R{ip}.g(:).';
        Emix = Emix + e_ref .* p_ip;
    end

    %Pure phase with no internal DOF
    if isempty(R{ip}.mu_c) || isempty(R{ip}.H_c) || isempty(R{ip}.Jac)
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

if need_F
    res = E_mat - Emix;
    F   = Gmix + 0.5 * eta_vec .* sum(res.^2,1);
else
    F   = [];
end

%Calculate dc = H^{-1}(J^T mu_e - mu_c)
for ip = 1:Np

    %If pure phase, no need to update c
    if isempty(R{ip}.mu_c) || isempty(R{ip}.H_c) || isempty(R{ip}.Jac)
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

%This is the final-response solver for original thermodynamic chi
function [dc_all,mu_mat,chi_page] = LE_Quadratic_Step_OriginalChi(pars,p,c,E,eta)
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
R = cell(1,Np);
for ip = 1:Np
    R{ip} = PhaseThermo(pars{ip}, c{ip});
end

%Calculate Bmix and Cmix
Bmix      = zeros(Ne, N);
Cmix      = zeros(Ne, Ne, N);

for ip = 1:Np

    e_ref = cell2mat(R{ip}.e(:));
    p_ip  = reshape(p(:,:,ip), 1, N);

    if isempty(R{ip}.mu_c) || isempty(R{ip}.H_c) || isempty(R{ip}.Jac)
        B_phase  = e_ref;
        C_phase  = zeros(Ne, Ne, N);
    else
        mu_c     = cell2mat(R{ip}.mu_c(:));
        Hc       = R{ip}.H_c;
        J        = R{ip}.Jac;
        Nc       = size(mu_c,1);
        Hreg     = SymHessian(Hc);
        Hinv_mu3 = pagemldivide(Hreg, reshape(mu_c, Nc, 1, N));
        Hinv_mu  = reshape(Hinv_mu3, Nc, N);
        JT       = permute(J, [2 1 3]);
        Hinv_JT  = pagemldivide(Hreg, JT);
        JHinvmu3 = pagemtimes(J, reshape(Hinv_mu, Nc, 1, N));
        JHinvmu  = reshape(JHinvmu3, Ne, N);
        B_phase  = e_ref - JHinvmu;
        C_phase  = pagemtimes(J, Hinv_JT);
    end

    Bmix = Bmix + B_phase .* p_ip;
    Cmix = Cmix + C_phase .* reshape(p_ip, 1, 1, N);
end

IpagesE  = repmat(eye(Ne), 1, 1, N);
chi_page = Cmix + IpagesE .* reshape(1 ./ eta_vec, 1, 1, N);
rhs      = E_mat - Bmix;
mu3      = pagemldivide(chi_page, reshape(rhs, Ne, 1, N));
mu_mat   = reshape(mu3, Ne, N);

for ip = 1:Np

    if isempty(R{ip}.mu_c) || isempty(R{ip}.H_c) || isempty(R{ip}.Jac)
        continue
    end

    mu_c    = cell2mat(R{ip}.mu_c(:));
    Hc      = R{ip}.H_c;
    J       = R{ip}.Jac;
    Nc      = size(mu_c,1);
    Hreg    = SymHessian(Hc);
    JT      = permute(J, [2 1 3]);
    JTmu3   = pagemtimes(JT, reshape(mu_mat, Ne, 1, N));
    JTmu    = reshape(JTmu3, Nc, N);
    rhs_dc  = JTmu - mu_c;
    dc3     = pagemldivide(Hreg, reshape(rhs_dc, Nc, 1, N));
    dc      = reshape(dc3, Nc, N);

    for ic = 1:length(c{ip})
        dc_all{ip}{ic} = dc(ic,:);
    end
end
end

%Evaluate the penalized energy only on selected nodes
function F = LE_Objective_Mask(pars,p,c,E,eta,mask)

Np    = length(c);
Ne    = length(E);
Nloc  = nnz(mask);

if Nloc == 0
    F = zeros(1,0);
    return
end

E_mat = cell2mat(E(:));
E_mat = E_mat(:,mask);

if isscalar(eta)
    eta_vec = eta * ones(1,Nloc);
else
    eta_vec = eta(:).';
    eta_vec = eta_vec(mask);
end

Gmix  = zeros(1,Nloc);
Emix  = zeros(Ne,Nloc);

try
    for ip = 1:Np

        c_ip = c{ip};
        c_m  = cell(size(c_ip));

        for ic = 1:numel(c_ip)
            c_m{ic} = reshape(c_ip{ic},1,[]);
            c_m{ic} = c_m{ic}(mask);
        end

        R    = PhaseThermo(pars{ip}, c_m);
        g    = R.g(:).';
        e    = cell2mat(R.e(:));
        p_ip = reshape(p(:,:,ip),1,[]);
        p_ip = p_ip(mask);

        Gmix = Gmix + p_ip .* g;
        Emix = Emix + e .* p_ip;
    end

    res = E_mat - Emix;
    F   = Gmix + 0.5 * eta_vec .* sum(res.^2,1);
catch
    F   = inf(1,Nloc);
end
end

%This function only symmetrizes Hessian and keeps the original curvature.
function Hs = SymHessian(Hc)
Hs = 0.5*(Hc + permute(Hc,[2 1 3]));
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

%This function simply update c with a step fraction alpha
function c_new = AddStep(c_old,dc_all,alpha_node)
c_new = c_old;
for ip = 1:length(c_old)
    for ic = 1:length(c_old{ip})
        c_new{ip}{ic} = c_old{ip}{ic} + alpha_node.*dc_all{ip}{ic};
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
