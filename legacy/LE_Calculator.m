function [c,mu_e,chi] = LE_Calculator(pars,p,c,E,eta,level)
%This function calculates the LE for any number of phase with damping
%pars:   input parameter {Np}
%p:      phase fraction 1*n*Np
%c:      endmember concentration {Np}{Nc}
%E:      bulk composition {Ne}
%eta:    penalty
%level: [alpha=damping factor, Miter=maximal iteration]

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

%Check whether any phase has internal composition degrees of freedom
has_dof      =  false(1,Np);
for ip = 1:Np
    R0          =   PhaseThermo(pars{ip}, c{ip});
    has_dof(ip) = ~(isempty(R0.mu_c) || isempty(R0.H_c) || isempty(R0.Jac));
end

%If all phases are pure/no-DOF phases, no c-iteration is needed
if ~any(has_dof)
    [~,mu_mat,chi_page] = LE_Quadratic_Step(pars,p,c,E,eta,0);
    for ie = 1:Ne
        mu_e{ie} = mu_mat(ie,:);
    end
    for i = 1:Ne
        for j = 1:Ne
            chi{i,j} = reshape(chi_page(i,j,:),1,[]);
        end
    end
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

        %True penalized energy before update
        F_old      = LE_Objective(pars,p,c_old,E,eta);

        %Analytical solution of local quadratic model
        dc_all     = LE_Quadratic_Step(pars,p,c_old,E,eta,lam_c);

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

            %Add c trial
            c_try = AddStep(c_old,dc_all,alpha_try);

            %Calculate true penalized energy
            F_try = LE_Objective(pars,p,c_try,E,eta);

            %Accept non-increasing energy within numerical tolerance
            good = isfinite(F_try) & (F_try <= F_old + energy_tol.*max(1,abs(F_old)));

            %Accept newly good nodes
            good_new            = good & ~good_node;
            alpha_acc(good_new) = alpha_try(good_new);
            good_node(good_new) = true;

            %Reduce bad nodes
            bad = ~good_node;
            alpha_try(bad) = 0.2 * alpha_try(bad);

            %Stop if all accepted or remaining steps are too small
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
%The iteration above still uses the legacy regularized Hessian.
%Only the returned chi is computed from the raw thermodynamic Hessian.
try
    [~,mu_mat,chi_page] = LE_Quadratic_Step_OriginalChi(pars,p,c,E,eta);
catch
    disp('LE_Calculator: raw chi failed, returning regularized chi.')
    [~,mu_mat,chi_page] = LE_Quadratic_Step(pars,p,c,E,eta,lam_c);
end

%If raw chi produced non-finite values, fall back to regularized chi.
if any(~isfinite(mu_mat(:))) || any(~isfinite(chi_page(:)))
    disp('LE_Calculator: non-finite raw response, returning regularized chi.')
    [~,mu_mat,chi_page] = LE_Quadratic_Step(pars,p,c,E,eta,lam_c);
end

%Diffusion potential
for ie = 1:Ne
    mu_e{ie} = mu_mat(ie,:);
end

%Susceptibility
for i = 1:Ne
    for j = 1:Ne
        chi{i,j} = reshape(chi_page(i,j,:),1,[]);
    end
end

end

%This is the core quadratic solver for the LE
function [dc_all,mu_mat,chi_page] = LE_Quadratic_Step(pars,p,c,E,eta,lam_c)
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
dc_all = c;
for ip = 1:Np
    for ic = 1:length(c{ip})
        dc_all{ip}{ic} = zeros(size(c{ip}{ic}));
    end
end
%Prepare thermodynamic for each phase
for ip = 1:Np
    R{ip} = PhaseThermo(pars{ip}, c{ip});
end
%Calculate Bmix and Cmix
Bmix      = zeros(Ne, N);
Cmix      = zeros(Ne, Ne, N);
%Loop through phases
for ip = 1:Np
    %Current e
    e_ref = cell2mat(R{ip}.e(:));
    p_ip  = reshape(p(:,:,ip), 1, N);
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
        %Regularized positive definite Hessian
        Hreg     = RegularizeHessian(Hc,lam_c);
        %H^{-1} * mu_c
        Hinv_mu3 = pagemldivide(Hreg, reshape(mu_c, Nc, 1, N));
        Hinv_mu  = reshape(Hinv_mu3, Nc, N);
        %H^{-1} * J^T
        JT       = permute(J, [2 1 3]);
        Hinv_JT  = pagemldivide(Hreg, JT);
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
    if isempty(R{ip}.mu_c) || isempty(R{ip}.H_c) || isempty(R{ip}.Jac)
        continue
    end
    mu_c    = cell2mat(R{ip}.mu_c(:));
    Hc      = R{ip}.H_c;
    J       = R{ip}.Jac;
    Nc      = size(mu_c,1);
    Hreg    = RegularizeHessian(Hc,lam_c);
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
dc_all = c;
for ip = 1:Np
    for ic = 1:length(c{ip})
        dc_all{ip}{ic} = zeros(size(c{ip}{ic}));
    end
end
%Prepare thermodynamic for each phase
for ip = 1:Np
    R{ip} = PhaseThermo(pars{ip}, c{ip});
end
%Calculate Bmix and Cmix
Bmix      = zeros(Ne, N);
Cmix      = zeros(Ne, Ne, N);
%Loop through phases
for ip = 1:Np
    %Current e
    e_ref = cell2mat(R{ip}.e(:));
    p_ip  = reshape(p(:,:,ip), 1, N);
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
        %Raw symmetric Hessian for original thermodynamic response
        Hreg     = SymHessian(Hc);
        %H^{-1} * mu_c
        Hinv_mu3 = pagemldivide(Hreg, reshape(mu_c, Nc, 1, N));
        Hinv_mu  = reshape(Hinv_mu3, Nc, N);
        %H^{-1} * J^T
        JT       = permute(J, [2 1 3]);
        Hinv_JT  = pagemldivide(Hreg, JT);
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


%This function evaluate the penalized energy 
function F = LE_Objective(pars,p,c,E,eta)
%Prepare
Np    = length(c);
Ne    = length(E);
N     = numel(E{1});
E_mat = cell2mat(E(:));
if isscalar(eta)
    eta_vec = eta * ones(1,N);
else
    eta_vec = eta(:).';
end
Gmix  = zeros( 1,N);
Emix  = zeros(Ne,N);
try
    for ip = 1:Np
        R    = PhaseThermo(pars{ip}, c{ip});
        g    = R.g(:).';
        e    = cell2mat(R.e(:));
        p_ip = reshape(p(:,:,ip),1,N);
        Gmix = Gmix + p_ip .* g;
        Emix = Emix + e .* p_ip;
    end
    res = E_mat - Emix;
    F   = Gmix + 0.5 * eta_vec .* sum(res.^2,1);
catch
    F   = inf(1,N);
end
end



%This function only symmetrizes Hessian and keeps the original curvature.
%It does not shift negative eigenvalues and therefore does not force chi
%to be positive. Use only for final returned chi, not for the iteration
%direction.
function Hs = SymHessian(Hc)
Hs = 0.5*(Hc + permute(Hc,[2 1 3]));
end

%This function regularize Hessian by adding a ridge and make it symmetric
function Hreg = RegularizeHessian(Hc,lam_c)
[Nc,~,N]        = size(Hc);
Hreg            = zeros(Nc,Nc,N);
for i = 1:N
    H           = 0.5 * (Hc(:,:,i) + Hc(:,:,i).');
    scale       = max(1,norm(H,'fro')/max(1,Nc));
    evmin       = min(eig(H));
    shift       = max(0,1e-10*scale - evmin);
    Hreg(:,:,i) = H + (shift + lam_c*scale) * eye(Nc);
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