function [c,mu_e,chi,DIAG] = LE_Calculator_SplitE(pars,p,c,c_ref,E,eta,level,NUM)
%LE_CALCULATOR_SPLITE E-conserving convex-split local equilibrium closure.
%
% Solves, at fixed E and p,
%
%   mu_c(c_ip) + S_ip*(c_ip-c_ref_ip) - J_ip' * mu_e = 0
%   E - sum_ip p_ip*e_ip(c_ip) - mu_e/eta = 0
%
% where S_ip is fixed from the accepted reference composition c_ref.  The
% matrix S is zero in locally convex regions and shifts negative curvature
% in a spinodal region into an explicit (old-composition) contribution.
% Thus chi returned here is the smooth implicit map tangent used by the
% AC-CH solver, whereas DIAG.lambda_min_phys reports raw physical curvature.
%
% Inputs follow LE_Calculator.m, with two additions:
%   c_ref : accepted old-time composition, same cell structure as c
%   NUM   : SplitE controls (all optional)
%
% Suggested controls:
%   NUM.split_floor_rel      = 1e-8;  % positive map eigenvalue floor
%   NUM.split_safety_fac     = 1.05;  % extra shift when raw H is negative
%   NUM.split_map_floor_rel  = 1e-10; % Newton-only conditioning floor
%   NUM.split_newton_ridge   = 1e-8;  % Newton-only ridge
%   NUM.split_c_tol          = 1e-7;
%   NUM.split_stat_tol       = 1e-7;
%
% The returned mu_e is a convex-split chemical potential during transient
% evolution.  When c == c_ref, the S contribution vanishes and it reduces
% to the physical local-equilibrium chemical potential.

if nargin < 8 || isempty(NUM)
    NUM = struct();
end

% Controls
alpha        = level(1);
Miter        = level(2);
c_tol        = GetControl(NUM,'split_c_tol',1e-7);
stat_tol     = GetControl(NUM,'split_stat_tol',1e-7);
MaxLS        = GetControl(NUM,'split_MaxLS',12);
amin         = GetControl(NUM,'split_amin',1e-9);
energy_tol   = GetControl(NUM,'split_energy_tol',1e-10);
lam_c        = GetControl(NUM,'split_newton_ridge',1e-8);

Np           = length(c);
Ne           = length(E);
N            = numel(E{1});

% Fixed convex-split shift from old accepted state. Never recompute this
% shift from the nonlinear iterate inside one attempted physical timestep.
[S,lambda_ref,shift_ref] = Build_Split_Shift(pars,c_ref,NUM);

% Detect internal phase composition degrees of freedom.
has_dof = false(1,Np);
for ip = 1:Np
    R0          = PhaseThermo(pars{ip},c{ip});
    has_dof(ip) = ~(isempty(R0.mu_c) || isempty(R0.H_c) || isempty(R0.Jac));
end

converged    = false;
failed_nodes = false(1,N);
it_done      = 0;

% Pure/no-DOF subsets need only the finite-eta map.
if ~any(has_dof)
    [~,mu_mat,chi_page,AUX] = Split_Quadratic_Step(pars,p,c,c_ref,E,eta,S,0,NUM);
    [mu_e,chi] = Mat_To_Output(mu_mat,chi_page,Ne);
    [stat_node,mass_node] = Split_Residual(pars,p,c,c_ref,E,mu_e,eta,S);
    converged = true;
else
    % Newton/Picard iteration with line search on the convex-split local
    % energy. Unlike ordinary LE minimization, the branch regularization is
    % continuous and is fixed by c_ref through the attempted timestep.
    for it = 1:Miter
        it_done = it;
        c_old   = c;
        F_old   = LE_Objective_Split(pars,p,c_old,c_ref,E,eta,S);

        [dc_all,~,~,~] = Split_Quadratic_Step(pars,p,c_old,c_ref,E,eta,S,lam_c,NUM);
        dc_all          = Limit_Nodewise_Step(dc_all,GetControl(NUM,'split_dc_step_max',inf));
        dcmax           = MaxAbsStep(dc_all);

        if dcmax < c_tol
            converged = true;
            break
        end

        good_node = false(1,N);
        alpha_try = alpha*ones(1,N);
        alpha_acc = zeros(1,N);

        for ils = 1:MaxLS
            c_try = AddStep(c_old,dc_all,alpha_try);
            F_try = LE_Objective_Split(pars,p,c_try,c_ref,E,eta,S);
            good  = isfinite(F_try) & ...
                    (F_try <= F_old + energy_tol.*max(1,abs(F_old)));

            good_new = good & ~good_node;
            alpha_acc(good_new) = alpha_try(good_new);
            good_node(good_new) = true;

            bad = ~good_node;
            alpha_try(bad) = 0.2*alpha_try(bad);
            if all(good_node | alpha_try < amin)
                break
            end
        end

        failed_nodes = failed_nodes | ~good_node;
        if any(good_node)
            c = AddStep(c_old,dc_all,alpha_acc);
        else
            % Keep the last finite state. The timestep controller should
            % reject/shrink if this yields a large residual.
            c = c_old;
            break
        end

        cchg = MaxAbsDifference(c,c_old,has_dof);
        if cchg < c_tol
            converged = true;
            break
        end
    end

    [~,mu_mat,chi_page,AUX] = Split_Quadratic_Step(pars,p,c,c_ref,E,eta,S,0,NUM);
    [mu_e,chi] = Mat_To_Output(mu_mat,chi_page,Ne);
    [stat_node,mass_node] = Split_Residual(pars,p,c,c_ref,E,mu_e,eta,S);
    if max(stat_node) <= stat_tol
        converged = true;
    end
end

% Diagnostics. lambda_min_phys is raw H_c at final c, not chi_map.
lambda_now = Calc_Raw_Min_Eigenvalue(pars,c);
DIAG.converged           = converged;
DIAG.iter                 = it_done;
DIAG.failed_nodes         = failed_nodes;
DIAG.max_stationarity     = max(stat_node);
DIAG.max_mass_residual    = max(mass_node);
DIAG.stat_node            = stat_node;
DIAG.mass_node            = mass_node;
DIAG.lambda_min_ref       = lambda_ref;
DIAG.lambda_min_phys      = lambda_now;
DIAG.split_shift          = shift_ref;
DIAG.extra_map_shift      = AUX.extra_shift;
DIAG.spinodal_ref         = lambda_ref < 0;
DIAG.spinodal_phys        = lambda_now < 0;
DIAG.max_split_shift      = max(shift_ref(:));
DIAG.max_extra_map_shift  = max(AUX.extra_shift(:));

end


function [dc_all,mu_mat,chi_page,AUX] = Split_Quadratic_Step(pars,p,c,c_ref,E,eta,S,lam_c,NUM)
%SPLIT_QUADRATIC_STEP Local tangent step for the E-primary SplitE map.

Np    = length(c);
Ne    = length(E);
N     = numel(E{1});
E_mat = cell2mat(E(:));

if isscalar(eta)
    eta_vec = eta*ones(1,N);
else
    eta_vec = eta(:).';
end

% Initialize composition step and calculate phase thermodynamics.
dc_all = c;
R      = cell(1,Np);
for ip = 1:Np
    for ic = 1:length(c{ip})
        dc_all{ip}{ic} = zeros(size(c{ip}{ic}));
    end
    R{ip} = PhaseThermo(pars{ip},c{ip});
end

Bmix        = zeros(Ne,N);
Cmix        = zeros(Ne,Ne,N);
AUX.extra_shift      = zeros(Np,N);
AUX.lambda_min_phys  = nan(Np,N);
Hmap_all    = cell(1,Np);
mu_eff_all  = cell(1,Np);

for ip = 1:Np
    e_ref = cell2mat(R{ip}.e(:));
    p_ip  = reshape(p(:,:,ip),1,N);

    if isempty(R{ip}.mu_c) || isempty(R{ip}.H_c) || isempty(R{ip}.Jac)
        B_phase = e_ref;
        C_phase = zeros(Ne,Ne,N);
        Hmap_all{ip}   = [];
        mu_eff_all{ip} = [];
    else
        mu_c    = cell2mat(R{ip}.mu_c(:));
        J       = R{ip}.Jac;
        Nc      = size(mu_c,1);
        dc_ref  = CellCompositionMatrix(c{ip}) - CellCompositionMatrix(c_ref{ip});
        Sdc3    = pagemtimes(S{ip},reshape(dc_ref,Nc,1,N));
        mu_eff  = mu_c + reshape(Sdc3,Nc,N);

        [Hmap,extra_shift,lambda_min] = Regularize_Map_Hessian( ...
            R{ip}.H_c,S{ip},lam_c,NUM);
        AUX.extra_shift(ip,:)     = extra_shift;
        AUX.lambda_min_phys(ip,:) = lambda_min;
        Hmap_all{ip}              = Hmap;
        mu_eff_all{ip}            = mu_eff;

        Hinv_mu3 = pagemldivide(Hmap,reshape(mu_eff,Nc,1,N));
        Hinv_mu  = reshape(Hinv_mu3,Nc,N);
        JT       = permute(J,[2 1 3]);
        Hinv_JT  = pagemldivide(Hmap,JT);
        JHinvmu3 = pagemtimes(J,reshape(Hinv_mu,Nc,1,N));
        JHinvmu  = reshape(JHinvmu3,Ne,N);
        B_phase  = e_ref - JHinvmu;
        C_phase  = pagemtimes(J,Hinv_JT);
    end

    Bmix = Bmix + B_phase .* p_ip;
    Cmix = Cmix + C_phase .* reshape(p_ip,1,1,N);
end

% Smooth E -> mu map tangent used by the CH and corrector solvers.
IpagesE   = repmat(eye(Ne),1,1,N);
chi_page  = Cmix + IpagesE .* reshape(1./eta_vec,1,1,N);
rhs       = E_mat - Bmix;
mu3       = pagemldivide(chi_page,reshape(rhs,Ne,1,N));
mu_mat    = reshape(mu3,Ne,N);

% Composition correction at the calculated split chemical potential.
for ip = 1:Np
    if isempty(Hmap_all{ip})
        continue
    end
    J       = R{ip}.Jac;
    Nc      = size(mu_eff_all{ip},1);
    JT      = permute(J,[2 1 3]);
    JTmu3   = pagemtimes(JT,reshape(mu_mat,Ne,1,N));
    JTmu    = reshape(JTmu3,Nc,N);
    rhs_dc  = JTmu - mu_eff_all{ip};
    dc3     = pagemldivide(Hmap_all{ip},reshape(rhs_dc,Nc,1,N));
    dc      = reshape(dc3,Nc,N);
    for ic = 1:length(c{ip})
        dc_all{ip}{ic} = dc(ic,:);
    end
end

end


function F = LE_Objective_Split(pars,p,c,c_ref,E,eta,S)
%LE_OBJECTIVE_SPLIT Convex-split, finite-eta local objective.
% Its stationarity gives mu_c + S(c-c_ref) - J'*mu_e = 0.

Np    = length(c);
Ne    = length(E);
N     = numel(E{1});
E_mat = cell2mat(E(:));

if isscalar(eta)
    eta_vec = eta*ones(1,N);
else
    eta_vec = eta(:).';
end

Gmix = zeros(1,N);
Emix = zeros(Ne,N);

try
    for ip = 1:Np
        R    = PhaseThermo(pars{ip},c{ip});
        g    = R.g(:).';
        e    = cell2mat(R.e(:));
        p_ip = reshape(p(:,:,ip),1,N);
        q    = zeros(1,N);
        if ~isempty(S{ip})
            dc_ref = CellCompositionMatrix(c{ip}) - CellCompositionMatrix(c_ref{ip});
            Nc     = size(dc_ref,1);
            Sdc3   = pagemtimes(S{ip},reshape(dc_ref,Nc,1,N));
            Sdc    = reshape(Sdc3,Nc,N);
            q      = 0.5*sum(dc_ref.*Sdc,1);
        end
        Gmix = Gmix + p_ip .* (g + q);
        Emix = Emix + e .* p_ip;
    end
    res = E_mat - Emix;
    F   = Gmix + 0.5*eta_vec.*sum(res.^2,1);
catch
    F = inf(1,N);
end

end


function [stat_node,mass_node] = Split_Residual(pars,p,c,c_ref,E,mu_e,eta,S)
%SPLIT_RESIDUAL Evaluate modified stationarity and finite-eta mass closure.

Np     = length(c);
Ne     = length(E);
N      = numel(E{1});
mu_mat = cell2mat(mu_e(:));
E_mat  = cell2mat(E(:));

if isscalar(eta)
    eta_vec = eta*ones(1,N);
else
    eta_vec = eta(:).';
end

stat_node = zeros(1,N);
Emix      = zeros(Ne,N);

for ip = 1:Np
    R    = PhaseThermo(pars{ip},c{ip});
    e    = cell2mat(R.e(:));
    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + e.*p_ip;

    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end
    mu_c   = cell2mat(R.mu_c(:));
    dc_ref = CellCompositionMatrix(c{ip}) - CellCompositionMatrix(c_ref{ip});
    Nc     = size(mu_c,1);
    Sdc3   = pagemtimes(S{ip},reshape(dc_ref,Nc,1,N));
    mu_eff = mu_c + reshape(Sdc3,Nc,N);
    JT     = permute(R.Jac,[2 1 3]);
    JTmu3  = pagemtimes(JT,reshape(mu_mat,Ne,1,N));
    JTmu   = reshape(JTmu3,Nc,N);
    stat_node = max(stat_node,max(abs(mu_eff-JTmu),[],1));
end

mass_res  = E_mat - Emix - mu_mat./eta_vec;
mass_node = max(abs(mass_res),[],1);

end


function [S,lambda_min,split_shift] = Build_Split_Shift(pars,c_ref,NUM)
%BUILD_SPLIT_SHIFT Calculate the fixed explicit-concavity shift from c_ref.

Np          = length(c_ref);
N           = numel(c_ref{1}{1});
floor_rel   = GetControl(NUM,'split_floor_rel',1e-8);
safety_fac  = GetControl(NUM,'split_safety_fac',1.05);
S           = cell(1,Np);
lambda_min  = nan(Np,N);
split_shift = zeros(Np,N);

for ip = 1:Np
    R = PhaseThermo(pars{ip},c_ref{ip});
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        S{ip} = [];
        continue
    end
    [Nc,~,~] = size(R.H_c);
    S{ip}    = zeros(Nc,Nc,N);
    for inode = 1:N
        H       = 0.5*(R.H_c(:,:,inode)+R.H_c(:,:,inode).');
        scale   = max(1,norm(H,'fro')/max(1,Nc));
        lam     = min(eig(H));
        target  = floor_rel*scale;
        shift   = max(0,target-lam);
        if shift > 0
            shift = safety_fac*shift;
        end
        lambda_min(ip,inode)  = lam;
        split_shift(ip,inode) = shift;
        S{ip}(:,:,inode)      = shift*eye(Nc);
    end
end

end


function [Hmap,extra_shift,lambda_min] = Regularize_Map_Hessian(Hc,S,lam_c,NUM)
%REGULARIZE_MAP_HESSIAN Condition the already split implicit tangent.
% Any extra_shift here is numerical Newton conditioning only and is reported
% separately; a persistent large value means split_safety_fac is too small.

[Nc,~,N]     = size(Hc);
floor_rel    = GetControl(NUM,'split_map_floor_rel',1e-10);
Hmap         = zeros(Nc,Nc,N);
extra_shift  = zeros(1,N);
lambda_min   = zeros(1,N);

for inode = 1:N
    Hphys     = 0.5*(Hc(:,:,inode)+Hc(:,:,inode).');
    lambda_min(inode) = min(eig(Hphys));
    H         = Hphys + S(:,:,inode);
    scale     = max(1,norm(H,'fro')/max(1,Nc));
    target    = floor_rel*scale;
    lam_map   = min(eig(H));
    extra     = max(0,target-lam_map) + lam_c*scale;
    extra_shift(inode) = extra;
    Hmap(:,:,inode)    = H + extra*eye(Nc);
end

end


function lambda_min = Calc_Raw_Min_Eigenvalue(pars,c)
%CALC_RAW_MIN_EIGENVALUE Raw H_c curvature of the final composition.

Np         = length(c);
N          = numel(c{1}{1});
lambda_min = nan(Np,N);

for ip = 1:Np
    R = PhaseThermo(pars{ip},c{ip});
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end
    for inode = 1:N
        H = 0.5*(R.H_c(:,:,inode)+R.H_c(:,:,inode).');
        lambda_min(ip,inode) = min(eig(H));
    end
end

end


function M = CellCompositionMatrix(c_ip)
%CELLCOMPOSITIONMATRIX independent endmember variables x spatial nodes.
if isempty(c_ip)
    M = [];
else
    M = cell2mat(c_ip(:));
end
end


function [mu_e,chi] = Mat_To_Output(mu_mat,chi_page,Ne)
mu_e = cell(1,Ne);
chi  = cell(Ne,Ne);
for ie = 1:Ne
    mu_e{ie} = mu_mat(ie,:);
end
for i = 1:Ne
    for j = 1:Ne
        chi{i,j} = reshape(chi_page(i,j,:),1,[]);
    end
end
end


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


function dmax = MaxAbsDifference(c,c_old,has_dof)
dmax = 0;
for ip = 1:length(c)
    if ~has_dof(ip)
        continue
    end
    for ic = 1:length(c{ip})
        dmax = max(dmax,max(abs(c{ip}{ic}-c_old{ip}{ic})));
    end
end
end


function dc_out = Limit_Nodewise_Step(dc_in,dc_limit)
%LIMIT_NODEWISE_STEP Optional trust radius preventing a one-call branch jump.
dc_out = dc_in;
if isempty(dc_limit) || ~isfinite(dc_limit) || dc_limit <= 0
    return
end
N = numel(dc_in{1}{1});
node_max = zeros(1,N);
for ip = 1:length(dc_in)
    for ic = 1:length(dc_in{ip})
        node_max = max(node_max,abs(dc_in{ip}{ic}));
    end
end
fac = min(1,dc_limit./max(node_max,eps));
for ip = 1:length(dc_in)
    for ic = 1:length(dc_in{ip})
        dc_out{ip}{ic} = dc_in{ip}{ic}.*fac;
    end
end
end


function value = GetControl(NUM,name,default_value)
if isstruct(NUM) && isfield(NUM,name) && ~isempty(NUM.(name))
    value = NUM.(name);
else
    value = default_value;
end
end
