function [c,mu_e,chi,DIAG] = LE_Calculator_PathE(pars,p,c,mu_e,E,eta,level,NUM)
%LE_CALCULATOR_PATHE Vectorized branch-continuation LE closure at fixed E,p.
%
% This is the vectorized replacement for the node-loop PathE version.
% There is NO loop over spatial grid nodes.  PhaseThermo is evaluated for
% all N spatial points at once and the coupled Newton systems are solved by
% page-wise linear algebra:
%
%   mu_c^ip(c_ip) - J_ip' * mu_e = 0
%   sum_ip p_ip e_ip(c_ip) + mu_e/eta - E = 0
%
% Unknown on each page/node:
%   z = [c_1 ; c_2 ; ... ; c_Np ; mu_e].
%
% The solver follows the incoming local branch because it starts from the
% incoming c and mu_e and minimizes the exact equation residual, rather than
% minimizing the local free energy.  Returned mu_e is reconstructed from the
% exact finite-eta mass relation at the returned c:
%
%   mu_e = eta .* (E - sum_ip p_ip e_ip(c_ip)).
%
% Loops retained in this function are only over thermodynamic phases,
% compositional coordinates, Newton iterations and backtracking iterations;
% those dimensions are small.  No loop scales with nx*ny.
%
% Required shapes, matching LE_Run_PathE:
%   p(1,N,Np), c{ip}{ic}(1,N), mu_e{ie}(1,N), E{ie}(1,N).
%
% Useful NUM controls:
%   path_stat_tol          stationarity tolerance, default 1e-8
%   path_mass_tol          mass residual tolerance, default 1e-10
%   path_MaxLS             line-search iterations, default 12
%   path_c_floor           composition interior floor, default 1e-12
%   path_dc_max            maximum |dc| in one LE call, default 5e-2
%   path_dmu_max           maximum |dmu| in one LE call, default inf
%   path_newton_ridge_rel  tiny page solve ridge, default 1e-12
%   path_H_ridge_rel       tiny chi-evaluation ridge, default 1e-12
%   path_use_chi_solver    0 returns physical chi; 1 returns PSD-conditioned
%
% For a physical spinodal calculation, retain path_use_chi_solver = 0.

if nargin < 8 || isempty(NUM); NUM = struct(); end
if nargin < 7 || isempty(level); level = [0.7,100]; end

alpha0      = level(1);
MaxIter     = level(2);
stat_tol    = GetControl(NUM,'path_stat_tol',1e-8);
mass_tol    = GetControl(NUM,'path_mass_tol',1e-10);
MaxLS       = GetControl(NUM,'path_MaxLS',12);
amin        = GetControl(NUM,'path_alpha_min',1e-10);
beta_ls     = GetControl(NUM,'path_ls_beta',0.25);
armijo      = GetControl(NUM,'path_armijo',1e-4);
c_floor     = GetControl(NUM,'path_c_floor',1e-12);
dc_max      = GetControl(NUM,'path_dc_max',5e-2);
dmu_max     = GetControl(NUM,'path_dmu_max',inf);
ridge_rel   = GetControl(NUM,'path_newton_ridge_rel',1e-12);
use_chi_sol = GetControl(NUM,'path_use_chi_solver',0);

Np = numel(c);
Ne = numel(E);
N  = numel(E{1});
E_mat  = cell2mat(E(:));
mu_mat = cell2mat(mu_e(:));
if isscalar(eta)
    eta_vec = eta*ones(1,N);
else
    eta_vec = reshape(eta,1,N);
end

[z,layout] = Pack_Unknown_Batch(c,mu_mat);
Nz = layout.Nz;

% Scaling affects only residual norm/step conditioning; it does not change
% the exact root. Scales are kept fixed over one LE call.
mu_scale = max(1,max(abs(mu_mat),[],1));
E_scale  = max(1,max(abs(E_mat),[],1));
if isfield(NUM,'path_mu_scale') && ~isempty(NUM.path_mu_scale)
    mu_scale(:) = NUM.path_mu_scale;
end
if isfield(NUM,'path_E_scale') && ~isempty(NUM.path_E_scale)
    E_scale(:) = NUM.path_E_scale;
end

conv_node   = false(1,N);
failed_node = false(1,N);
iter_node   = zeros(1,N);
ls_node     = zeros(1,N);

for it = 1:MaxIter
    [r,J,AUX] = Residual_Jacobian_Batch(pars,p,z,E_mat,eta_vec,layout);
    [rs,Js]   = Scale_Residual_Jacobian(r,J,layout,mu_scale,E_scale);
    merit0    = max(abs(rs),[],1);

    new_conv  = AUX.stat_node <= stat_tol & AUX.mass_node <= mass_tol;
    conv_node = conv_node | new_conv;
    work      = ~conv_node & ~failed_node;
    iter_node(work) = it;
    if ~any(work)
        break
    end

    dz = Page_Newton_Direction(Js,rs,ridge_rel);
    dz(:,~work) = 0;
    dz = Limit_Step_Batch(dz,layout,dc_max,dmu_max);

    a_bound = Fraction_To_Boundary_Batch(z,dz,layout,c_floor);
    a_try   = zeros(1,N);
    a_try(work) = min(alpha0,a_bound(work));
    accepted = false(1,N);
    a_acc    = zeros(1,N);

    for ils = 1:MaxLS
        pending = work & ~accepted & a_try >= amin;
        if ~any(pending)
            break
        end
        z_try = z + dz .* a_try;
        [r_try,~,AUX_try] = Residual_Jacobian_Batch(pars,p,z_try,E_mat,eta_vec,layout);
        rs_try    = Scale_Residual_Only(r_try,layout,mu_scale,E_scale);
        merit_try = max(abs(rs_try),[],1);
        good = pending & isfinite(merit_try) & ...
            (merit_try <= merit0.*(1-armijo*a_try) | ...
             (AUX_try.stat_node <= stat_tol & AUX_try.mass_node <= mass_tol));
        a_acc(good)    = a_try(good);
        accepted(good) = true;
        a_try(pending & ~good) = beta_ls*a_try(pending & ~good);
        ls_node(pending) = ils;
    end

    if any(accepted)
        z(:,accepted) = z(:,accepted) + dz(:,accepted).*a_acc(accepted);
    end
    failed_node(work & ~accepted) = true;
end

% Recover c, then enforce the finite-eta mass relation exactly.  This is
% required even if a Newton/page solve has failed at a subset of nodes.
[c,~] = Unpack_Unknown_Batch(z,layout,c);
[Emix,chi_phys_page,lambda_min,lambda_exact] = ...
    Exact_Emix_Chi_Batch(pars,p,c,eta_vec,Ne,NUM);
mu_mat = (E_mat - Emix) .* eta_vec;
z(layout.mu_ids,:) = mu_mat;

% Final residual is always measured with exact mass-compatible mu_e.
[r_final,~,AUX] = Residual_Jacobian_Batch(pars,p,z,E_mat,eta_vec,layout);
rs_final = Scale_Residual_Only(r_final,layout,mu_scale,E_scale);
merit_node = max(abs(rs_final),[],1);
conv_node = AUX.stat_node <= stat_tol & AUX.mass_node <= mass_tol;
failed_node = failed_node | ~conv_node;

[chi_solver_page,chi_shift_node] = Condition_Chi_Pages(chi_phys_page,NUM);
if use_chi_sol == 1
    chi_page = chi_solver_page;
    chi_mode = 'conditioned_solver';
else
    chi_page = chi_phys_page;
    chi_mode = 'physical';
end

for ie = 1:Ne
    mu_e{ie} = mu_mat(ie,:);
end
chi = cell(Ne,Ne);
for ie = 1:Ne
    for je = 1:Ne
        chi{ie,je} = reshape(chi_page(ie,je,:),1,N);
    end
end

DIAG.converged          = all(conv_node);
DIAG.converged_node     = conv_node;
DIAG.failed_node        = failed_node;
DIAG.iter_node          = iter_node;
DIAG.max_iter           = max(iter_node);
DIAG.ls_node            = ls_node;
DIAG.stat_node          = AUX.stat_node;
DIAG.mass_node          = AUX.mass_node;
DIAG.merit_node         = merit_node;
DIAG.max_stationarity   = max(AUX.stat_node);
DIAG.max_mass_residual  = max(AUX.mass_node);
DIAG.max_merit          = max(merit_node);
DIAG.lambda_min_phys    = lambda_min;
DIAG.lambda_exact       = lambda_exact;
DIAG.min_lambda_phys    = min(lambda_min(:),[],'omitnan');
DIAG.spinodal_node      = any(lambda_min < 0,1);
DIAG.spinodal_fraction  = mean(DIAG.spinodal_node);
DIAG.chi_phys_page      = chi_phys_page;
DIAG.chi_solver_page    = chi_solver_page;
DIAG.chi_shift_node     = chi_shift_node;
DIAG.max_chi_shift      = max(chi_shift_node);
DIAG.returned_chi_mode  = chi_mode;
DIAG.vectorized         = true;
DIAG.N_pages            = N;

end


function [r,J,AUX] = Residual_Jacobian_Batch(pars,p,z,E_mat,eta_vec,layout)
% Exact root residual/Jacobian on all spatial pages simultaneously.
Np = layout.Np;
Ne = layout.Ne;
N  = layout.N;
Nz = layout.Nz;
mu = z(layout.mu_ids,:);
r  = zeros(Nz,N);
J  = zeros(Nz,Nz,N);
Emix = zeros(Ne,N);

for ip = 1:Np
    ci = Cell_From_Z_Phase(z,layout,ip);
    R  = PhaseThermo(pars{ip},ci);
    ei = cell2mat(R.e(:));
    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + ei.*p_ip;

    ids = layout.c_ids{ip};
    if isempty(ids) || isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end
    mu_c = cell2mat(R.mu_c(:));
    Hc   = 0.5*(R.H_c + permute(R.H_c,[2 1 3]));
    Jac  = R.Jac;
    JT   = permute(Jac,[2 1 3]);
    JTmu = reshape(pagemtimes(JT,reshape(mu,Ne,1,N)),numel(ids),N);
    r(ids,:) = mu_c - JTmu;
    J(ids,ids,:)          = Hc;
    J(ids,layout.mu_ids,:)= -JT;
    J(layout.E_ids,ids,:) = J(layout.E_ids,ids,:) + ...
        Jac.*reshape(p_ip,1,1,N);
end

r(layout.E_ids,:) = Emix + mu./eta_vec - E_mat;
J(layout.E_ids,layout.mu_ids,:) = repmat(eye(Ne),1,1,N) ./ reshape(eta_vec,1,1,N);

if layout.Nc_tot > 0
    AUX.stat_node = max(abs(r(1:layout.Nc_tot,:)),[],1);
else
    AUX.stat_node = zeros(1,N);
end
AUX.mass_node = max(abs(r(layout.E_ids,:)),[],1);
end


function [Emix,chi_page,lambda_min,lambda_exact] = Exact_Emix_Chi_Batch(pars,p,c,eta_vec,Ne,NUM)
% Compute mass mixture and physical susceptibility without looping over nodes.
Np = numel(c);
N  = numel(eta_vec);
Emix = zeros(Ne,N);
chi_page = repmat(eye(Ne),1,1,N) ./ reshape(eta_vec,1,1,N);
lambda_min   = nan(Np,N);
lambda_exact = false(Np,1);
Hridge_rel   = GetControl(NUM,'path_H_ridge_rel',1e-12);

for ip = 1:Np
    R = PhaseThermo(pars{ip},c{ip});
    ei = cell2mat(R.e(:));
    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + ei.*p_ip;
    if isempty(R.H_c) || isempty(R.Jac)
        continue
    end
    H   = 0.5*(R.H_c + permute(R.H_c,[2 1 3]));
    Jac = R.Jac;
    Nc  = size(H,1);
    [lambda_min(ip,:),lambda_exact(ip)] = Min_Eig_Sym_Pages(H);
    hscale = max(1,reshape(sqrt(sum(H.^2,[1 2])),1,N)/max(Nc,1));
    Hsolve = H + repmat(eye(Nc),1,1,N).*reshape(Hridge_rel*hscale,1,1,N);
    HinvJT = pagemldivide(Hsolve,permute(Jac,[2 1 3]));
    Cphase = pagemtimes(Jac,HinvJT);
    chi_page = chi_page + Cphase.*reshape(p_ip,1,1,N);
end
chi_page = 0.5*(chi_page + permute(chi_page,[2 1 3]));
end


function [lam,is_exact] = Min_Eig_Sym_Pages(H)
% Fast exact eigenvalue for Nc<=2; use pageeig when available for Nc>2.
Nc = size(H,1);
N  = size(H,3);
is_exact = true;
if Nc == 1
    lam = reshape(H(1,1,:),1,N);
elseif Nc == 2
    a = reshape(H(1,1,:),1,N);
    b = reshape(H(1,2,:),1,N);
    d = reshape(H(2,2,:),1,N);
    lam = 0.5*(a+d-sqrt((a-d).^2+4*b.^2));
else
    try
        [~,D] = pageeig(H);
        eigvals = zeros(Nc,N);
        for ic = 1:Nc
            eigvals(ic,:) = reshape(real(D(ic,ic,:)),1,N);
        end
        lam = min(eigvals,[],1);
    catch
        % No node loop fallback: return a vectorized Gershgorin lower bound.
        % It is conservative for detecting possible negative curvature but is
        % not an exact eigenvalue; DIAG.lambda_exact reports this condition.
        lower = zeros(Nc,N);
        for ic = 1:Nc
            center = reshape(H(ic,ic,:),1,N);
            rowabs = reshape(sum(abs(H(ic,:,:)),2),1,N);
            lower(ic,:) = center - (rowabs-abs(center));
        end
        lam = min(lower,[],1);
        is_exact = false;
    end
end
end


function [chi_solver,shift] = Condition_Chi_Pages(chi_phys,NUM)
% Conservative PSD conditioning for optional solver diagnostics. Uses a
% Gershgorin lower bound and therefore never needs a grid-node loop.
Ne = size(chi_phys,1);
N  = size(chi_phys,3);
floor_rel = GetControl(NUM,'path_chi_floor_rel',1e-12);
scale = max(1,reshape(sqrt(sum(chi_phys.^2,[1 2])),1,N)/max(Ne,1));
lower = zeros(Ne,N);
for ie = 1:Ne
    center = reshape(chi_phys(ie,ie,:),1,N);
    rowabs = reshape(sum(abs(chi_phys(ie,:,:)),2),1,N);
    lower(ie,:) = center - (rowabs-abs(center));
end
lb = min(lower,[],1);
shift = max(0,floor_rel*scale-lb);
chi_solver = chi_phys + repmat(eye(Ne),1,1,N).*reshape(shift,1,1,N);
chi_solver = 0.5*(chi_solver + permute(chi_solver,[2 1 3]));
end


function dz = Page_Newton_Direction(Js,rs,ridge_rel)
% Batched damped Newton direction.  A tiny ridge protects nearly singular
% pages without changing the residual whose root is tested by line search.
Nz = size(Js,1);
N  = size(Js,3);
scale = max(1,reshape(sqrt(sum(Js.^2,[1 2])),1,N)/max(Nz,1));
ridge = ridge_rel*scale;
A = Js + repmat(eye(Nz),1,1,N).*reshape(ridge,1,1,N);
dz3 = pagemldivide(A,-reshape(rs,Nz,1,N));
dz = reshape(dz3,Nz,N);

bad = any(~isfinite(dz),1);
if any(bad)
    % Vectorized LM rescue on only failed pages, still no spatial loop.
    Jb = Js(:,:,bad);
    rb = rs(:,bad);
    Nb = sum(bad);
    JT = permute(Jb,[2 1 3]);
    lam = max(1e-10,ridge(bad));
    A2 = pagemtimes(JT,Jb) + repmat(eye(Nz),1,1,Nb).*reshape(lam,1,1,Nb);
    b2 = -pagemtimes(JT,reshape(rb,Nz,1,Nb));
    dz(:,bad) = reshape(pagemldivide(A2,b2),Nz,Nb);
end
end


function [rs,Js] = Scale_Residual_Jacobian(r,J,layout,mu_scale,E_scale)
rs = Scale_Residual_Only(r,layout,mu_scale,E_scale);
Js = J;
if layout.Nc_tot > 0
    Js(1:layout.Nc_tot,:,:) = Js(1:layout.Nc_tot,:,:) ./ reshape(mu_scale,1,1,layout.N);
end
Js(layout.E_ids,:,:) = Js(layout.E_ids,:,:) ./ reshape(E_scale,1,1,layout.N);
end


function rs = Scale_Residual_Only(r,layout,mu_scale,E_scale)
rs = r;
if layout.Nc_tot > 0
    rs(1:layout.Nc_tot,:) = rs(1:layout.Nc_tot,:) ./ mu_scale;
end
rs(layout.E_ids,:) = rs(layout.E_ids,:) ./ E_scale;
end


function dz = Limit_Step_Batch(dz,layout,dc_max,dmu_max)
N = size(dz,2);
fac = ones(1,N);
if layout.Nc_tot > 0 && ~isempty(dc_max) && isfinite(dc_max) && dc_max > 0
    dmax = max(abs(dz(1:layout.Nc_tot,:)),[],1);
    fac = min(fac,dc_max./max(dmax,eps));
end
if ~isempty(dmu_max) && isfinite(dmu_max) && dmu_max > 0
    mmax = max(abs(dz(layout.mu_ids,:)),[],1);
    fac = min(fac,dmu_max./max(mmax,eps));
end
dz = dz.*min(1,fac);
end


function a = Fraction_To_Boundary_Batch(z,dz,layout,c_floor)
% Keep independent and dependent endmember fractions inside the simplex.
N   = layout.N;
tau = 0.995;
a   = ones(1,N);
for ip = 1:layout.Np
    ids = layout.c_ids{ip};
    if isempty(ids); continue; end
    cv = z(ids,:);
    dv = dz(ids,:);
    rat = inf(size(cv));
    neg = dv < 0;
    rat(neg) = (cv(neg)-c_floor)./(-dv(neg));
    a = min(a,tau*min(rat,[],1));
    cdep = 1-sum(cv,1);
    ddep = -sum(dv,1);
    negdep = ddep < 0;
    ratdep = inf(1,N);
    ratdep(negdep) = (cdep(negdep)-c_floor)./(-ddep(negdep));
    a = min(a,tau*ratdep);
end
a = max(min(a,1),0);
end


function [z,layout] = Pack_Unknown_Batch(c,mu)
Np = numel(c);
Ne = size(mu,1);
N  = size(mu,2);
layout.Np = Np;
layout.Ne = Ne;
layout.N  = N;
layout.Nc = zeros(1,Np);
layout.c_ids = cell(1,Np);
counter = 0;
for ip = 1:Np
    layout.Nc(ip) = numel(c{ip});
    layout.c_ids{ip} = counter + (1:layout.Nc(ip));
    counter = counter + layout.Nc(ip);
end
layout.Nc_tot = counter;
layout.mu_ids = counter + (1:Ne);
layout.E_ids  = layout.mu_ids;
layout.Nz     = counter + Ne;
z = zeros(layout.Nz,N);
for ip = 1:Np
    for ic = 1:layout.Nc(ip)
        z(layout.c_ids{ip}(ic),:) = reshape(c{ip}{ic},1,N);
    end
end
z(layout.mu_ids,:) = mu;
end


function [c,mu] = Unpack_Unknown_Batch(z,layout,c)
for ip = 1:layout.Np
    for ic = 1:layout.Nc(ip)
        c{ip}{ic} = reshape(z(layout.c_ids{ip}(ic),:),1,layout.N);
    end
end
mu = z(layout.mu_ids,:);
end


function ci = Cell_From_Z_Phase(z,layout,ip)
ids = layout.c_ids{ip};
ci = cell(1,numel(ids));
for ic = 1:numel(ids)
    ci{ic} = reshape(z(ids(ic),:),1,layout.N);
end
end


function value = GetControl(NUM,name,default_value)
if isstruct(NUM) && isfield(NUM,name) && ~isempty(NUM.(name))
    value = NUM.(name);
else
    value = default_value;
end
end
