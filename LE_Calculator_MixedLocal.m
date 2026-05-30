function [c,mu_e,chi,DIAG] = LE_Calculator_MixedLocal(pars,p,c,E,mu_e,eta,level)
%LE_CALCULATOR_MIXEDLOCAL Local mixed [dmu; dc] LE calculator.
%
% Minimal local implementation of option 1:
% solve dmu and dc together instead of eliminating dc first.
%
% Unknown per node:
%   x = [dmu_e ; dc_1 ; dc_2 ; ...]
%
% Equations:
%   (1/eta)dmu + sum_a p_a J_a dc_a = E - sum_a p_a e_a - mu/eta
%   H_a dc_a - J_a' dmu = -(mu_c,a - J_a' mu)
%
% This keeps negative H in the block system. No positive Hessian
% regularization is applied.
%
% level = [alpha, max_iter, h_shift, chi_mode]
%   alpha    : damping, e.g. 0.2
%   max_iter : e.g. 50-200
%   h_shift  : optional tiny diagonal shift in H, default 0
%   chi_mode : 0 signed output chi, 1 positive-clipped output chi

Np = numel(c);
Ne = numel(E);
N  = numel(E{1});

alpha = level(1);
Miter = level(2);

if numel(level) >= 3 && ~isempty(level(3))
    h_shift = level(3);
else
    h_shift = 0;
end

if numel(level) >= 4 && ~isempty(level(4))
    chi_mode = level(4);
else
    chi_mode = 0;
end

tol_dx  = 1e-10;
tol_res = 1e-8;

if isscalar(eta)
    eta_vec = eta*ones(1,N);
else
    eta_vec = eta(:).';
end

E_mat  = cell2mat(E(:));
mu_mat = cell2mat(mu_e(:));

DIAG.res_hist = nan(1,Miter);
DIAG.dx_hist  = nan(1,Miter);
DIAG.dmu_hist = nan(1,Miter);

for it = 1:Miter

    [dmu,dc,STEP] = MixedLocalStep(pars,p,c,E_mat,mu_mat,eta_vec,h_shift);

    dxmax  = max(MaxAbsCell(dc),max(abs(dmu(:))));
    resmax = MixedResidualNorm(pars,p,c,E_mat,mu_mat,eta_vec);

    DIAG.res_hist(it) = resmax;
    DIAG.dx_hist(it)  = dxmax;
    DIAG.dmu_hist(it) = max(abs(dmu(:)));

    if dxmax < tol_dx || resmax < tol_res
        break
    end

    mu_mat = mu_mat + alpha*dmu;
    c      = AddCellStep(c,dc,alpha);

    if STEP.has_nonfinite
        break
    end
end

[dmu,dc,STEP] = MixedLocalStep(pars,p,c,E_mat,mu_mat,eta_vec,h_shift);
mu_mat = mu_mat + alpha*dmu;
c      = AddCellStep(c,dc,alpha);

mu_e = cell(1,Ne);
for ie = 1:Ne
    mu_e{ie} = reshape(mu_mat(ie,:),size(E{ie}));
end

chi_page = EffectiveChi(pars,p,c,eta_vec,chi_mode);

chi = cell(Ne,Ne);
for ie = 1:Ne
    for je = 1:Ne
        chi{ie,je} = reshape(chi_page(ie,je,:),size(E{1}));
    end
end

Cdiag = ChiDiagnostics(chi_page);

DIAG.iter               = it;
DIAG.max_residual       = MixedResidualNorm(pars,p,c,E_mat,mu_mat,eta_vec);
DIAG.max_dmu            = max(abs(dmu(:)));
DIAG.max_dc             = MaxAbsCell(dc);
DIAG.has_nonfinite      = STEP.has_nonfinite;
DIAG.h_shift            = h_shift;
DIAG.chi_mode           = chi_mode;
DIAG.min_chi_eig        = Cdiag.min_chi_eig;
DIAG.negative_chi_nodes = Cdiag.negative_chi_nodes;
DIAG.max_abs_chi        = Cdiag.max_abs_chi;
DIAG.has_nonfinite_chi  = Cdiag.has_nonfinite_chi;
DIAG.res_hist           = DIAG.res_hist(1:it);
DIAG.dx_hist            = DIAG.dx_hist(1:it);
DIAG.dmu_hist           = DIAG.dmu_hist(1:it);

end


function [dmu,dc_all,DIAG] = MixedLocalStep(pars,p,c,E_mat,mu_mat,eta_vec,h_shift)

Np = numel(c);
Ne = size(E_mat,1);
N  = size(E_mat,2);

R = cell(1,Np);
Nc_phase = zeros(1,Np);

for ip = 1:Np
    R{ip} = PhaseThermo(pars{ip},c{ip});
    if ~isempty(R{ip}.H_c) && ~isempty(R{ip}.Jac)
        Nc_phase(ip) = size(R{ip}.H_c,1);
    end
end

K = Ne + sum(Nc_phase);

A = zeros(K,K,N);
b = zeros(K,N);

Emix = zeros(Ne,N);

Ipg = repmat(eye(Ne),1,1,N);
A(1:Ne,1:Ne,:) = Ipg .* reshape(1./eta_vec,1,1,N);

col0 = Ne;

for ip = 1:Np

    e_ip = cell2mat(R{ip}.e(:));
    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + e_ip.*p_ip;

    Nc = Nc_phase(ip);
    if Nc == 0
        continue
    end

    cols = col0 + (1:Nc);
    J = R{ip}.Jac;

    A(1:Ne,cols,:) = J .* reshape(p_ip,1,1,N);

    col0 = col0 + Nc;
end

b(1:Ne,:) = E_mat - Emix - mu_mat./reshape(eta_vec,1,N);

row0 = Ne;
col0 = Ne;

for ip = 1:Np

    Nc = Nc_phase(ip);
    if Nc == 0
        continue
    end

    rows = row0 + (1:Nc);
    cols = col0 + (1:Nc);

    J    = R{ip}.Jac;
    H    = 0.5*(R{ip}.H_c + permute(R{ip}.H_c,[2 1 3]));
    mu_c = cell2mat(R{ip}.mu_c(:));

    if h_shift ~= 0
        H = H + repmat(eye(Nc),1,1,N).*h_shift;
    end

    JTmu3 = pagemtimes(permute(J,[2 1 3]),reshape(mu_mat,Ne,1,N));
    JTmu  = reshape(JTmu3,Nc,N);
    r     = mu_c - JTmu;

    A(rows,1:Ne,:) = -permute(J,[2 1 3]);
    A(rows,cols,:) = H;
    b(rows,:)      = -r;

    row0 = row0 + Nc;
    col0 = col0 + Nc;
end

x3 = pagemldivide(A,reshape(b,K,1,N));
x  = reshape(x3,K,N);

DIAG.has_nonfinite = any(~isfinite(x(:))) || any(~isfinite(A(:))) || any(~isfinite(b(:)));

dmu = x(1:Ne,:);

dc_all = c;
for ip = 1:Np
    for ic = 1:numel(c{ip})
        dc_all{ip}{ic} = zeros(size(c{ip}{ic}));
    end
end

row0 = Ne;
for ip = 1:Np

    Nc = Nc_phase(ip);
    if Nc == 0
        continue
    end

    rows = row0 + (1:Nc);
    dc_ip = x(rows,:);

    for ic = 1:numel(c{ip})
        dc_all{ip}{ic} = reshape(dc_ip(ic,:),size(c{ip}{ic}));
    end

    row0 = row0 + Nc;
end

end


function resmax = MixedResidualNorm(pars,p,c,E_mat,mu_mat,eta_vec)

Np = numel(c);
Ne = size(E_mat,1);
N  = size(E_mat,2);

Emix = zeros(Ne,N);
res_node = zeros(1,N);

for ip = 1:Np

    R = PhaseThermo(pars{ip},c{ip});
    e_ip = cell2mat(R.e(:));
    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + e_ip.*p_ip;

    if ~isempty(R.mu_c) && ~isempty(R.Jac)
        mu_c = cell2mat(R.mu_c(:));
        JTmu3 = pagemtimes(permute(R.Jac,[2 1 3]),reshape(mu_mat,Ne,1,N));
        JTmu  = reshape(JTmu3,size(mu_c,1),N);
        res_node = max(res_node,max(abs(mu_c-JTmu),[],1));
    end
end

mass = max(abs(E_mat - Emix - mu_mat./reshape(eta_vec,1,N)),[],1);
res_node = max(res_node,mass);
res_node(any(~isfinite(mu_mat),1)) = inf;

resmax = max(res_node);

end


function chi_page = EffectiveChi(pars,p,c,eta_vec,chi_mode)

Np = numel(c);
N  = numel(eta_vec);

R0 = PhaseThermo(pars{1},c{1});
Ne = numel(R0.e);

chi_page = repmat(eye(Ne),1,1,N).*reshape(1./eta_vec,1,1,N);

for ip = 1:Np

    R = PhaseThermo(pars{ip},c{ip});
    if isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    H = 0.5*(R.H_c + permute(R.H_c,[2 1 3]));
    J = R.Jac;

    Hinv_JT = pagemldivide(H,permute(J,[2 1 3]));
    S = pagemtimes(J,Hinv_JT);

    if chi_mode == 1
        S = PositiveClipPages(S);
    end

    p_ip = reshape(p(:,:,ip),1,1,N);
    chi_page = chi_page + S.*p_ip;
end

end


function S = PositiveClipPages(S)

[m,~,N] = size(S);

for q = 1:N
    A = 0.5*(S(:,:,q)+S(:,:,q).');
    [V,D] = eig(A);
    lam = diag(D);
    scale = max(1,max(abs(lam)));
    floor_val = 1e-10*scale;
    lam(lam < floor_val) = floor_val;
    S(:,:,q) = V*diag(lam)*V.';
end

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

for q = 1:N
    A = 0.5*(chi_page(:,:,q)+chi_page(:,:,q).');
    if any(~isfinite(A(:)))
        bad = true;
        continue
    end
    ev = eig(A);
    mineig(q) = min(ev);
    maxabs = max(maxabs,max(abs(A(:))));
end

valid = mineig(isfinite(mineig));

if isempty(valid)
    D.min_chi_eig = nan;
else
    D.min_chi_eig = min(valid);
end

D.negative_chi_nodes = nnz(mineig < 0);
D.max_abs_chi = maxabs;
D.has_nonfinite_chi = bad;

end
