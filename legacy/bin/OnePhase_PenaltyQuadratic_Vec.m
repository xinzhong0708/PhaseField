function [c_new,mu_e,chi] = OnePhase_PenaltyQuadratic_Vec(Fph, c_old, Etar, eta, reg, alpha, Miter, tol)

% Hyperparameters
if nargin < 4 || isempty(eta),   eta   = 1000;  end
if nargin < 5 || isempty(reg),   reg   = 1e-8;  end
if nargin < 6 || isempty(alpha), alpha = 1.0;   end
if nargin < 7 || isempty(Miter), Miter = 1;     end
if nargin < 8 || isempty(tol),   tol   = 1e-7;  end
pars  = Fph;

% sizes
Ne    = size(Etar,1);
N     = size(Etar,2);

% eta as 1 x N
eta_vec = normalize_eta_vec(eta, size(Etar,2), N);

% initialize solution
c_new   = c_old;

for it = 1:Miter

    % Evaluate at current c
    R           = PhaseThermo(pars, c_new);
    e_ref       = cell2mat(R.e(:));                    % Ne x N

    % No internal composition DOF (e.g. pure quartz)
    mu_c        = cell2mat(R.mu_c(:));                 % Nc x N
    Hc          = R.H_c;                               % Nc x Nc x N
    J           = R.Jac;                               % Ne x Nc x N

    Nc          = size(mu_c,1);
    IpagesC     = repmat(eye(Nc), 1, 1, N);
    IpagesE     = repmat(eye(Ne), 1, 1, N);

    % Regularize Hessian in c-space
    Hreg        = Hc + reg * IpagesC;

    % Build C = J H^{-1} J^T and B = e - J H^{-1} mu_c
    JT          = permute(J, [2 1 3]);                 % Nc x Ne x N
    Hinv_JT     = pagemldivide(Hreg, JT);              % Nc x Ne x N
    C_page      = pagemtimes(J, Hinv_JT);              % Ne x Ne x N

    Hinv_mu3    = pagemldivide(Hreg, reshape(mu_c, Nc, 1, N));
    Hinv_mu     = reshape(Hinv_mu3, Nc, N);            % Nc x N
    JHinvmu3    = pagemtimes(J, reshape(Hinv_mu, Nc, 1, N));
    JHinvmu     = reshape(JHinvmu3, Ne, N);            % Ne x N

    B_mat       = e_ref - JHinvmu;                     % Ne x N

    % Effective one-phase susceptibility
    chi_page    = C_page + IpagesE .* reshape(1 ./ eta_vec, 1, 1, N);

    % Solve for elemental diffusion potential:
    % (I/eta + C) mu_e = E_tar - B
    rhs_mu      = Etar - B_mat;                        % Ne x N
    mu3         = pagemldivide(chi_page, reshape(rhs_mu, Ne, 1, N));
    mu_mat      = reshape(mu3, Ne, N);                 % Ne x N

    % One-phase tangent update:
    % dc = H^{-1}(J^T mu_e - mu_c)
    JTmu3       = pagemtimes(JT, reshape(mu_mat, Ne, 1, N));
    JTmu        = reshape(JTmu3, Nc, N);               % Nc x N
    try
        rhs_dc  = JTmu - mu_c;  % For full phase
    catch
        rhs_dc  = [];           % For pure endmembers
    end
    dc3         = pagemldivide(Hreg, reshape(rhs_dc, Nc, 1, N));
    dc          = reshape(dc3, Nc, N);                 % Nc x N

    % Damped update
    c_old_it    = c_new;
    for ic = 1:Nc
        c_new{ic} = c_new{ic} + alpha * reshape(dc(ic,:), size(c_new{ic}));
    end

    % Check change of c
    cchg = 0;
    for ic = 1:Nc
        cchg = max(cchg, max(abs(c_new{ic}(:) - c_old_it{ic}(:))));
    end
    if cchg < tol
        break
    end
end

% Recompute at final updated state and return effective mu_e / chi_eff
R           = PhaseThermo(pars, c_new);
e_ref       = cell2mat(R.e(:));                        % Ne x N
has_dof     = ~(isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac));

IpagesE     = repmat(eye(Ne), 1, 1, N);

if ~has_dof
    % Pure rigid phase case
    chi_page = IpagesE .* reshape(1 ./ eta_vec, 1, 1, N);
    rhs_mu   = Etar - e_ref;
    mu3      = pagemldivide(chi_page, reshape(rhs_mu, Ne, 1, N));
    mu_mat   = reshape(mu3, Ne, N);
else
    mu_c     = cell2mat(R.mu_c(:));                    % Nc x N
    Hc       = R.H_c;                                  % Nc x Nc x N
    J        = R.Jac;                                  % Ne x Nc x N
    Nc       = size(mu_c,1);

    IpagesC  = repmat(eye(Nc), 1, 1, N);
    Hreg     = Hc + reg * IpagesC;

    JT       = permute(J, [2 1 3]);                    % Nc x Ne x N
    Hinv_JT  = pagemldivide(Hreg, JT);                 % Nc x Ne x N
    C_page   = pagemtimes(J, Hinv_JT);                 % Ne x Ne x N

    Hinv_mu3 = pagemldivide(Hreg, reshape(mu_c, Nc, 1, N));
    Hinv_mu  = reshape(Hinv_mu3, Nc, N);               % Nc x N
    JHinvmu3 = pagemtimes(J, reshape(Hinv_mu, Nc, 1, N));
    JHinvmu  = reshape(JHinvmu3, Ne, N);               % Ne x N

    B_mat    = e_ref - JHinvmu;                        % Ne x N

    chi_page = C_page + IpagesE .* reshape(1 ./ eta_vec, 1, 1, N);
    rhs_mu   = Etar - B_mat;                           % Ne x N

    mu3      = pagemldivide(chi_page, reshape(rhs_mu, Ne, 1, N));
    mu_mat   = reshape(mu3, Ne, N);
end

% Output as cells in your original style
mu_e = cell(1,Ne);
chi  = cell(Ne,Ne);
for i = 1:Ne
    mu_e{i} = reshape(mu_mat(i,:), size(Etar(i,:)));
    mu_e{i} = reshape(mu_mat(i,:), size(Etar(1,:)));
    for j = 1:Ne
        chi{i,j} = reshape(chi_page(i,j,:), 1, []);
    end
end

% reshape chi and mu_e to match field shape
for i = 1:Ne
    mu_e{i} = reshape(mu_mat(i,:), size(Etar(1,:)));
    for j = 1:Ne
        chi{i,j} = reshape(chi_page(i,j,:), size(Etar(1,:)));
    end
end

end


function eta_vec = normalize_eta_vec(eta, ncol, N)
% Return eta as 1 x N

if isscalar(eta)
    eta_vec = eta * ones(1,N);
    return
end

if isvector(eta) && numel(eta) == N
    eta_vec = reshape(eta, 1, []);
    return
end

if isequal(size(eta), [1 N]) || isequal(size(eta), [N 1])
    eta_vec = reshape(eta, 1, []);
    return
end

if isequal(size(eta), [1 ncol])
    eta_vec = reshape(eta, 1, []);
    return
end

error('normalize_eta_vec: eta has incompatible size.');
end