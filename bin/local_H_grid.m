function Hpg = local_H_grid(F,JacFun,HmapFun,ccell,lambda)
% H_cell: r×r cell, each ny×nx, with H = ∂μ_e/∂E (pagewise inverse of S)
% lambda (optional): small Tikhonov for robustness (e.g., 1e-12)

if nargin < 5, lambda = 0; end

% --- dims & flattening ---
d     = numel(ccell);
sz    = size(ccell{1});
N     = numel(ccell{1});
Xcols = cellfun(@(M) M(:).', ccell, 'uni', 0);

% --- μ_e (use its size to define r robustly) ---
X = zeros(N,d);
for i = 1:d, X(:,i) = ccell{i}(:); end
cols  = mat2cell(X, size(X,1), ones(1,d));
mu_e  = cell2mat(F.mu_e(cols))';           % r×N
r     = size(mu_e,1);

% --- J: r×d×N (from JacFun) ---
Jflat  = JacFun(Xcols{:});                 % (r*d)×N
Jpages = reshape(Jflat, r, d, N);

% --- H_c: d×d×N ---
Hc     = F.H_c(ccell);

% --- curvature tensors H^{(k)}: d×d×N each ---
Hk = cell(1,r);
for k = 1:r
    Hkflat = HmapFun{k}(Xcols{:});         % (d*d)×N
    Hk{k}  = reshape(Hkflat, d, d, N);
end

% --- A = Hc - Σ_k μ_e(k,:)*H^{(k)} ---
A = Hc;
for k = 1:r
    A = A - Hk{k} .* reshape(mu_e(k,:),1,1,N);
end

% --- S = J * A^{-1} * J^T (same as your local_S_grid) ---
Jt  = permute(Jpages,[2 1 3]);            % d×r×N
Xsol= pagemldivide(A, Jt);                % d×r×N   (A\J^T)
Spg = pagemtimes(Jpages, Xsol);           % r×r×N

% --- H = S^{-1} (pagewise), with optional λI regularization ---
if lambda > 0
    I = eye(r); I = I(:,:,ones(1,N));
    Spg = Spg + lambda * I;
end
I = eye(r); I = I(:,:,ones(1,N));
Hpg = pagemldivide(Spg, I);               % r×r×N

% (optional) symmetrize tiny asymmetries
Hpg = 0.5*(Hpg + permute(Hpg,[2 1 3]));

end
