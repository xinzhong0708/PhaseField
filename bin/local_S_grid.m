function Spg = local_S_grid(F, JacFun, HmapFun, ccell)
%Dimensions and flattening
d     = numel(ccell);
N     = numel(ccell{1});
Xcols = cellfun(@(M) M(:).', ccell, 'uni', 0);

%Obtain mu_e
X = zeros(N,d);
for i = 1:d, X(:,i) = ccell{i}(:); end
cols = mat2cell(X, size(X,1), ones(1,d));
mu_e = cell2mat(F.mu_e(cols))';              % r×N

%J: r×d×N
Jflat    = JacFun(Xcols{:});                 % (r*d)×N
[Mrow,~] = size(Jflat);
r        = Mrow / d;                         % number of independent element
Jpages   = reshape(Jflat, r, d, N);

%H_c: d×d×N
Hc       = F.H_c(ccell);

%Curvature H^{(k)}: each d×d×N
Hk       = cell(1,r);
for k = 1:r
    Hkflat = HmapFun{k}(Xcols{:});          % (d*d)×N
    Hk{k}  = reshape(Hkflat, d, d, N);
end

% Assemble A = Hc - Σ mu_e(k,:)*H^{(k)} (pagewise)
A = Hc;
for k = 1:r
    wk = reshape(mu_e(k,:),1,1,N);
    A  = A - Hk{k} .* wk;
end

% S = J * A^{-1} * J^T without per-page loop
Jt  = permute(Jpages,[2 1 3]);          % d×r×N
X   = pagemldivide(A, Jt);              % d×r×N
Spg = pagemtimes(Jpages, X);            % r×r×N

end