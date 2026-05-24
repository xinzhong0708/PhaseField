function mu = localmue_component(mu_fun, Jac, ccell)
%Prepare
d     = numel(ccell);         % # independent endmembers
sz    = size(ccell{1});       % grid shape
N     = numel(ccell{1});      % # points

% Flatten input matrices to N×d
X = zeros(N, d);
for i = 1:d, X(:,i) = ccell{i}(:); end

% Chemical potential in c-space (N x d)
cols  = mat2cell(X, size(X,1), ones(1,d));
mu_c  = mu_fun(cols{:});                      % N x d

% Jacobian J = de/dc, provided flattened as (n_e*d) x N
args  = mat2cell(X, N, ones(1,d));
args  = cellfun(@(v) v.', args, 'UniformOutput', false);   % each 1 x N
Jflat = Jac(args{:});                                      % (n_e*d) x N

% Infer n_e and reshape to n_e x d x N
[M, Ncheck] = size(Jflat);
assert(Ncheck==N, 'Jac returned N mismatch.');
assert(mod(M, d)==0, 'Jac size incompatible with d.');
n_e  = M / d;
Jval = reshape(Jflat, n_e, d, N);                          % n_e x d x N

% Solve J^T * mu_e^T = mu_c^T pagewise  -> mu_e = mu_c / J
JT   = permute(Jval, [2 1 3]);                              % d x n_e x N
gcT  = permute(mu_c, [2 3 1]);                              % d x 1 x N
geT  = pagemldivide(JT, gcT);                               % n_e x 1 x N
mu_e = permute(geT, [3 1 2]);                               % N x n_e

% Pack to cells on the grid
mu = cell(1, n_e);
for i = 1:n_e
    mu{i} = (reshape(mu_e(:,i), sz));
end
end