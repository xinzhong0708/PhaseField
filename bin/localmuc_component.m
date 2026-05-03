function mu = localmuc_component(mu_fun,ccell)
%Prepare
d     = numel(ccell);         % number of components
sz    = size(ccell{1});       % grid shape
%Flatten input matrices to N×d ---
for i = 1:d
    X(:,i) = ccell{i}(:);
end
%Chemical potential
cols  = mat2cell(X, size(X,1), ones(1,d));
mu_c  = mu_fun(cols{:});                      % N x d
for i = 1:d
    mu{i} = (reshape(mu_c(:,i),sz(1),sz(2)));
end
end
