function [f, mu, S] = Fit_Energy_Chebyshev(g, T, X)
% Fit energy using multidimensional Chebyshev polynomials
% g: energy vector
% X: cell array of composition arrays
% Returns:
%   f: energy function handle
%   mu: cell array of chemical potentials
%   S: inverse Hessian function handle

ORD = 5;  % Chebyshev polynomial order

%% Flatten compositions
X_Fit = [];
g0 = 0;
for i = 1:length(X)
    g0 = g0 + 8.314*T(1)*(X{i}.*logo(X{i}));
    X_Fit = [X_Fit, X{i}(:)];
end

%% Select valid points
id = find(((imag(g)<1e-10)==1).*(g<1e10));
X_Fit = X_Fit(id,1:end-1); % drop last component
Y_Fit = g(id) - g0(id);
id = isfinite(Y_Fit);
X_Fit = X_Fit(id,:);
Y_Fit = real(Y_Fit(id));

%% Map to [-1,1] for Chebyshev
d = size(X_Fit,2);
Xmin = min(X_Fit);
Xmax = max(X_Fit);
Xscaled = 2*(X_Fit - Xmin)./(Xmax - Xmin) - 1;

%% Build Chebyshev basis
terms = multiIndex(d, ORD);  % all multi-indices
Nterms = size(terms,1);
A = zeros(size(Xscaled,1), Nterms);
for k = 1:Nterms
    tmp = ones(size(Xscaled,1),1);
    for j = 1:d
        tmp = tmp .* chebyshevT(terms(k,j), Xscaled(:,j));
    end
    A(:,k) = tmp;
end

%% Solve least squares
coeff = A\Y_Fit;

%% Construct symbolic energy function
c = sym('c',[1,length(X)]);
c(end) = 1 - sum(c(1:end-1));
f_sym = 0;

% Add back entropy
for i = 1:length(X)
    f_sym = f_sym + 8.314*T(1)*c(i)*log(c(i));
end

% Map symbolic c to Chebyshev domain
c_scaled = 2*(c(1:end-1) - Xmin) ./ (Xmax - Xmin) - 1;

% Add Chebyshev polynomial term
for k = 1:Nterms
    prodTerm = 1;
    for j = 1:d
        prodTerm = prodTerm * chebyshevT(terms(k,j), c_scaled(j));
    end
    f_sym = f_sym + coeff(k)*prodTerm;
end

%% Chemical potential
mu = cell(1,d);
for i = 1:d
    mu_fun = matlabFunction(diff(f_sym,c(i)), 'Vars',{c});
    mu{i} = @(cval) mu_fun(cval);
end
mu = @(cval) cellfun(@(f) f(cval), mu, 'UniformOutput', false);

%% Hessian and inverse
H = hessian(f_sym, c(1:d));
Hfun = matlabFunction(H,'Vars',{c});
S = @(C) localHinv_component(Hfun,C);

%% Energy function handle
f_fun = matlabFunction(f_sym,'Vars',{c});
f = @(cval) f_fun(cval);

disp('Chebyshev fit done.');
end

%% -------------------------
function logo = logo(x)
x(x<=0) = x(x<=0) + 1e-12;
x(x>=1) = x(x>=1) - 1e-12;
logo = log(x);
end

%% -------------------------
function terms = multiIndex(d, order)
% Generate all multi-indices for d variables up to total order
terms = [];
for total = 0:order
    tmp = nchoosek(repmat(0:total,1,d), d);
    tmp = unique(tmp(sum(tmp,2)==total,:), 'rows');
    terms = [terms; tmp];
end
end

%% -------------------------
function Hij = localHinv_component(Hfun,C)
% C is a cell array of matrices
if ~iscell(C), error('Input must be a cell array'); end
sz = size(C{1});
nPts = numel(C{1});
Nc = numel(C);

% Flatten inputs
args_flat = cellfun(@(x) x(:), C, 'UniformOutput', false);

% Evaluate Hessians
Hvals = zeros(Nc,Nc,nPts);
for k = 1:nPts
    args_k = cellfun(@(x) x(k), args_flat, 'UniformOutput', false);
    Hvals(:,:,k) = feval(Hfun, args_k{:});
end

% Invert Hessians
I = repmat(eye(Nc),1,1,nPts);
Hinv_all = pagemldivide(Hvals,I);

% Return as cell
for i = 1:Nc
    for j = 1:Nc
        Hij{i,j} = reshape(Hinv_all(i,j,:), sz);
    end
end
end
