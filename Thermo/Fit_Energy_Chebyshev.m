function [f, mu, S] = Fit_Energy_Chebyshev(convex_yes,g, T, X)
% Fit energy using multidimensional Chebyshev polynomials
%Order of energy function
ORD           =  4; 
Nc            =  length(X)-1;

%Subtract configurational entropy (will need to be added back)
X_Fit         = [];
g0            =  0;
for i = 1:length(X)
    g0        =  g0 + 8.314*T(1).*(X{i}.*logo(X{i}));
    X_Fit     = [X_Fit,X{i}(:)];
end

%Select valid points
id            =  find( ((imag(g)<1e-10)==1).*(g<1e10)  );
X_Fit         =  X_Fit(id,1:end-1);
Y_Fit         =  g(id) - g0(id);
id            =  isfinite(Y_Fit);
X_Fit         =  X_Fit(id,:);
Y_Fit         =  real(Y_Fit(id));

%Map to [-1,1] for Chebyshev
Xmin          =  min(X_Fit);
Xmax          =  max(X_Fit);
Xscaled       =  2*(X_Fit - Xmin)./(Xmax - Xmin) - 1;

%Build Chebyshev basis
%Precompute Chebyshev values for each dimension
terms         =  multiIndex(Nc,ORD); 
Nterms        =  size(terms,1);
Npts          = size(Xscaled,1);
ChebVals      = cell(1,Nc);
for j = 1:Nc
    maxOrd = max(terms(:,j));
    ChebVals{j} = zeros(Npts, maxOrd+1);
    % recurrence: T_0(x) = 1, T_1(x) = x, T_{n+1}(x) = 2x*T_n(x) - T_{n-1}(x)
    ChebVals{j}(:,1) = 1;
    ChebVals{j}(:,2) = Xscaled(:,j);
    for n = 2:maxOrd
        ChebVals{j}(:,n+1) = 2*Xscaled(:,j).*ChebVals{j}(:,n) - ChebVals{j}(:,n-1);
    end
end

% Now build A
A = ones(Npts, Nterms);
for k = 1:Nterms
    for j = 1:Nc
        A(:,k) = A(:,k) .* ChebVals{j}(:, terms(k,j)+1);
    end
end

%Solve least squares
coeff  = A\Y_Fit;

%Construct symbolic energy function
c      = sym('c',[1,length(X)]);
c(end) = 1 - sum(c(1:end-1));
f_sym  = 0;

%Add back entropy
for i = 1:length(X)
    f_sym = f_sym + 8.314*T(1)*c(i)*log(1e-32+c(i));
end

%Map symbolic c to Chebyshev domain
c_scaled = 2*(c(1:end-1) - Xmin) ./ (Xmax - Xmin) - 1;

%Add Chebyshev polynomial term
for k = 1:Nterms
    prodTerm = 1;
    for j = 1:Nc
        prodTerm = prodTerm * chebyshevT(terms(k,j), c_scaled(j));
    end
    f_sym = simplify(f_sym + coeff(k)*prodTerm);
end

%Chemical potential
c             =  c(1:end-1);
mu            =  cell(1,Nc);
for i = 1:Nc
    mu_fun    =  matlabFunction(diff(f_sym,c(i)));
    mu{i}     =  @(c)  mu_fun(c{:});
end
mu            =  @(c) cellfun(@(f) f(c), mu, 'UniformOutput', false);

%Hessian and inverse
H             =  simplify(hessian(f_sym, c(1:Nc)));
Hfun          =  matlabFunction(H,'Vars',num2cell(c));
S             =  @(C) localHinv_component(Hfun,C);

%Energy function
f             =  matlabFunction(f_sym);
f             =  @(c)   f(c{:});

%Convex hull
if convex_yes==1
    %Penalty
    lambda     =  1e9;
    %Hull
    K          =  convhulln(X_Fit,{'QJ'});      % indices of points forming each facet
    hullV      =  unique(K(:));                 % unique indices of vertices on the hull
    V          =  X_Fit(hullV,:);               % convex hull vertices
    % inside     =  inhull(X_test,X_Fit,[],1e-8);
    
    %Prepare convex hull
    facets     =  build_hull_facets(X_Fit);
    % [facets,proj,meanX] = build_hull_facets(X_Fit);
    %Make function handle
    fc         = @(c) hull_penalty_value(c, V, facets, lambda);
    muc        = @(c) hull_grad_value(c, V, facets, lambda);
    Sc         = @(c) hull_sinv_value(c, V, facets, lambda);
    %Add on f
    f_tot      = @(c) real(f(c)+fc(c));
    %Add on mu
    mu_tot     = @(c) cellfun(@(a,b) a + b, mu(c), muc(c), 'UniformOutput', false);
    %Add on mu
    S_tot      = @(c) cellfun(@(a,b) a + b, S(c) , Sc(c) , 'UniformOutput', false);
    %Replace
    f          = f_tot;
    mu         = mu_tot;
    S          = S_tot;

    fc({X_Fit(:,1),X_Fit(:,2),X_Fit(:,3),X_Fit(:,4)})

    %Display
    disp('Choose to use convex hull')
    disp(['Number of vertices: ',num2str(size(V,1))])
    disp(['Number of facets: ',num2str(size(facets.n,1))])
end


disp('Chebyshev fit done.');
end

%-------------------------
function logo = logo(x)
x(x<=0) = x(x<=0) + 1e-12;
x(x>=1) = x(x>=1) - 1e-12;
logo = log(x);
end


%-------------------------
function terms = multiIndex(d, order)
% Generate all multi-indices for d variables up to total order
terms = [];
for total = 0:order
    tmp = nchoosek(repmat(0:total,1,d), d);
    tmp = unique(tmp(sum(tmp,2)==total,:), 'rows');
    terms = [terms; tmp];
end
end



%-------------------------
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




%% Penalty
function P = hull_penalty_value(ccell,V,facets,lambda)
% Vectorized convex hull penalty
% Inputs:
%   ccell  : 1 x d cell array, each cell is matrix of any size
%   V      : n x d convex hull vertices
%   lambda : penalty strength
% Output:
%   P      : penalty values (same size as ccell{1})

d  = numel(ccell);
sz = size(ccell{1});
N  = numel(ccell{1});

% Flatten input
X = zeros(N,d);
for i = 1:d
    X(:,i) = ccell{i}(:);
end

% Convex hull facets
Svals   = X*facets.n' - facets.b';
outside = max(Svals,[],2) > 1e-6;

%Initialize P
P       = zeros(N,1);

if any(outside)
    Xout  = X(outside,:);
    n_out = sum(outside);

    % Closest vertex per outside point
    diffs = permute(V,[1,3,2]) - reshape(Xout,1,n_out,d);
    dist2 = sum(diffs.^2,3);
    [~, idx] = min(dist2,[],1);
    nearest  = V(idx,:);

    delta = Xout - nearest;
    P(outside) = lambda * sum(delta.^2,2);
end

% Reshape to original input size
P = reshape(P, sz);
end







%% Grad
function mu = hull_grad_value(ccell,V,facets, lambda)
% Gradient of convex hull penalty
d  = numel(ccell);
sz = size(ccell{1});
N  = numel(ccell{1});

% Flatten input
X = zeros(N,d);
for i = 1:d
    X(:,i) = ccell{i}(:);
end

% Convex hull facets
Svals   = X*facets.n' - facets.b';
outside = max(Svals,[],2) > 1e-10;

%Initialize mu
mu      = cell(1,d);
for i = 1:d
    mu{i} = zeros(sz);
end

if any(outside)
    Xout = X(outside,:);
    n_out = sum(outside);

    diffs = permute(V,[1,3,2]) - reshape(Xout,1,n_out,d);
    dist2 = sum(diffs.^2,3);
    [~, idx] = min(dist2,[],1);
    nearest = V(idx,:);

    delta = Xout - nearest;

    for i = 1:d
        temp = zeros(N,1);
        temp(outside) = 2*lambda*delta(:,i);
        mu{i} = reshape(temp, sz);
    end
end
end





%% Sinv
function S = hull_sinv_value(ccell,V,facets,lambda)
% Inverse Hessian of convex hull penalty (diagonal)
d  = numel(ccell);
sz = size(ccell{1});
N  = numel(ccell{1});

% Flatten input
X = zeros(N,d);
for i = 1:d
    X(:,i) = ccell{i}(:);
end

% Convex hull facets
Svals   = X*facets.n' - facets.b';
outside = max(Svals,[],2) > 1e-10;

%Initialize S
S       = cell(d,d);
for i = 1:d
    for j = 1:d
        temp = zeros(N,1);
        if i==j && any(outside)
            temp(outside) = 1/(2*lambda);
        end
        S{i,j} = reshape(temp, sz);
    end
end
end



function facets = build_hull_facets(X)
% Precompute hull facet equations from full dataset X
% facets: struct with fields n (F x d) normals, b (F x 1) offsets

K = convhulln(X);
F = size(K,1);
d = size(X,2);

n = zeros(F,d);
b = zeros(F,1);

for jf = 1:F
    verts = X(K(jf,:),:);
    Mdiff = (verts(2:end,:) - verts(1,:))';
    nvec  = null(Mdiff.');
    if isempty(nvec)
        nvec = cross(verts(2,:) - verts(1,:), ...
                     verts(3,:) - verts(1,:)).';
    end
    nn = nvec(:,1); nn = nn/norm(nn);

    % orient inward
    if (nn' * mean(X,1)' - nn' * verts(1,:)') > 0
        nn = -nn;
    end

    n(jf,:) = nn;
    b(jf)   = nn'*verts(1,:)';
end

facets.n = n;   % facet normals
facets.b = b;   % offsets
end


function outside = hull_outside_mask(X, facets)
tol = 1e-10;
Svals = X*facets.n' - facets.b';
outside = max(Svals,[],2) > tol;
end
