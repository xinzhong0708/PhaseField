clear;clf;hold on

x3         = [0 2 2 0 0 2 2 0 2.1 2*rand(1,5)];
y3         = [0 0 2 2 0 0 2 2 1   2*rand(1,5)];
z3         = [0 0 0 0 2 2 2 2 1.2 2*rand(1,5)];
plot3(x3,y3,z3,'.','markersize',30);axis equal;grid on
X_Fit      = [x3(:), y3(:), z3(:)]; % Nx3 points

%Hull
K          =  convhulln(X_Fit);      % indices of points forming each facet
hullV      =  unique(K(:));          % unique indices of vertices on the hull
V          =  X_Fit(hullV,:);        % convex hull vertices
lambda     =  1e9;

%Prepare convex hull
facets     =  build_hull_facets(X_Fit);

%Take reference
xo         = {[5,rand*3] [4,rand*5] [-7,rand*4]}; % cell input for consistency
f1         =  hull_penalty_value(xo,V,facets,lambda);
mu1        =  hull_grad_value(xo,V,facets,lambda);
S          =  hull_sinv_value(xo,V,facets,lambda);
dx         =  1e-7;
%
for i = 1:3
    for j = 1:3
        x          =  xo;
        x{j}       =  xo{j}+dx;
        f2         =  hull_penalty_value(x,V,facets,lambda);
        mu2        =  hull_grad_value(x,V,facets,lambda);
        S_num{i,j} = (mu2{i} - mu1{i})/dx;
        mu_num{j}  = (f2-f1)/dx;
    end
end
%invert
for l = 1:length(S_num{1})
    M = [S_num{1,1}(l) S_num{1,2}(l) S_num{1,3}(l)
         S_num{2,1}(l) S_num{2,2}(l) S_num{2,3}(l)
         S_num{3,1}(l) S_num{3,2}(l) S_num{3,3}(l)];
    M =  inv(M);
    for i = 1:3
        for j = 1:3
            SS{i,j}(l) = M(i,j);
        end
    end
end

S
SS

hold on
plot3(xo{1}(1),xo{2}(1),xo{3}(1),'r.','markersize',40)
plot3(xo{1}(2),xo{2}(2),xo{3}(2),'k.','markersize',40)

mu

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
outside = max(Svals,[],2) > 1e-10;

%Initialize P
P       = zeros(N,1);

if any(outside)
    Xout = X(outside,:);
    n_out = sum(outside);

    % Closest vertex per outside point
    diffs = permute(V,[1,3,2]) - reshape(Xout,1,n_out,d);
    dist2 = sum(diffs.^2,3);
    [~, idx] = min(dist2,[],1);
    nearest = V(idx,:);

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
