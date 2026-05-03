function [f,mu,S] = Fit_Energy(convex_yes,g,T,X)
%Here neural network is used for fitting
%Subtract configurational entropy (will need to be added back)
X_Fit                      =  zeros(size(X{1},1),length(X));
g0                         =  0;
for i = 1:length(X)
    g0                     =  g0 + 8.314*T(1).*(X{i}.*logo(X{i}));
    X_Fit(:,i)             =  X{i};
end

%Fit the function
id                         =  find( ((imag(g)<1e-10)==1).*(g<1e10)  );
X_Fit                      =  X_Fit(id,1:end-1);
Y_Fit                      =  real(g(id) - g0(id));
id                         =  isfinite(Y_Fit);
X_Fit                      =  X_Fit(id,:);
Y_Fit                      =  Y_Fit(id  );

%Training data
[Xs, inputSettings]        = mapminmax(X_Fit', -1, 1);   % scale inputs
[Ys, outputSettings]       = mapminmax(Y_Fit', -1, 1);   % scale outputs

%Training neural network
hiddenLayerSize            = [35,35];
net                        =  fitnet(hiddenLayerSize,'trainlm');
net.trainParam.epochs      =  1000;
net.trainParam.goal        =  1e-9;
net.divideParam.trainRatio =  0.8;
net.divideParam.valRatio   =  0.1;
net.divideParam.testRatio  =  0.1;
net                        =  train(net, Xs, Ys);

%Function handle for energy (cell-array input)
f                          =  @(ccell)    f_eval(ccell,net,inputSettings,outputSettings,T);

%Function handle for mu (cell-array input)
mu                         =  @(ccell)   mu_eval(ccell, net, inputSettings, outputSettings, T);

%Fit S=dcdmu by inverting dmudc matrix
S                          =  @(ccell) Sinv_eval(ccell, net, inputSettings, outputSettings, T);

%Convex hull
if convex_yes==1
    %Penalty
    lambda     =  1e9;
    %Hull
    K          =  convhulln(X_Fit,{'QJ'});      % indices of points forming each facet
    hullV      =  unique(K(:));                 % unique indices of vertices on the hull
    V          =  X_Fit(hullV,:);               % convex hull vertices
    
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

%Display fitting R2
disp(['Goodness of fit in energy:',num2str( sqrt(mean((g - f(X(1:end-1))).^2)) )])

end



%% Function helpers below

function [logo] = logo(x)
%Add buffer
x(x<=0) = x(x<=0)+1e-12;
x(x>=1) = x(x>=1)-1e-12;
logo    = log(x);
end


function Y = f_eval(ccell, net,inputSettings,outputSettings, T)
%Flatten cell input into [d × N]
d = numel(ccell);
N = numel(ccell{1});
X = zeros(d, N);
for i = 1:d
    X(i,:)  = ccell{i}(:)';
end
%ScaleX
Xs          =  mapminmax('apply', X, inputSettings);
%Evaluate
Y           =  net(Xs);
Y           =  mapminmax('reverse', Y, outputSettings);
%Add entropy
X           = [X ; 1-sum(X)];
Y           =  Y + 8.314*T(1)*sum(X.*logo(X));
% Reshape back to original input shape
Y = reshape(real(Y), size(ccell{1}));
end



function mu_cell = mu_eval(ccell, net, inputSettings, outputSettings, T)
% Flatten cell input into [d × N]
d = numel(ccell);
N = numel(ccell{1});      % flatten all cells
X = zeros(d, N);
for i = 1:d
    X(i,:) = ccell{i}(:)';
end

% Scale inputs
Xs = mapminmax('apply', X, inputSettings);

% Compute derivatives w.r.t. scaled inputs
% Extract weights from fitnet
W1 = net.IW{1,1}; b1 = net.b{1};
W2 = net.LW{2,1}; b2 = net.b{2};
W3 = net.LW{3,2}; b3 = net.b{3};

% Hidden layer activations (tanh)
Z1 = tanh(W1*Xs + b1);    % hidden1 × N
Z2 = tanh(W2*Z1 + b2);    % hidden2 × N

% Compute Jacobian (df/dX) using chain rule
mu_mat = zeros(d, N);
for i = 1:N
    J1          = diag(1 - Z1(:,i).^2) * W1;   % hidden1 × d
    J2          = diag(1 - Z2(:,i).^2) * W2;   % hidden2 × hidden1
    mu_mat(:,i) = (W3 * J2 * J1)';    % 1 × d, transpose to column
end

% Scale back to original X and Y
% Input scaling factor
Xrange =  inputSettings.xmax - inputSettings.xmin;
Xscale =  1      ./ Xrange(:);
mu_mat =  mu_mat .* Xscale;

% Output scaling factor
Yscale =  outputSettings.xmax - outputSettings.xmin;
mu_mat =  mu_mat * Yscale;

% Add configurational entropy term: d/dc (R*T*c*log(c))
cd     =  1-sum(X);
mu_mat =  mu_mat + 8.314*T(1)*(logo(X)-logo(cd));

% Reshape back into cell array matching input
mu_cell = cell(1,d);
for i = 1:d
    mu_cell{i} = real(reshape(mu_mat(i,:), size(ccell{i})));
end
end



function Sinv_cell = Sinv_eval(ccell, net, inputSettings, outputSettings, T)
% Compute inverse Hessian (Sinv) using finite difference on mu
d      = numel(ccell);
N      = numel(ccell{1});
% Reference mu
mu_ref = mu_eval(ccell, net, inputSettings, outputSettings, T);
% Allocate Hessians: d × d × N
H      = zeros(d, d, N);
% Step size for finite difference
dc     = 1e-6;
% Loop over perturbations in each concentration
for j = 1:d
    ccell_div    = ccell;
    ccell_div{j} = ccell{j} + dc;
    mu_div       = mu_eval(ccell_div, net, inputSettings, outputSettings, T);
    % Build column j of Hessian for all points
    for i = 1:d
        % Flatten difference into [1×N]
        diff_vec = (mu_div{i}(:)' - mu_ref{i}(:)') / dc;
        H(i,j,:) = diff_vec;
    end
end
% Invert Hessian for each point
I        = repmat(eye(d),1,1,N);
Hinv_all = pagemldivide(H, I);
% Reshape back into cell components
Sinv_cell = cell(d,d);
for i = 1:d
    for j = 1:d
        Sinv_cell{i,j} = reshape(Hinv_all(i,j,:), size(ccell{1}));
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
