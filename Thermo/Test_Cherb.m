clear;clf;hold on

x1         = linspace(-1,1,50);
[c1,c2,c3] = ndgrid(x1,x1,x1);
X_Fit      = [c1(:) c2(:) c3(:)];
g_fun      =  @(c1,c2,c3) c1.^2 + c2.^3 + 3*c3 + sin(c1) + tan(c2) + 1 + exp(c1+c2*2) - exp(c2.^2+c3.^3) - 10*log(2+c1.^2+c3.^3);
g          =  g_fun(c1,c2,c3);


order      =  5;

%Polynomial fitting
pp         =  polyfitn(X_Fit,g(:),order);

%Cherbishev fitting
[coeffs,fCheb] = chebyshev_fitND(X_Fit, g(:), order);



%Evaluate at new points
Xnew       = rand(1000,3);

%Correct solution
g_correct  = g_fun(Xnew(:,1),Xnew(:,2),Xnew(:,3));

%Polyeval
g_poly     = polyvaln(pp,Xnew);

%Cherb
gnew       = fCheb(Xnew);


plot(g_correct,'ko-')
plot(gnew,'b.','markersize',20)
plot(g_poly,'r.','markersize',20)

sum((g_correct-g_poly).^2)
sum((g_correct-gnew).^2)




function [coeffs, basisFun, scaleX] = chebyshev_fitND(X_Fit, g, order)
% N-D Chebyshev tensor-product fit
% Inputs:
%   X_Fit : N x d matrix of input variables
%   g     : N x 1 target values
%   order : maximum Chebyshev order
% Outputs:
%   coeffs  : coefficients of tensor-product Chebyshev basis
%   basisFun: function handle @(X) evaluating fitted function at points X
%   scaleX  : struct with min/max for scaling

[N,d] = size(X_Fit);

% Scale each variable to [-1,1]
scaleX.min = min(X_Fit,[],1);
scaleX.max = max(X_Fit,[],1);
X_scaled = 2*(X_Fit - scaleX.min) ./ (scaleX.max - scaleX.min) - 1;

% Precompute Chebyshev polynomials for each dimension
Tx = cell(d,1);
for dim = 1:d
    Tx{dim} = zeros(N, order+1);
    x = X_scaled(:,dim);
    Tx{dim}(:,1) = 1;           % T0(x) = 1
    if order >= 1
        Tx{dim}(:,2) = x;       % T1(x) = x
    end
    for k = 2:order
        Tx{dim}(:,k+1) = 2*x.*Tx{dim}(:,k) - Tx{dim}(:,k-1); % recurrence
    end
end

% Construct tensor-product design matrix
nTerms = (order+1)^d;
A = ones(N, nTerms);
% Generate all multi-indices for tensor-product
idx = cell(1,d);
for dim = 1:d
    idx{dim} = 0:order;
end
[grids{1:d}] = ndgrid(idx{:});
multiIdx = zeros(nTerms,d);
for dim = 1:d
    multiIdx(:,dim) = grids{dim}(:);
end

% Fill design matrix
for col = 1:nTerms
    prodVec = ones(N,1);
    for dim = 1:d
        prodVec = prodVec .* Tx{dim}(:, multiIdx(col,dim)+1);
    end
    A(:,col) = prodVec;
end

% Solve least-squares
coeffs = A \ g;

% Function handle for evaluation
basisFun = @(Xnew) eval_chebND(Xnew, coeffs, order, scaleX);

end

%% Evaluation function
function gfit = eval_chebND(Xnew, coeffs, order, scaleX)
% Evaluate N-D Chebyshev tensor-product function
[N,d] = size(Xnew);

% Scale to [-1,1]
X_scaled = 2*(Xnew - scaleX.min) ./ (scaleX.max - scaleX.min) - 1;

% Precompute Chebyshev polynomials
Tx = cell(d,1);
for dim = 1:d
    Tx{dim} = zeros(N, order+1);
    x = X_scaled(:,dim);
    Tx{dim}(:,1) = 1;
    if order >= 1
        Tx{dim}(:,2) = x;
    end
    for k = 2:order
        Tx{dim}(:,k+1) = 2*x.*Tx{dim}(:,k) - Tx{dim}(:,k-1);
    end
end

% Construct tensor-product design matrix
nTerms = (order+1)^d;
A = ones(N, nTerms);
idx = cell(1,d);
for dim = 1:d
    idx{dim} = 0:order;
end
[grids{1:d}] = ndgrid(idx{:});
multiIdx = zeros(nTerms,d);
for dim = 1:d
    multiIdx(:,dim) = grids{dim}(:);
end

for col = 1:nTerms
    prodVec = ones(N,1);
    for dim = 1:d
        prodVec = prodVec .* Tx{dim}(:, multiIdx(col,dim)+1);
    end
    A(:,col) = prodVec;
end

% Evaluate
gfit = A * coeffs;

end
