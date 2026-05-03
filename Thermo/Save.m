function [f,mu,S] = Fit_Energy(convex_yes,g,T,X)

%Order of energy function
ORD           = 4;

%Subtract configurational entropy (will need to be added back)
X_Fit         = [];
g0            =  0;
for i = 1:length(X)
    g0        =  g0 + 8.314*T(1).*(X{i}.*log(X{i}));
    X_Fit     = [X_Fit,X{i}(:)];
end

%Fit the function
id            =  find( ((imag(g)<1e-10)==1).*(g<1e10)  );
X_Fit         =  X_Fit(id,1:end-1);
Y_Fit         =  g(id) - g0(id);
id            =  isfinite(Y_Fit);
pp            =  polyfitn(X_Fit(id,:),real(Y_Fit(id)),ORD);
pv            =  pp.Coefficients;
pw            =  pp.ModelTerms;
c             =  sym('c',[1,length(X)]);
c(end)        =  1-sum(c(1:end-1));
f             =  0;

%Add back entropy
for i = 1:length(X)
    f         =  f + 8.314*T(1).*c(i).*log(c(i));
end

%Add polynomial
c             =  c(1:end-1);
for ip  = 1:length(pv)
    f         =  f + pv(ip).*prod(c.^pw(ip,:));
end

%Fit S=dcdmu by inverting dmudc matrix
Nc            =  length(c);
H             =  hessian(f,c);
Hfun          =  matlabFunction(H,'Vars',num2cell(c));
S             =  @(C) localHinv_component(Hfun,C);  % S is the inverted Hessian

%Chemical potential
mu   = cell(1,Nc);
for i = 1:length(X)-1
    mu_fun    =  matlabFunction(diff(f,c(i)));
    mu{i}     =  @(c)  mu_fun(c{:});
end

%Energy function
f             =  matlabFunction(f);
f             =  @(c)   f(c{:});

%Convex hull
if convex_yes == 1
    K         =  convhulln(X_Fit);
    [A,b]     =  facet_planes(X_Fit,K);

    lambda    =  100;
    eta       =  100;

    [f,mu,S]  = penalized_free_energy({rand(1,10), rand(1,10), rand(1,10)},f, A, b, lambda, eta);

end


%Display fitting R2
disp(['Goodness of fit:',num2str(pp.R2)])

end






% --------------------------
function Hij = localHinv_component(Hfun,C)
% C is a cell array of matrices
if ~iscell(C)
    error('Input must be a cell array of matrices.');
end

sz   = size(C{1});
nPts = numel(C{1});
Nc   = numel(C);

% Flatten inputs
args_flat = cellfun(@(x) x(:), C, 'UniformOutput', false);

% Evaluate Hessians
Hvals = zeros(Nc,Nc,nPts);
for k = 1:nPts
    args_k       = cellfun(@(x) x(k), args_flat, 'UniformOutput', false);
    Hvals(:,:,k) = feval(Hfun, args_k{:});
end

% Invert Hessians
I        = repmat(eye(Nc),1,1,nPts);
Hinv_all = pagemldivide(Hvals,I);

% Return requested component and reshape
for i = 1:Nc
    for j = 1:Nc
        Hij{i,j} = reshape(Hinv_all(i,j,:), sz);
    end
end
end



function [dist, y_proj, alpha] = project_to_hull(V, c)
% PROJECT_TO_HULL  Project composition c onto convex hull of vertices V
%
% V : m x d matrix (rows are hull vertices in d dimensions)
% c : cell array {c1, c2, ..., cd} or numeric 1 x d vector
%
% Returns:
%   dist   : Euclidean distance ||c - y_proj||
%   y_proj : 1 x d projected point on hull
%   alpha  : m x 1 convex weights (sum(alpha)=1, alpha>=0)

[m,d] = size(V);

% --- Convert c to numeric vector ---
if iscell(c)
    x = cell2mat(c(:));      % flatten cell contents into numeric vector
else
    x = c(:);
end

if numel(x) ~= d
    error('Point must have same dimension as V columns');
end

% Quadratic program:
H = 2*(V*V.');     
f = -2*(V*x);      

Aeq = ones(1,m);
beq = 1;
lb  = zeros(m,1);
ub  = [];           

opts = optimoptions('quadprog','Display','off', ...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-12);

[alpha, ~, exitflag] = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], opts);

if exitflag <= 0
    warning('quadprog did not converge (exitflag=%d). Using fallback.', exitflag);
    alpha = V \ x;       
    alpha(alpha < 0) = 0;
    if sum(alpha) == 0
        alpha = ones(m,1)/m;
    else
        alpha = alpha / sum(alpha);
    end
end

y_proj = (alpha.' * V);      
dist   = norm(x(:) - y_proj(:));

end




function [fval, muval, Hval] = penalized_free_energy(c, f_poly, A, b, lambda, eta)
% c: {c1, c2, c3} each N1xN2... array
% f_poly: function handle f_poly(c1,c2,c3)
% A: Fxd matrix of normalized facet normals (rows = a_j')
% b: Fx1 vector of facet offsets
% lambda: penalty weight
% eta: sharpness parameter (e.g. 100)

c1  = c{1}; c2 = c{2}; c3 = c{3};
sz  = size(c1);
X   = [c1(:), c2(:), c3(:)];   % Qxd

% --- Base free energy ---
f_base  = f_poly(c);
mu_base = zeros([sz,3]);
H_base  = zeros([sz,3,3]);

% For base gradient/Hessian you can get them from symbolic matlabFunction
% For now assume mu_base, H_base are available
% (I can show how to generate them automatically if needed)

% --- Penalty computations ---
S    = X*A.' - b.';       % QxF signed facet distances
E    = exp(eta*S);           
den  = 1 + sum(E,2);      % Qx1
dval = (1/eta)*log(den);  % Qx1 smooth distance
Pval = lambda*dval.^2;    % penalty value

% Probabilities
Pj = E ./ den;            % QxF
grad_d = Pj*A;            % Qxd
% Hessian of d: η * sum_j p_j (a_j - grad_d)(a_j - grad_d)^T
Q = size(X,1);
H_d = zeros(Q,3,3);
for i=1:Q
    gi = grad_d(i,:); 
    Hi = zeros(3,3);
    for j=1:size(A,1)
        aj = A(j,:);
        diff = (aj - gi);
        Hi = Hi + Pj(i,j)*(diff.'*diff);
    end
    H_d(i,:,:) = eta*Hi;
end

% Gradient and Hessian of penalty
gradP = zeros(Q,3);
Hpen  = zeros(Q,3,3);
for i=1:Q
    gi = grad_d(i,:);
    di = dval(i);
    Hi = squeeze(H_d(i,:,:));
    gradP(i,:)   = 2*lambda*di*gi;
    Hpen(i,:,:)  = 2*lambda*(gi.'*gi + di*Hi);
end

% --- Combine with base ---
fval  = f_base  + reshape(Pval, sz);
muval = mu_base + reshape(gradP,[sz,3]);
Hval  = H_base  + reshape(Hpen,[sz,3,3]);

end


function [A,b] = facet_planes(M,K)
F = size(K,1); d = size(M,2);
A = zeros(F,d); b = zeros(F,1);
for j=1:F
    verts = M(K(j,:),:);
    N = null(verts(2:end,:) - verts(1,:));
    n = N(:,1); n = n/norm(n);
    A(j,:) = n';
    b(j)   = n'*verts(1,:)';
end
end
