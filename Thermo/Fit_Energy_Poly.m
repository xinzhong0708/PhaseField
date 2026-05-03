function [f,mu,S,H] = Fit_Energy_Poly(g,T,X)
%Here polynomial fitting is used.
%Order of energy function
ORD           =  7;

%Subtract configurational entropy (will need to be added back)
X_Fit         = [];
g0            =  0;
for i = 1:length(X)
    g0        =  g0 + 8.314*T(1).*(X{i}.*logo(X{i}));
    X_Fit     = [X_Fit,X{i}(:)];
end

%Fit the function
id            =  find( ((imag(g)<1e-10)==1).*(g<1e10)  );
X_Fit         =  X_Fit(id,1:end-1);
Y_Fit         =  g(id) - g0(id);
id            =  isfinite(Y_Fit);
X_Fit         =  X_Fit(id,:);
Y_Fit         =  real(Y_Fit(id));
pp            =  polyfitn(X_Fit,Y_Fit,ORD);
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
H             =  simplify(hessian(f,c));
Hflat         =  H(:);

%Create function handle for dcdmu
Hfun          =  matlabFunction(Hflat, 'Vars', num2cell(c));
S             =  @(c) localHinv_component(Hfun,c);

%Create function handle for dcdmu
Hfun          =  matlabFunction(Hflat, 'Vars', num2cell(c));
H             =  @(c) localH_component(Hfun,c);

%Chemical potential
mu            =  cell(1,Nc);
for i = 1:Nc
    mu_fun    =  matlabFunction(diff(f,c(i)));
    mu{i}     =  @(c)  mu_fun(c{:});
end
mu            =  @(c) cellfun(@(f) f(c), mu, 'UniformOutput', false);

%Energy function
f             =  matlabFunction(f);
f             =  @(c)   f(c{:});


%Display fitting R2
disp(['Goodness of fit:',num2str(pp.R2)])

end



function [logo] = logo(x)
%Add buffer
x(x<=0) = x(x<=0)+1e-12;
x(x>=1) = x(x>=1)-1e-12;
logo    = log(x);
end


function Hij = localHinv_component(Hfun, ccell, reg_eps)
    if nargin < 3 || isempty(reg_eps)
        reg_eps = 0;
    end

    d  = numel(ccell);         % number of components
    sz = size(ccell{1});       % grid shape
    N  = numel(ccell{1});      % number of points

    % --- Flatten input matrices to N×d ---
    X = zeros(N, d);
    for i = 1:d
        X(:,i) = ccell{i}(:);
    end

    % --- Call Hfun with variable number of args ---
    args  = num2cell(X.',2);
    H_val = Hfun(args{:});

    % --- SAFE reshape to d×d×N ---
    H_val = reshape(H_val, [], 1);    % ensure column vector
    H_val = reshape(H_val, d^2, N);   % d^2 × N
    H_val = reshape(H_val, d, d, N);  % d × d × N

    % --- Optional small regularization for numerical stability ---
    if reg_eps ~= 0
        H_val = H_val + reg_eps * repmat(eye(d), 1, 1, N);
    end

    % --- Invert each d×d page ---
    Ipages    = repmat(eye(d), 1, 1, N);
    Hinv_all  = pagemldivide(H_val, Ipages);   % d × d × N

    % --- Convert back to d×d cell, each same shape as input ---
    Hij = cell(d,d);
    for i = 1:d
        for j = 1:d
            Hij{i,j} = reshape(Hinv_all(i,j,:), sz);
        end
    end
end



function H_val = localH_component(Hfun, ccell, reg_eps)
    if nargin < 3 || isempty(reg_eps)
        reg_eps = 0;
    end

    d  = numel(ccell);         % number of components
    sz = size(ccell{1});       % grid shape
    N  = numel(ccell{1});      % number of points

    % --- Flatten input matrices to N×d ---
    X = zeros(N, d);
    for i = 1:d
        X(:,i) = ccell{i}(:);
    end

    % --- Call Hfun with variable number of args ---
    args  = num2cell(X.',2);
    H_val = Hfun(args{:});

    % --- SAFE reshape to d×d×N ---
    H_val = reshape(H_val, [], 1);    % ensure column vector
    H_val = reshape(H_val, d^2, N);   % d^2 × N
    H_val = reshape(H_val, d, d, N);  % d × d × N

    % --- Optional small regularization for numerical stability ---
    if reg_eps ~= 0
        H_val = H_val + reg_eps * repmat(eye(d), 1, 1, N);
    end

end
