function Hij = localHinv_component(Hfun, ccell)

d  = numel(ccell);         % number of components
sz = size(ccell{1});       % grid shape
N  = numel(ccell{1});      % number of points

%Flatten input matrices to N×d ---
X = zeros(N, d);
for i = 1:d
    X(:,i) = ccell{i}(:);
end

%Call Hfun with variable number of args ---
args  = num2cell(X.',2);
H_val = Hfun(args{:});

%SAFE reshape to d×d×N ---
H_val = reshape(H_val, [], 1);    % ensure column vector
H_val = reshape(H_val, d^2, N);   % d^2 × N
H_val = reshape(H_val, d, d, N);  % d × d × N

%Invert each d×d page ---
Ipages    = repmat(eye(d), 1, 1, N);
Hinv_all  = pagemldivide(H_val, Ipages);   % d × d × N

%Convert back to d×d cell, each same shape as input ---
Hij = cell(d,d);
for i = 1:d
    for j = 1:d
        Hij{i,j} = reshape(Hinv_all(i,j,:), sz);
    end
end
end