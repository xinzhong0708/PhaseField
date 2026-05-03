function H_val = localH_component(Hfun, ccell)
% Hfun returns flattened d^2 x 1 Hessian at one point
% ccell is {c1_grid, c2_grid, ..., cd_grid}
% output H_val is d x d x N

d = numel(ccell);
N = numel(ccell{1});

X = zeros(N,d);
for i = 1:d
    X(:,i) = ccell{i}(:);
end

args  = num2cell(X',2);   % one point only
h     = Hfun(args{:});
h     = double(h(:));
H_val = reshape(h,d,d,[]);

end