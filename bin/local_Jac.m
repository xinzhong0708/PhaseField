function [Jpages] = local_Jac(JacFun,ccell)
d        = numel(ccell);
N        = numel(ccell{1});
Xcols    = cellfun(@(M) M(:).', ccell, 'uni', 0);
Jflat    = JacFun(Xcols{:});                 % (r*d)×N
[Mrow,~] = size(Jflat);
r        = Mrow / d;                         % number of independent element
Jpages   = reshape(Jflat, r, d, N);
end