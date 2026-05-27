function [n_new,INFO] = Perturb_n_LE(n,eps_rel,method)
%PERTURB_N_LE Small deterministic perturbation of endmember compositions.
%
% PURPOSE
%   Add trace compositional flexibility to the endmember-composition
%   matrix used by the phase-field LE calculation.  This can remove exact
%   rank degeneracy in the phase compositional Jacobian.
%
% INPUT
%   n        Nem x Ne matrix. Rows are endmembers; the last endmember is
%            the dependent endmember used by PhaseThermo.
%   eps_rel  Relative perturbation amplitude. A good diagnostic sweep is:
%               [0 1e-10 1e-8 1e-6 1e-4]
%   method   'trace' : positive deterministic trace additions (default).
%                      This is more interpretable: it introduces tiny
%                      pseudo-solubilities without making negative n.
%            'svd'   : floor deficient singular values of the LE Jacobian.
%                      This guarantees numerical rank when possible, but
%                      may introduce tiny negative entries and should be
%                      used only as a numerical diagnostic.
%
% OUTPUT
%   n_new    Perturbed endmember matrix.
%   INFO     Diagnostics of the LE-relevant Jacobian:
%               J = n(1:end-1,:) - n(end,:)
%
% IMPORTANT
%   The rank relevant to LE is rank(J), not rank(n).  A phase with Nc
%   independent endmember variables cannot have rank larger than Nc.
%   Therefore full elemental rank is possible only when Nc >= Ne.
%
% EXAMPLE
%   [n,INFO] = Perturb_n_LE(td.n_em(:,1:end-1),1e-8,'trace');
%   fprintf('%s: rank(J) %d -> %d, rcond %.3e -> %.3e\n', ...
%           phs_name{1},INFO.rank0,INFO.rank1,INFO.rcond0,INFO.rcond1)

if nargin < 2 || isempty(eps_rel), eps_rel = 0;       end
if nargin < 3 || isempty(method),  method  = 'trace'; end

if eps_rel < 0
    error('eps_rel must be non-negative.')
end

[Nem,Ne] = size(n);
Nc       = Nem-1;

if Nc < 1
    error('n must contain at least two endmembers.')
end

J0       = n(1:Nc,:) - n(end,:);
scale    = max([1; abs(n(:)); abs(J0(:))]);
eps_abs  = eps_rel*scale;

INFO.method                 = method;
INFO.eps_rel                = eps_rel;
INFO.eps_abs                = eps_abs;
INFO.Nem                    = Nem;
INFO.Nc                     = Nc;
INFO.Ne                     = Ne;
INFO.max_rank_possible      = min(Nc,Ne);
INFO.full_element_possible  = Nc >= Ne;
INFO.rank0                  = rank(J0);
INFO.rcond0                 = Local_Rcond(J0);
INFO.min_n0                 = min(n(:));

n_new = n;

if eps_abs == 0
    INFO.rank1              = INFO.rank0;
    INFO.rcond1             = INFO.rcond0;
    INFO.min_n1             = INFO.min_n0;
    INFO.max_abs_change     = 0;
    INFO.changed            = false;
    return
end

switch lower(method)

    case 'trace'

        %Positive deterministic pattern. Each endmember receives a
        %different trace amount of each element. Consequently the row
        %differences entering J are perturbed rather than cancelled.
        x = (1:Nem)'/(Nem+1);
        P = zeros(Nem,Ne);

        for ie = 1:Ne
            P(:,ie) = x.^ie;
        end

        P     = P/max(P(:));
        n_new = n + eps_abs.*P;

    case 'svd'

        %Minimum-style numerical rank regularization of the LE Jacobian.
        %The dependent endmember is kept unchanged. This can produce tiny
        %negative entries and is intended only to diagnose rank problems.
        [U,S,V] = svd(J0,'full');
        ns      = min(size(S));

        for i = 1:ns
            if S(i,i) < eps_abs
                S(i,i) = eps_abs;
            end
        end

        J1               = U*S*V';
        n_new(1:Nc,:)    = n(end,:) + J1;

    otherwise

        error('Unknown method. Use ''trace'' or ''svd''.')

end

J1                       = n_new(1:Nc,:) - n_new(end,:);
INFO.rank1               = rank(J1);
INFO.rcond1              = Local_Rcond(J1);
INFO.min_n1              = min(n_new(:));
INFO.max_abs_change      = max(abs(n_new(:)-n(:)));
INFO.changed             = INFO.max_abs_change > 0;

if strcmpi(method,'svd') && INFO.min_n1 < -100*eps_abs
    warning(['Perturb_n_LE:svdNegativeEntry: SVD perturbation introduced ', ...
             'negative n entries (minimum %.3e). Use only as a numerical ', ...
             'diagnostic, or use method ''trace''.'],INFO.min_n1)
end

end


function r = Local_Rcond(J)

s = svd(J,'econ');

if isempty(s) || max(s) == 0
    r = 0;
else
    r = min(s)/max(s);
end

end
