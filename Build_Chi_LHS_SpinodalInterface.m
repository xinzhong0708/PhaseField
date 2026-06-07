function chi_lhs = Build_Chi_LHS_SpinodalInterface(chi_phys,PARAM,NUM)
%BUILD_CHI_LHS_SPINODALINTERFACE Solver-side chi stabilization.
%
% This function does not change thermodynamics.
% It only builds a temporary chi_lhs for ACCH/corrector matrices.
%
% Rule:
%   outside PARAM.chi_lhs_mask:
%       chi_lhs = chi_phys
%
%   inside PARAM.chi_lhs_mask:
%       if min(eig(chi_phys)) < chi_floor, project to positive definite
%       otherwise keep chi_phys unchanged
%
% Therefore good phases with positive chi are never changed, even if they
% accidentally fall inside the mask.

Ne = size(chi_phys,1);
[ny,nx] = size(chi_phys{1,1});
N = ny*nx;

chi_floor = 1e-10;

if nargin >= 3 && isfield(NUM,'chi_lhs_floor')
    chi_floor = NUM.chi_lhs_floor;
end

mask = zeros(1,N);

if isfield(PARAM,'chi_lhs_mask') && ~isempty(PARAM.chi_lhs_mask)
    mask = reshape(PARAM.chi_lhs_mask,1,N);
end

chi_lhs = chi_phys;

for n = 1:N

    if mask(n) <= 0
        continue
    end

    A = zeros(Ne,Ne);

    for i = 1:Ne
        for j = 1:Ne
            tmp = chi_phys{i,j};
            A(i,j) = tmp(n);
        end
    end

    A = 0.5*(A + A.');

    if any(~isfinite(A(:)))
        A2 = chi_floor*eye(Ne);

    else
        [V,D] = eig(A);
        lam = diag(D);

        %Only change truly spinodal/negative nodes.
        if min(lam) >= chi_floor
            continue
        end

        lam = max(lam,chi_floor);

        A2 = V*diag(lam)*V.';
        A2 = 0.5*(A2 + A2.');
    end

    for i = 1:Ne
        for j = 1:Ne
            tmp = chi_lhs{i,j};
            tmp(n) = A2(i,j);
            chi_lhs{i,j} = tmp;
        end
    end
end

end
