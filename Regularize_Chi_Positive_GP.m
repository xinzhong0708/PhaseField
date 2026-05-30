function [STATE,DIAG] = Regularize_Chi_Positive_GP(STATE,PARAM)
%REGULARIZE_CHI_POSITIVE_GP Convexify STATE.chi for mu-primary/GP-like runs.
%
% Use this after LE_Run_FromMu when you do NOT want spinodal decomposition.
% It clips eigenvalues of chi=dE/dmu so the CH solve sees a positive
% response everywhere.

Ne = size(STATE.chi,1);
sz = size(STATE.chi{1,1});
N  = numel(STATE.chi{1,1});

floor_abs = GetParamLocal(PARAM,'GP_chi_floor_abs',1e-12);
floor_rel = GetParamLocal(PARAM,'GP_chi_floor_rel',1e-8);
chi_cap   = GetParamLocal(PARAM,'GP_chi_cap',inf);

Hpg = zeros(Ne,Ne,N);
for i = 1:Ne
    for j = 1:Ne
        Hpg(i,j,:) = reshape(STATE.chi{i,j},1,1,N);
    end
end

min_before = nan(1,N);
min_after  = nan(1,N);
n_clipped  = 0;
n_capped   = 0;
has_bad    = false;

for inode = 1:N
    H = 0.5*(Hpg(:,:,inode) + Hpg(:,:,inode).');

    if any(~isfinite(H(:)))
        has_bad = true;
        H(~isfinite(H)) = 0;
    end

    [V,D] = eig(H);
    lam   = diag(D);
    min_before(inode) = min(lam);

    scale = max(1,max(abs(lam)));
    lam_floor = max(floor_abs,floor_rel*scale);

    bad = lam < lam_floor;
    if any(bad)
        lam(bad) = lam_floor;
        n_clipped = n_clipped + nnz(bad);
    end

    if isfinite(chi_cap)
        too_large = lam > chi_cap;
        if any(too_large)
            lam(too_large) = chi_cap;
            n_capped = n_capped + nnz(too_large);
        end
    end

    Hnew = V*diag(lam)*V.';
    Hnew = 0.5*(Hnew+Hnew.');
    Hpg(:,:,inode) = Hnew;
    min_after(inode) = min(eig(Hnew));
end

for i = 1:Ne
    for j = 1:Ne
        STATE.chi{i,j} = reshape(Hpg(i,j,:),sz);
    end
end

DIAG.min_before        = reshape(min_before,sz);
DIAG.min_after         = reshape(min_after,sz);
DIAG.global_min_before = min(min_before,[],'omitnan');
DIAG.global_min_after  = min(min_after,[],'omitnan');
DIAG.n_negative_nodes  = nnz(min_before < 0);
DIAG.n_clipped_eigs    = n_clipped;
DIAG.n_capped_eigs     = n_capped;
DIAG.has_bad_input     = has_bad;

end

function v = GetParamLocal(S,name,default)
if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
    v = S.(name);
else
    v = default;
end
end
