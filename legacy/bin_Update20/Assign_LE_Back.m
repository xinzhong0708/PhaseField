function [c, mu_e, chi] = Assign_LE_Back(c, mu_e, chi, c_sub, mu_sub, chi_sub, phase_ids, mask)
% ---- assign c back for selected phases ----
if numel(phase_ids) ~= numel(c_sub)
    error('numel(phase_ids) must match numel(c_sub).');
end

for k = 1:numel(phase_ids)
    ip = phase_ids(k);
    for ic = 1:numel(c_sub{k})
        tmp = c{ip}{ic};
        tmp(mask) = c_sub{k}{ic};
        c{ip}{ic} = tmp;
    end
end

% ---- assign mu_e back ----
for ie = 1:numel(mu_sub)
    tmp = mu_e{ie};
    tmp(mask) = mu_sub{ie};
    mu_e{ie} = tmp;
end

% ---- assign chi back ----
[r1, r2] = size(chi_sub);
for i = 1:r1
    for j = 1:r2
        tmp = chi{i,j};
        tmp(mask) = chi_sub{i,j};
        chi{i,j} = tmp;
    end
end

end