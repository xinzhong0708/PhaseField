function E = EnforceMeanE_Local(E,E_old)
Ne = numel(E);
for ie = 1:Ne
    target_mean = mean(E_old{ie}(:));
    new_mean    = mean(E{ie}(:));
    E{ie}       = E{ie} + target_mean - new_mean;
end
end