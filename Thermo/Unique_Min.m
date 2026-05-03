function [M_unique,g_unique,rows_to_keep] = Unique_Min(M,g)
M_round       =  round(M*1e12)/1e12;  

g_flat        =  g(:);
[~,~,ic]      =  unique(M_round,'rows');

% Find minimal g per unique composition
g_min         =  accumarray(ic, g_flat, [], @min);

% Now keep only the rows with minimal g
rows_to_keep  =  false(size(g_flat));
for k = 1:numel(g_min)
    idx = find(ic==k);                % all indices with this composition
    [~,imin] = min(g_flat(idx));      % index of minimal g
    rows_to_keep(idx(imin)) = true;   % keep only that one
end

M_unique = M(rows_to_keep,:);
g_unique = g_flat(rows_to_keep);

end