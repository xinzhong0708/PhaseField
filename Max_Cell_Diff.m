function dmax = Max_Cell_Diff(A,B)
%MAX_CELL_DIFF maximum absolute difference between two cell arrays of fields.
%
% A and B are usually:
%   STATE.E
%   STATE.mu_e
%   STATE.c{ip}
%
% Each cell contains a ny x nx array.

dmax = 0;

for i = 1:numel(A)
    if isempty(A{i}) || isempty(B{i})
        continue
    end

    d = max(abs(A{i}(:) - B{i}(:)));

    if d > dmax
        dmax = d;
    end
end

end