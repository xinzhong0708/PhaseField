function e = PhaseE_Local(pars,c)

if iscell(c)
    c = cell2mat(c');
end

n_real = pars.n;

if isfield(pars,'nN') && ~isempty(pars.nN)
    n_pen = pars.nN;
else
    n_pen = zeros(0,size(n_real,2));
end

n_all = [n_real; n_pen];

nAll = size(n_all,1);

if size(c,1) == nAll
    c_full = c;
elseif size(c,1) == nAll-1
    c_full = [c; 1-sum(c,1)];
else
    error('PhaseE_Local: input c must have %d or %d rows, got %d.', ...
          nAll-1,nAll,size(c,1));
end

e = (c_full.' * n_all);
e = (e ./ sum(e,2)).';
e = e(1:end-1,:);

end