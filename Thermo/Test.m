clear

% Example: 2D demo
c1 = linspace(0,1,100);
c2 = linspace(0,2,200);
[C1,C2] = ndgrid(c1,c2);
fData = sin(C1) + cos(C2);

fh = makeLookupND({c1,c2}, fData);

surf(C1,C2,fData);shading interp

% Evaluate like a function handle
X_test = rand(3,3);
Y_test = rand(3,3);

val1 = fh({X_test,Y_test});
val2 = griddata(C1,C2,fData,X_test,Y_test);

val1-val2

% save('temp','fh')





function fh = makeLookupND(cGrid, fData)
% makeLookupND  Create a piecewise linear lookup function handle in N-D
%
%   fh = makeLookupND(cGrid, fData)
%
%   INPUTS:
%     cGrid : cell array {1..Nc}, each containing vector of grid points
%     fData : N-D array of function values at the grid points
%
%   OUTPUT:
%     fh    : function handle @(Ccell) returning interpolated values
%             where Ccell = {C1,C2,...}, each Ci is matrix (same shape)

Nc = numel(cGrid);
gridSize = cellfun(@numel, cGrid);

% Precompute bases and slopes
bases = zeros(gridSize-1);
slopes = cell(1,Nc);
for d = 1:Nc
    slopes{d} = zeros(gridSize-1);
end

allSubs = arrayfun(@(n) 1:n-1, gridSize, 'UniformOutput', false);
[idxList{1:Nc}] = ndgrid(allSubs{:});

for k = 1:numel(idxList{1})
    I = cellfun(@(X) X(k), idxList, 'UniformOutput', false);
    cornerVals = arrayfun(@(d) cGrid{d}(I{d}), 1:Nc);
    f0 = fData(I{:});
    bases(I{:}) = f0;

    for d = 1:Nc
        J = I;
        J{d} = J{d}+1;
        fd = fData(J{:});
        dx = cGrid{d}(J{d}) - cGrid{d}(I{d});
        slopes{d}(I{:}) = (fd - f0)/dx;
    end
end

% Return function handle
fh = @(Ccell) evalLookup(Ccell, cGrid, bases, slopes);
end

% ================== Internal evaluator ==================
function vals = evalLookup(Ccell, cGrid, bases, slopes)
Nc = numel(cGrid);
sz = size(Ccell{1});
vals = zeros(sz);

% Flatten everything for looping
nPts = numel(Ccell{1});
Cmat = zeros(nPts,Nc);
for d = 1:Nc
    Cmat(:,d) = Ccell{d}(:);
end

outVals = zeros(nPts,1);
for p = 1:nPts
    c_query = Cmat(p,:);

    % find bin in each dim
    idx = zeros(1,Nc);
    cornerVals = zeros(1,Nc);
    for d = 1:Nc
        [~,~,bin] = histcounts(c_query(d), cGrid{d});
        if bin==0
            error('Query out of bounds in dimension %d',d);
        end
        idx(d) = bin;
        cornerVals(d) = cGrid{d}(bin);
    end

    subs = num2cell(idx);
    f0 = bases(subs{:});
    val = f0;
    for d = 1:Nc
        slope_d = slopes{d}(subs{:});
        val = val + slope_d*(c_query(d) - cornerVals(d));
    end
    outVals(p) = val;
end

vals(:) = outVals;
end
