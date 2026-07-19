function STATE = Horizon_Ave(STATE,MODEL)
%HORIZON_AVE Make STATE horizontally averaged but keep original ny,nx size.
%
% It replaces every row by the y-average:
%   A(ny,nx,...) -> repmat(mean(A,1),ny,1,...)
%
% This keeps the model 2D in array size, but makes it pseudo-1D.

if nargin < 2
    MODEL = [];
end

% Get original ny
if isfield(STATE,'phi')
    ny = size(STATE.phi,1);
else
    ny = [];
end

fields = fieldnames(STATE);

for i = 1:numel(fields)

    name = fields{i};

    if strcmp(name,'LE_state')
        STATE.LE_state = [];
        continue
    end

    STATE.(name) = AveY_KeepSize_Local(STATE.(name),ny);

end

% Normalize phi again
if isfield(STATE,'phi') && ~isempty(STATE.phi)
    S = sum(STATE.phi,3);
    STATE.phi = STATE.phi ./ max(S,eps);
end

% Recompute p from averaged phi
if ~isempty(MODEL) && isfield(STATE,'phi')
    STATE.p = Calc_p(MODEL,STATE.phi);
end

% Recompute e from averaged c
if ~isempty(MODEL) && isfield(STATE,'c')
    STATE.e = Calc_e(MODEL.pars,STATE.c);
end

% Active-set memory is no longer valid after averaging
STATE.LE_state = [];

end


function A = AveY_KeepSize_Local(A,ny)

if isnumeric(A) || islogical(A)

    if isempty(A)
        return
    end

    if isempty(ny)
        ny = size(A,1);
    end

    if ndims(A) >= 2 && size(A,1) == ny && ny > 1

        Amean = mean(double(A),1);

        reps = ones(1,ndims(A));
        reps(1) = ny;

        A = repmat(Amean,reps);

    end

elseif iscell(A)

    for i = 1:numel(A)
        A{i} = AveY_KeepSize_Local(A{i},ny);
    end

elseif isstruct(A)

    fn = fieldnames(A);

    for i = 1:numel(fn)
        A.(fn{i}) = AveY_KeepSize_Local(A.(fn{i}),ny);
    end

end

end