function STATE = Snap_PureP(STATE,PARAM,MODEL)
%SNAP_PUREP Snap nearly pure p to exactly pure and remove tiny p tails.

p_cut  = 1 - 1e-12;
p_zero = 1e-12;

if isfield(PARAM,'p_snap_cut')
    p_cut = PARAM.p_snap_cut;
end

if isfield(PARAM,'p_zero_cut')
    p_zero = PARAM.p_zero_cut;
end

STATE.p = Calc_p(MODEL,STATE.phi);

[pmax,idmax] = max(STATE.p,[],3);
mask_pure = pmax > p_cut;

% Snap nearly pure grids to exactly one grain
for ig = 1:size(STATE.phi,3)
    A = STATE.phi(:,:,ig);
    A(mask_pure & idmax ~= ig) = 0;
    A(mask_pure & idmax == ig) = 1;
    STATE.phi(:,:,ig) = A;
end

% Recompute p after pure snapping
STATE.p = Calc_p(MODEL,STATE.phi);

% Remove tiny p tails
mask_zero = STATE.p < p_zero;

for ig = 1:size(STATE.phi,3)
    A = STATE.phi(:,:,ig);
    A(mask_zero(:,:,ig)) = 0;
    STATE.phi(:,:,ig) = A;
end

% Renormalize phi where needed
s = sqrt(sum(STATE.phi.^2,3));
mask = s > 1e-30;

for ig = 1:size(STATE.phi,3)
    A = STATE.phi(:,:,ig);
    A(mask) = A(mask)./s(mask);
    STATE.phi(:,:,ig) = A;
end

STATE.p = Calc_p(MODEL,STATE.phi);

end



