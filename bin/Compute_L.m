function [PARAM] = Compute_L(STATE,PARAM)
%COMPUTE_L Compute local Allen-Cahn mobility from current phase compositions.
%
% Formula:
%   L = L_fac * 4*m/(3*kap*zeta)
%
% with
%   zeta = sum_i (eA_i - eB_i)^2 / M_i
%
% STATE.e must already be updated from STATE.c.

[ny,nx,Np] = size(STATE.p);
Ne         = numel(STATE.E);

% ------------------------------------------------------------
% Parameters
% ------------------------------------------------------------
L_fac   = 1.0;
p_cut   = 1e-4;
M_min   = 1e-30;
z_min   = 1e-30;
ngrow   = 2;

if isfield(PARAM,'L_fac'),   L_fac = PARAM.L_fac;   end
if isfield(PARAM,'L_p_cut'), p_cut = PARAM.L_p_cut; end
if isfield(PARAM,'L_ngrow'), ngrow = PARAM.L_ngrow; end

% ------------------------------------------------------------
% Get m and kappa
% ------------------------------------------------------------
if isfield(PARAM,'m')
    m0 = PARAM.m;
else
    tmp = PARAM.Lm./max(PARAM.L,eps);
    tmp = tmp(isfinite(tmp));
    m0  = median(tmp(:));
end

if isfield(PARAM,'kap')
    kap0 = PARAM.kap;
else
    tmp  = PARAM.LK./max(PARAM.L,eps);
    tmp  = tmp(isfinite(tmp));
    kap0 = median(tmp(:));
end

% ------------------------------------------------------------
% Old L as fallback
% ------------------------------------------------------------
if isfield(PARAM,'L') && ~isempty(PARAM.L)
    if isscalar(PARAM.L)
        L_old = PARAM.L*ones(ny,nx);
    else
        L_old = PARAM.L;
    end
else
    L_old = ones(ny,nx);
end

% ------------------------------------------------------------
% Pairwise local L, expanded around interface
% ------------------------------------------------------------
L_sum = zeros(ny,nx);
W_sum = zeros(ny,nx);

PARAM.L_pair = cell(Np,Np);

for ia = 1:Np-1
    for ib = ia+1:Np

        zeta = zeros(ny,nx);

        for ie = 1:Ne

            de = STATE.e{ia}{ie} - STATE.e{ib}{ie};

            Mloc = PARAM.M{ie};
            if isscalar(Mloc)
                Mloc = Mloc*ones(ny,nx);
            end

            zeta = zeta + de.^2 ./ max(abs(Mloc),M_min);

        end

        zeta  = max(zeta,z_min);
        Lpair = L_fac * 4*m0 ./ (3*kap0.*zeta);

        if isfield(PARAM,'L_min')
            Lpair = max(Lpair,PARAM.L_min);
        end
        if isfield(PARAM,'L_max')
            Lpair = min(Lpair,PARAM.L_max);
        end

        PARAM.L_pair{ia,ib} = Lpair;
        PARAM.L_pair{ib,ia} = Lpair;

        % Original two-phase interface mask
        mask = STATE.p(:,:,ia) > p_cut & STATE.p(:,:,ib) > p_cut;

        % Expand interface mask by ngrow grid cells
        if ngrow > 0
            ker       = ones(2*ngrow+1,2*ngrow+1);
            mask_grow = conv2(double(mask),ker,'same') > 0;
        else
            mask_grow = mask;
        end

        % Weight. Use grown mask, but weight strongest near real interface.
        w = STATE.p(:,:,ia).*STATE.p(:,:,ib);
        w(mask_grow & w == 0) = max(w(mask),[],'all');

        if isempty(w(mask)) || max(w(mask),[],'all') == 0
            w(mask_grow) = 1;
        end

        L_sum(mask_grow) = L_sum(mask_grow) + w(mask_grow).*Lpair(mask_grow);
        W_sum(mask_grow) = W_sum(mask_grow) + w(mask_grow);

    end
end

% ------------------------------------------------------------
% Construct final L map
% ------------------------------------------------------------
L_new           = L_old;
mask_L          = W_sum > 0;
L_new(mask_L)   = L_sum(mask_L)./W_sum(mask_L);

% Sharp initial two-phase map may have no overlap.
% Then assign pair value near phase boundary using gradient-like mask.
if ~any(mask_L(:)) && Np == 2 && ~isempty(PARAM.L_pair{1,2})

    p1 = STATE.p(:,:,1);
    p2 = STATE.p(:,:,2);

    mask = false(ny,nx);
    mask(:,1:end-1) = mask(:,1:end-1) | abs(p1(:,2:end)-p1(:,1:end-1)) > 0;
    mask(:,2:end)   = mask(:,2:end)   | abs(p1(:,2:end)-p1(:,1:end-1)) > 0;
    mask(1:end-1,:) = mask(1:end-1,:) | abs(p1(2:end,:)-p1(1:end-1,:)) > 0;
    mask(2:end,:)   = mask(2:end,:)   | abs(p1(2:end,:)-p1(1:end-1,:)) > 0;

    ker       = ones(2*ngrow+1,2*ngrow+1);
    mask_grow = conv2(double(mask),ker,'same') > 0;

    L_new(mask_grow) = PARAM.L_pair{1,2}(mask_grow);

end

% ------------------------------------------------------------
% Save fields used by AC solver
% ------------------------------------------------------------
PARAM.L  = L_new;
PARAM.Lm = L_new*m0;
PARAM.LK = L_new*kap0;
PARAM.LL = PARAM.L;
PARAM.LM = PARAM.Lm;

end