function [STATE,PARAM] = Calc_Kappa_Chi_ACCH_Ramp(STATE,PARAM,MODEL,PHYS,NUM)
%CALC_KAPPA_CHI_ACCH_RAMP Build kappa ramp with nonzero interface floor.
%
% This version does NOT modify chi.
%
% Purpose:
%   PARAM.kappa_eff is reduced near phase interfaces, but it is not forced
%   all the way to zero inside a phase. This avoids the dangerous case:
%
%       chi < 0 but kappa_eff = 0
%
%   which corresponds to unregularized spinodal diffusion.
%
% Rule:
%   outside/absent phase        : kappa = 0
%   interface part of phase     : kappa = kappa_min_frac * phase_kappa
%   transition to core          : smooth ramp
%   phase core                  : kappa = phase_kappa
%
% Typical settings:
%
%   PARAM.kappa_phase       = PHYS.kappa .* cellfun(@(x) size(x.n,1) > 1, MODEL.pars);
%   PARAM.kappa_p_interface = 0.999;
%   PARAM.kappa_p_presence  = 1e-5;
%   PARAM.kappa_zero_width  = 1;
%   PARAM.kappa_ramp_width  = 4;
%   PARAM.kappa_min_frac    = 0.02;
%
%   [STATE,PARAM] = Calc_Kappa_Chi_ACCH_Ramp(STATE,PARAM,MODEL,PHYS,NUM);
%
% Notes:
%   STATE.chi is left untouched.
%   Old fields STATE.chi_acch and STATE.chi_phys are removed if present.

[ny,nx,Ng] = size(STATE.p);

% Remove old solver-chi fields to avoid accidental reuse
if isfield(STATE,'chi_acch')
    STATE = rmfield(STATE,'chi_acch');
end
if isfield(STATE,'chi_phys')
    STATE = rmfield(STATE,'chi_phys');
end

% Defaults
p_interface    = 0.999;
p_presence     = 1e-5;
zero_width     = 2;
ramp_width     = 4;
kappa_min_frac = 0.1;

if isfield(PARAM,'kappa_p_interface')
    p_interface = PARAM.kappa_p_interface;
end
if isfield(PARAM,'kappa_p_presence')
    p_presence = PARAM.kappa_p_presence;
end
if isfield(PARAM,'kappa_zero_width')
    zero_width = PARAM.kappa_zero_width;
end
if isfield(PARAM,'kappa_ramp_width')
    ramp_width = PARAM.kappa_ramp_width;
end
if isfield(PARAM,'kappa_min_frac')
    kappa_min_frac = PARAM.kappa_min_frac;
end

zero_width     = max(0,round(zero_width));
ramp_width     = max(1,round(ramp_width));
kappa_min_frac = min(max(kappa_min_frac,0),1);

% Collapse grains to thermodynamic phases
phase_index = MODEL.phase_index(:).';
phase_id    = unique(phase_index,'stable');
Nphase      = numel(phase_id);

p_phase = zeros(ny,nx,Nphase);
grain_to_phase = zeros(1,Ng);

for iph = 1:Nphase
    grains = find(phase_index == phase_id(iph));
    grain_to_phase(grains) = iph;
    p_phase(:,:,iph) = sum(STATE.p(:,:,grains),3);
end

% Determine phase kappa after collapse
kappa_phase = zeros(1,Nphase);

if isfield(PARAM,'kappa_phase') && ~isempty(PARAM.kappa_phase)

    kp_in = PARAM.kappa_phase(:).';

    if numel(kp_in) == Ng
        for iph = 1:Nphase
            grains = find(grain_to_phase == iph);
            kappa_phase(iph) = max(kp_in(grains));
        end

    elseif numel(kp_in) == Nphase
        kappa_phase = kp_in;

    elseif isscalar(kp_in)
        kappa_phase(:) = kp_in;
    end

elseif isfield(PARAM,'kappa_phase_name') && isfield(MODEL,'phase_name')

    for iph = 1:Nphase
        grains = find(grain_to_phase == iph);
        ig0 = grains(1);

        if strcmpi(MODEL.phase_name{ig0},PARAM.kappa_phase_name)
            kappa_phase(iph) = PHYS.kappa;
        end
    end

else

    kappa_phase(:) = PHYS.kappa;

end

% Build weights and kappa_eff
weight_phase     = zeros(ny,nx,Nphase);
weight_phase_raw = zeros(ny,nx,Nphase);
kappa_eff        = zeros(ny,nx);

for iph = 1:Nphase

    pcur = p_phase(:,:,iph);
    kp   = kappa_phase(iph);

    if kp == 0
        continue
    end

    present = pcur > p_presence;

    % Cells below p_interface are treated as phase boundary/outside.
    zero_mask = pcur < p_interface;

    if ~any(zero_mask(:))

        % This phase fills the whole model.
        w_raw = ones(ny,nx);

    else

        max_dist = zero_width + ramp_width + 1;

        dist  = inf(ny,nx);
        known = zero_mask;
        dist(known) = 0;

        for ir = 1:max_dist
            known_new = DilateNoWrap(known);
            add = known_new & ~known;
            dist(add) = ir;
            known = known_new;
        end

        dist(~isfinite(dist)) = max_dist;

        x = (dist - zero_width) ./ ramp_width;
        x = min(max(x,0),1);

        % Smoothstep 0 -> 1
        w_raw = x.^2 .* (3 - 2*x);
    end

    % Outside/absent phase should have no kappa from this phase.
    w_raw(~present) = 0;

    % Nonzero interface floor inside the phase.
    % At interface: w_raw = 0 -> w = kappa_min_frac
    % At core:      w_raw = 1 -> w = 1
    w = w_raw;
    w(present) = kappa_min_frac + (1-kappa_min_frac).*w_raw(present);

    weight_phase_raw(:,:,iph) = w_raw;
    weight_phase(:,:,iph)     = w;

    kappa_eff = max(kappa_eff,kp.*w);

end

% Scalar weight for diagnostics
active_kappa_phase = kappa_phase ~= 0;

if any(active_kappa_phase)
    w_core = max(weight_phase(:,:,active_kappa_phase),[],3);
    w_raw_core = max(weight_phase_raw(:,:,active_kappa_phase),[],3);
else
    w_core = max(weight_phase,[],3);
    w_raw_core = max(weight_phase_raw,[],3);
end

% Store fields
PARAM.kappa_eff             = kappa_eff;
PARAM.kappa_weight_phase    = weight_phase;
PARAM.kappa_weight_phase_raw= weight_phase_raw;
PARAM.kappa_chi_weight      = w_core;
PARAM.kappa_raw_weight      = w_raw_core;
PARAM.kappa_phase_collapsed = kappa_phase;
PARAM.kappa_p_phase         = p_phase;
PARAM.kappa_min_frac_used   = kappa_min_frac;

end


function B = DilateNoWrap(A)
%One-cell 4-neighbour dilation without periodic wrapping.

[ny,nx] = size(A);
B = A;

if ny > 1
    B(2:ny,:) = B(2:ny,:) | A(1:ny-1,:);
    B(1:ny-1,:) = B(1:ny-1,:) | A(2:ny,:);
end

if nx > 1
    B(:,2:nx) = B(:,2:nx) | A(:,1:nx-1);
    B(:,1:nx-1) = B(:,1:nx-1) | A(:,2:nx);
end

end
