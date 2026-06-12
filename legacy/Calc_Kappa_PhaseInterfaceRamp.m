function PARAM = Calc_Kappa_PhaseInterfaceRamp(STATE,PARAM,MODEL,PHYS)
%CALC_KAPPA_PHASEINTERFACERAMP Smooth kappa from interface to phase core.
%
% Rule:
%   interface / outside phase : kappa_eff = 0
%   transition band           : smooth ramp from 0 to full kappa
%   phase core                : kappa_eff = phase kappa
%
% This is intended to avoid applying the fourth-order CH regularization
% directly across phase interfaces, where it can push E outside the local
% equilibrium basin.
%
% Typical use:
%   PARAM.kappa_phase = PHYS.kappa .* cellfun(@(x) size(x.n,1) > 1, pars);
%   PARAM.kappa_zero_width = 1;
%   PARAM.kappa_ramp_width = 4;
%   PARAM.kappa_p_interface = 0.999;
%   PARAM = Calc_Kappa_PhaseInterfaceRamp(STATE,PARAM,MODEL,PHYS);
%
% If PARAM.kappa_phase is nonzero for all phases, this still sets kappa
% to zero at every phase interface and ramps it into each phase core.

[ny,nx,Ng] = size(STATE.p);

%Defaults
p_interface = 0.999;
p_presence  = 1e-12;
zero_width  = 1;
ramp_width  = 4;

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

zero_width = max(0,round(zero_width));
ramp_width = max(1,round(ramp_width));

%Collapse grains to thermodynamic phases
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

%Determine phase kappa after collapse
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
    %Fallback: apply PHYS.kappa to all thermodynamic phases.
    %This is safe because the interface ramp still forces kappa to zero at
    %phase boundaries.
    kappa_phase(:) = PHYS.kappa;
end

%Build effective kappa
kappa_eff = zeros(ny,nx);
weight_phase = zeros(ny,nx,Nphase);

for iph = 1:Nphase

    kp = kappa_phase(iph);

    if kp == 0
        continue
    end

    pcur = p_phase(:,:,iph);

    %Zero set for this phase:
    %  - outside the phase
    %  - diffuse-interface cells where p is not close enough to 1
    %
    %Then distance is propagated into the phase core.
    zero_mask = pcur < p_interface;

    %If the whole domain is this phase, there is no interface for this phase.
    if ~any(zero_mask(:))
        w = ones(ny,nx);

    else
        max_dist = zero_width + ramp_width + 1;
        dist = inf(ny,nx);
        known = zero_mask;
        dist(known) = 0;

        for ir = 1:max_dist
            known_new = DilateNoWrap(known);
            add = known_new & ~known;
            dist(add) = ir;
            known = known_new;
        end

        %Cells not reached within max_dist are deep core.
        dist(~isfinite(dist)) = max_dist;

        %Smooth ramp
        x = (dist - zero_width) ./ ramp_width;
        x = min(max(x,0),1);
        w = x.^2 .* (3 - 2*x);
    end

    %Do not give this phase kappa in places where it is truly absent.
    w(pcur <= p_presence) = 0;

    weight_phase(:,:,iph) = w;

    %If more than one phase has kappa, take the largest local contribution.
    %This keeps interface zero because all relevant phase weights are zero
    %near their boundaries.
    kappa_eff = max(kappa_eff,kp.*w);
end

PARAM.kappa_eff = kappa_eff;
PARAM.kappa_weight_phase = weight_phase;
PARAM.kappa_phase_collapsed = kappa_phase;
PARAM.kappa_p_phase = p_phase;

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
