function [PARAM] = Calc_Kappa_WScale_InterfaceRamp(STATE,PARAM,MODEL,PHYS,NUM)
%CALC_KAPPA_WSCALE_INTERFACERAMP Coordinate kappa and excess-energy scaling.
%
% This function does NOT change chi directly.
%
% It builds:
%   PARAM.kappa_eff      : spatial kappa map
%   PARAM.w_scale_phase  : excess-energy scale map for LE/GP thermodynamics
%
% Interpretation:
%   interface/protected layer:
%       kappa is reduced but nonzero inside a present selected phase
%       excess energy can be weakened for selected phases only
%
%   grain core:
%       kappa is full
%       excess energy is full
%
% Typical settings:
%   PARAM.kappa_phase       = PHYS.kappa .* cellfun(@(x) size(x.n,1) > 1, MODEL.pars);
%   PARAM.ramp_p_interface  = 0.999;
%   PARAM.ramp_p_presence   = 1e-5;
%   PARAM.ramp_zero_width   = 1;
%   PARAM.ramp_width        = 4;
%   PARAM.kappa_min_frac    = 0.02;
%   PARAM.w_scale_min       = 0.75;
%
% WScale phase selection:
%   default: damp all phases, preserving old behavior
%
%   PARAM.WScale_phase_name = 'Cpx';
%   PARAM.WScale_phase_name = {'Cpx','Opx'};
%   PARAM.WScale_phase_id   = [1 3];       % collapsed thermodynamic phase id
%   PARAM.WScale_grain_id   = [2 5];       % original grain id
%   PARAM.WScale_phase_mask = logical/numeric vector, length Ng or Nphase
%
% Optional expensive check:
%   PARAM.WScale_auto_posdef = 1;
%   PARAM.WScale_h_tol       = 1e-10;
%
% Required LE_Run_Mode_New behavior:
%   use PARAM.w_scale_phase to set pars{iph}.w_scale after phase collapse.
%
% NUM is currently unused but kept in the input list for compatibility.

[ny,nx,Ng] = size(STATE.p);

% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------
p_interface    = 0.999;
p_presence     = 1e-5;
zero_width     = 4;
ramp_width     = 5;
kappa_min_frac = 0.1;
w_scale_min    = 1;

if isfield(PARAM,'ramp_p_interface')
    p_interface = PARAM.ramp_p_interface;
elseif isfield(PARAM,'kappa_p_interface')
    p_interface = PARAM.kappa_p_interface;
end

if isfield(PARAM,'ramp_p_presence')
    p_presence = PARAM.ramp_p_presence;
elseif isfield(PARAM,'kappa_p_presence')
    p_presence = PARAM.kappa_p_presence;
end

if isfield(PARAM,'ramp_zero_width')
    zero_width = PARAM.ramp_zero_width;
elseif isfield(PARAM,'kappa_zero_width')
    zero_width = PARAM.kappa_zero_width;
end

if isfield(PARAM,'ramp_width')
    ramp_width = PARAM.ramp_width;
elseif isfield(PARAM,'kappa_ramp_width')
    ramp_width = PARAM.kappa_ramp_width;
end

if isfield(PARAM,'kappa_min_frac')
    kappa_min_frac = PARAM.kappa_min_frac;
end

if isfield(PARAM,'w_scale_min')
    w_scale_min = PARAM.w_scale_min;
end

zero_width     = max(0,round(zero_width));
ramp_width     = max(1,round(ramp_width));
kappa_min_frac = min(max(kappa_min_frac,0),1);
w_scale_min    = min(max(w_scale_min,0),1);

% -------------------------------------------------------------------------
% Collapse grains to thermodynamic phases
% -------------------------------------------------------------------------
phase_index = MODEL.phase_index(:).';
phase_id    = unique(phase_index,'stable');
Nphase      = numel(phase_id);

p_phase        = zeros(ny,nx,Nphase);
grain_to_phase = zeros(1,Ng);

for iph = 1:Nphase
    grains = find(phase_index == phase_id(iph));
    grain_to_phase(grains) = iph;
    p_phase(:,:,iph) = sum(STATE.p(:,:,grains),3);
end

% -------------------------------------------------------------------------
% Determine phase kappa after collapse
% -------------------------------------------------------------------------
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

    else

        error('PARAM.kappa_phase must be scalar, length Ng, or length Nphase.')

    end

elseif isfield(PARAM,'kappa_phase_name') && isfield(MODEL,'phs_name')

    for iph = 1:Nphase
        grains = find(grain_to_phase == iph);
        ig0 = grains(1);

        if strcmpi(MODEL.phs_name{ig0},PARAM.kappa_phase_name)
            kappa_phase(iph) = PHYS.kappa;
        end
    end

else

    kappa_phase(:) = PHYS.kappa;

end

% -------------------------------------------------------------------------
% Determine which phases receive WScale damping
%
% Default: damp all phases, preserving old behavior.
% If any WScale selector is given, only selected phases are damped.
% -------------------------------------------------------------------------
wscale_active = true(1,Nphase);

has_wscale_select = isfield(PARAM,'WScale_phase_name') || ...
                    isfield(PARAM,'WScale_phase_id')   || ...
                    isfield(PARAM,'WScale_grain_id')   || ...
                    isfield(PARAM,'WScale_phase_mask');

if has_wscale_select

    wscale_active = false(1,Nphase);

    % Select by phase name
    if isfield(PARAM,'WScale_phase_name')

        if ~isfield(MODEL,'phs_name')
            error('PARAM.WScale_phase_name requires MODEL.phase_name.')
        end

        names = PARAM.WScale_phase_name;

        if ischar(names) || isstring(names)
            names = cellstr(names);
        end

        for iph = 1:Nphase

            grains = find(grain_to_phase == iph);
            ig0    = grains(1);

            for iname = 1:numel(names)
                if strcmpi(MODEL.phs_name{ig0},names{iname})
                    wscale_active(iph) = true;
                end
            end

        end
    end

    % Select by collapsed thermodynamic phase id
    if isfield(PARAM,'WScale_phase_id')

        ids = PARAM.WScale_phase_id(:).';

        for ii = 1:numel(ids)
            if ids(ii) >= 1 && ids(ii) <= Nphase
                wscale_active(ids(ii)) = true;
            end
        end
    end

    % Select by original grain id
    if isfield(PARAM,'WScale_grain_id')

        gids = PARAM.WScale_grain_id(:).';

        for ii = 1:numel(gids)
            if gids(ii) >= 1 && gids(ii) <= Ng
                iph = grain_to_phase(gids(ii));
                wscale_active(iph) = true;
            end
        end
    end

    % Select by logical/numeric mask
    if isfield(PARAM,'WScale_phase_mask')

        mask_in = PARAM.WScale_phase_mask(:).';

        if isscalar(mask_in)

            wscale_active(:) = mask_in ~= 0;

        elseif numel(mask_in) == Ng

            for iph = 1:Nphase
                grains = find(grain_to_phase == iph);
                wscale_active(iph) = any(mask_in(grains) ~= 0);
            end

        elseif numel(mask_in) == Nphase

            wscale_active = mask_in ~= 0;

        else

            error('PARAM.WScale_phase_mask must be scalar, length Ng, or length Nphase.')

        end
    end
end

% -------------------------------------------------------------------------
% Build ramp maps
% -------------------------------------------------------------------------
weight_raw    = zeros(ny,nx,Nphase);
weight_kappa  = zeros(ny,nx,Nphase);
w_scale_phase = ones(ny,nx,Nphase);
kappa_eff     = zeros(ny,nx);

for iph = 1:Nphase

    pcur = p_phase(:,:,iph);
    kp   = kappa_phase(iph);

    present = pcur > p_presence;

    % Cells with p below p_interface define interface/outside.
    zero_mask = pcur < p_interface;

    if ~any(zero_mask(:))

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

    w_raw(~present) = 0;

    % Kappa: reduced but not zero inside present phase.
    wk = zeros(ny,nx);

    if kp ~= 0
        wk(present) = kappa_min_frac + (1-kappa_min_frac).*w_raw(present);
    end

    % Excess-energy scale:
    %   selected phase     -> weakened near interface, full in core
    %   non-selected phase -> raw thermodynamics everywhere
    ws = ones(ny,nx);

    if wscale_active(iph)
        ws(present) = w_scale_min + (1-w_scale_min).*w_raw(present);
    end

    weight_raw(:,:,iph)    = w_raw;
    weight_kappa(:,:,iph)  = wk;
    w_scale_phase(:,:,iph) = ws;

    if kp ~= 0
        kappa_eff = max(kappa_eff,kp.*wk);
    end
end

% -------------------------------------------------------------------------
% Optional expensive positivity search
% Only selected phases can be changed because non-selected phases have ws=1.
% -------------------------------------------------------------------------
if isfield(PARAM,'WScale_auto_posdef') && PARAM.WScale_auto_posdef == 1
    w_scale_phase = Auto_Positive_WScale(STATE,PARAM,MODEL,p_phase,w_scale_phase,grain_to_phase,wscale_active);
end

% -------------------------------------------------------------------------
% Store diagnostics and output
% -------------------------------------------------------------------------
active_kappa_phase = kappa_phase ~= 0;

if any(active_kappa_phase)
    PARAM.kappa_raw_weight = max(weight_raw(:,:,active_kappa_phase),[],3);
    PARAM.kappa_chi_weight = max(weight_kappa(:,:,active_kappa_phase),[],3);
else
    PARAM.kappa_raw_weight = max(weight_raw,[],3);
    PARAM.kappa_chi_weight = max(weight_kappa,[],3);
end

PARAM.kappa_eff              = kappa_eff;
PARAM.kappa_weight_phase     = weight_kappa;
PARAM.kappa_weight_phase_raw = weight_raw;
PARAM.kappa_phase_collapsed  = kappa_phase;
PARAM.kappa_p_phase          = p_phase;

PARAM.w_scale_phase          = w_scale_phase;
PARAM.w_scale_active_phase   = wscale_active;
PARAM.w_scale_min_used       = w_scale_min;

PARAM.kappa_min_frac_used    = kappa_min_frac;
PARAM.use_WScale             = 1;

end


% =========================================================================
% Optional automatic WScale reduction
% =========================================================================
function w_scale_phase = Auto_Positive_WScale(STATE,PARAM,MODEL,p_phase,w_scale_phase,grain_to_phase,wscale_active)
%AUTO_POSITIVE_WSCALE Reduce w_scale until H_c is positive in protected layer.

[ny,nx,Nphase] = size(w_scale_phase);
N = ny*nx;

h_tol  = 1e-10;
max_it = 8;
shrink = 0.5;

if isfield(PARAM,'WScale_h_tol')
    h_tol = PARAM.WScale_h_tol;
end
if isfield(PARAM,'WScale_auto_iter')
    max_it = PARAM.WScale_auto_iter;
end
if isfield(PARAM,'WScale_auto_shrink')
    shrink = PARAM.WScale_auto_shrink;
end

c_phase = CollapseCToPhase(STATE,p_phase,grain_to_phase);

for iph = 1:Nphase

    if ~wscale_active(iph)
        continue
    end

    ig0 = find(grain_to_phase == iph,1,'first');

    if isempty(ig0) || isempty(c_phase{iph})
        continue
    end

    pars = MODEL.pars{ig0};
    ws   = w_scale_phase(:,:,iph);

    active = ws < 0.999999;

    if ~any(active(:))
        continue
    end

    for it = 1:max_it

        pars_tmp = pars;
        pars_tmp.w_scale = reshape(ws,1,N);

        R = PhaseThermo(pars_tmp,c_phase{iph});

        if isempty(R.H_c)
            break
        end

        mineig = reshape(MinEigPages(R.H_c),ny,nx);
        bad = mineig < h_tol;
        bad = bad & active;

        if ~any(bad(:))
            break
        end

        ws(bad) = shrink .* ws(bad);
    end

    w_scale_phase(:,:,iph) = ws;
end

end


% =========================================================================
% Collapse grain compositions to thermodynamic phases
% =========================================================================
function c_phase = CollapseCToPhase(STATE,p_phase,grain_to_phase)

[ny,nx,~] = size(STATE.p);
Nphase = size(p_phase,3);
N = ny*nx;

c_phase = cell(1,Nphase);

for iph = 1:Nphase

    grains = find(grain_to_phase == iph);

    if isempty(grains)
        c_phase{iph} = {};
        continue
    end

    ig0 = grains(1);
    Nc = numel(STATE.c{ig0});
    c_phase{iph} = cell(1,Nc);

    den = reshape(p_phase(:,:,iph),1,N);
    good = den > eps;

    for ic = 1:Nc

        num = zeros(1,N);

        for ig = grains
            num = num + reshape(STATE.p(:,:,ig),1,N).*reshape(STATE.c{ig}{ic},1,N);
        end

        tmp = reshape(STATE.c{ig0}{ic},1,N);
        tmp(good) = num(good)./den(good);

        c_phase{iph}{ic} = tmp;
    end
end

end


% =========================================================================
% Small helpers
% =========================================================================
function mineig = MinEigPages(H)

[~,~,N] = size(H);
mineig = zeros(1,N);

for i = 1:N
    A = 0.5*(H(:,:,i) + H(:,:,i).');

    if any(~isfinite(A(:)))
        mineig(i) = -inf;
    else
        mineig(i) = min(eig(A));
    end
end

end


function B = DilateNoWrap(A)

[ny,nx] = size(A);
B = A;

if ny > 1
    B(2:ny,:)   = B(2:ny,:)   | A(1:ny-1,:);
    B(1:ny-1,:) = B(1:ny-1,:) | A(2:ny,:);
end

if nx > 1
    B(:,2:nx)   = B(:,2:nx)   | A(:,1:nx-1);
    B(:,1:nx-1) = B(:,1:nx-1) | A(:,2:nx);
end

end