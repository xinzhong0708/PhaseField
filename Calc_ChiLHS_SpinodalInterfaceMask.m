function PARAM = Calc_ChiLHS_SpinodalInterfaceMask(STATE,PARAM,MODEL)
%CALC_CHILHS_SPINODALINTERFACEMASK Build interface mask for chi_lhs.
%
% This is only for solver-side chi stabilization in ACCH/corrector.
% It does not change STATE.chi and does not change LE/GP thermodynamics.
%
% Rule:
%   selected spinodal phase near interface : mask = 1
%   selected spinodal phase core           : mask = 0
%   all other good phases                  : mask = 0
%
% Typical use:
%   PARAM.chi_lhs_phase_name  = 'Cpx';
%   PARAM.chi_lhs_p_interface = 0.999;
%   PARAM.chi_lhs_p_presence  = 1e-5;
%   PARAM.chi_lhs_width       = 3;
%   PARAM = Calc_ChiLHS_SpinodalInterfaceMask(STATE_OLD,PARAM,MODEL);

[ny,nx,Ng] = size(STATE.p);

%Defaults
p_interface = 0.999;
p_presence  = 1e-5;
width        = 3;

if isfield(PARAM,'chi_lhs_p_interface')
    p_interface = PARAM.chi_lhs_p_interface;
end

if isfield(PARAM,'chi_lhs_p_presence')
    p_presence = PARAM.chi_lhs_p_presence;
end

if isfield(PARAM,'chi_lhs_width')
    width = PARAM.chi_lhs_width;
end

width = max(0,round(width));

%Collapse grains to thermodynamic phases
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

%Choose which thermodynamic phases are allowed to be protected
use_phase = false(1,Nphase);

if isfield(PARAM,'chi_lhs_phase') && ~isempty(PARAM.chi_lhs_phase)

    phase_in = PARAM.chi_lhs_phase(:).';

    if numel(phase_in) == Nphase
        use_phase = phase_in ~= 0;

    elseif numel(phase_in) == Ng
        for iph = 1:Nphase
            grains = find(grain_to_phase == iph);
            use_phase(iph) = any(phase_in(grains) ~= 0);
        end
    end

elseif isfield(PARAM,'chi_lhs_phase_name') && isfield(MODEL,'phase_name')

    names = PARAM.chi_lhs_phase_name;

    if ischar(names) || isstring(names)
        names = cellstr(names);
    end

    for iph = 1:Nphase
        grains = find(grain_to_phase == iph);
        ig0 = grains(1);

        for iname = 1:numel(names)
            if strcmpi(MODEL.phase_name{ig0},names{iname})
                use_phase(iph) = true;
            end
        end
    end

else

    %No phase was specified. Do not protect anything by default.
    use_phase(:) = false;

end

%Build mask only for selected phases
mask_phase = zeros(ny,nx,Nphase);

for iph = 1:Nphase

    if ~use_phase(iph)
        continue
    end

    pcur = p_phase(:,:,iph);
    present = pcur > p_presence;

    %Interface/outside set for this phase
    interface = pcur < p_interface;

    if ~any(interface(:)) || all(interface(:))
        protect = false(ny,nx);

    else
        protect = interface;

        for ir = 1:width
            protect = DilateNoWrap(protect);
        end

        %Only inside the selected phase where it is actually present
        protect = protect & present;
    end

    mask_phase(:,:,iph) = double(protect);
end

PARAM.chi_lhs_mask_phase = mask_phase;
PARAM.chi_lhs_mask       = max(mask_phase,[],3);
PARAM.chi_lhs_p_phase    = p_phase;
PARAM.chi_lhs_use_phase  = use_phase;

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
