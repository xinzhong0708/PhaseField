function PARAM = Calc_WScale_InterfaceStable_AllPhases(STATE,PARAM,MODEL,PHYS,NUM)
%CALC_WSCALE_INTERFACESTABLE_ALLPHASES Smooth phase-wise W scaling at interfaces.
%
% Goal:
%   Scale solution-model excess energy W only near phase interfaces.
%   Keep WScale = 1 in pure interiors.
%   Use phase-wise damping factors from metadata, not pair-wise local choices.
%
% Output:
%   PARAM.WScale(:,:,ig) for every grain ig.
%
% Suggested metadata/parameter fields:
%   PARAM.Interface_damp_W_factor   one value per thermodynamic phase
%   PARAM.WScale_interface_factor   same meaning
%   PARAM.WScale_min                fallback scalar or vector
%
% Optional controls:
%   PARAM.WScale_p_tail        default 1e-5
%   PARAM.WScale_p_pure        default 1e-3
%   PARAM.WScale_halo          default 3
%   PARAM.WScale_smooth        default 1
%   PARAM.WScale_temporal      default 0.35
%   PARAM.WScale_enable        default 1
%
% Notes:
%   This function only constructs PARAM.WScale. PhaseThermo/LE must actually
%   use PARAM.WScale when building the phase thermodynamics. Do not multiply
%   omega directly after LE.

if isfield(PARAM,'WScale_enable') && PARAM.WScale_enable == 0
    PARAM.WScale = ones(size(STATE.p));
    return
end

[ny,nx,Ng] = size(STATE.p);

phase_index = MODEL.phase_index(:).';
phase_ids   = unique(phase_index,'stable');
Nphase      = numel(phase_ids);

p_tail   = 1e-5;
p_pure   = 1e-3;
nhalo    = 3;
nsmooth  = 1;
temporal = 0.35;

if isfield(PARAM,'WScale_p_tail'),   p_tail  = PARAM.WScale_p_tail;   end
if isfield(PARAM,'WScale_p_pure'),   p_pure  = PARAM.WScale_p_pure;   end
if isfield(PARAM,'WScale_halo'),     nhalo   = PARAM.WScale_halo;     end
if isfield(PARAM,'WScale_smooth'),   nsmooth = PARAM.WScale_smooth;   end
if isfield(PARAM,'WScale_temporal'), temporal = PARAM.WScale_temporal; end

f_phase = Get_W_Factor_Local(PARAM,Nphase,phase_ids);

% Safety bounds
f_phase = max(min(f_phase,1),0);

% Collapse p from grains to thermodynamic phases
p_phase = zeros(ny,nx,Nphase);

for iph = 1:Nphase
    grains = find(phase_index == phase_ids(iph));
    p_phase(:,:,iph) = sum(STATE.p(:,:,grains),3);
end

Wnew = ones(ny,nx,Ng);

for iph = 1:Nphase

    pph = p_phase(:,:,iph);

    % Interface core for this thermodynamic phase against all other phases.
    % This does not respond to same-phase grain boundaries.
    core = pph > p_tail & pph < 1 - p_pure;

    % Also catch sharp transitions that may be missed by a very small core.
    jump = Local_Jump_Mask(pph,p_tail);
    core = core | jump;

    % Smooth halo around interface.
    h = Local_Smooth_Halo(core,nhalo,nsmooth);

    % Stronger damping only at/near interface, exactly 1 far away.
    Wphase = 1 - (1 - f_phase(iph)).*h;
    Wphase = max(min(Wphase,1),f_phase(iph));

    grains = find(phase_index == phase_ids(iph));

    for ig = grains
        Wnew(:,:,ig) = Wphase;
    end
end

% Optional temporal smoothing to avoid WScale jumping when the interface moves.
if isfield(PARAM,'WScale') && isequal(size(PARAM.WScale),size(Wnew)) && temporal > 0
    temporal = max(min(temporal,1),0);
    Wnew = temporal*Wnew + (1-temporal)*PARAM.WScale;
end

PARAM.WScale = Wnew;

end


function f_phase = Get_W_Factor_Local(PARAM,Nphase,phase_ids)

f_phase = ones(1,Nphase);

names = {'Interface_damp_W_factor', ...
         'WScale_interface_factor', ...
         'WScale_min', ...
         'W_damp_factor', ...
         'Wfactor_interface'};

val = [];

for i = 1:numel(names)
    if isfield(PARAM,names{i}) && ~isempty(PARAM.(names{i}))
        val = PARAM.(names{i});
        break
    end
end

if isempty(val)
    return
end

val = val(:).';

if isscalar(val)

    f_phase(:) = val;

elseif numel(val) == Nphase

    f_phase = val;

elseif numel(val) >= max(phase_ids)

    for iph = 1:Nphase
        f_phase(iph) = val(phase_ids(iph));
    end

else

    n = min(numel(val),Nphase);
    f_phase(1:n) = val(1:n);

end

end


function jump = Local_Jump_Mask(q,tol)

[ny,nx] = size(q);

if nx == 1
    qL = q;
    qR = q;
else
    qL = q(:,[2,1:nx-1]);
    qR = q(:,[2:nx,nx-1]);
end

if ny == 1
    qU = q;
    qD = q;
else
    qU = q([2,1:ny-1],:);
    qD = q([2:ny,ny-1],:);
end

jump = abs(q-qL) > tol | abs(q-qR) > tol | ...
       abs(q-qU) > tol | abs(q-qD) > tol;

end


function h = Local_Smooth_Halo(core,nhalo,nsmooth)

if nhalo <= 0
    h = double(core);
else
    ker = ones(2*nhalo+1,2*nhalo+1);
    h   = conv2(double(core),ker,'same') > 0;
    h   = double(h);
end

% Smooth the sharp halo edge. Keep simple box smoothing.
for i = 1:max(0,nsmooth)
    h = conv2(h,ones(3,3)/9,'same');
end

h = max(min(h,1),0);

end
