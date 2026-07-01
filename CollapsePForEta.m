function p_phase = CollapsePForEta(p_grain,phase_index)
%COLLAPSEPFORETA Collapse grain phase fields into thermodynamic phase fields.
%
% p_grain    : ny x nx x Ngrain
% phase_index: 1 x Ngrain, thermodynamic phase id for each grain

[ny,nx,Ng] = size(p_grain);

phase_index = phase_index(:).';
phase_id    = unique(phase_index,'stable');
Nphase      = numel(phase_id);

p_phase = zeros(ny,nx,Nphase);

for iph = 1:Nphase
    grains = find(phase_index == phase_id(iph));
    p_phase(:,:,iph) = sum(p_grain(:,:,grains),3);
end

end