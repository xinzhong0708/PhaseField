function [phase_id,phase_prop] = PhaseProportion_Local(p,MODEL)

[ny,nx,Ngrain] = size(p);
Nnode          = ny*nx;

if isfield(MODEL,'phase_index') && numel(MODEL.phase_index) == Ngrain
    phase_index = MODEL.phase_index(:).';
else
    phase_index = 1:Ngrain;
end

phase_id   = unique(phase_index,'stable');
phase_prop = zeros(1,numel(phase_id));

for iph = 1:numel(phase_id)

    grains = find(phase_index == phase_id(iph));

    tmp = zeros(ny,nx);

    for ig = grains
        tmp = tmp + p(:,:,ig);
    end

    phase_prop(iph) = sum(tmp(:))/Nnode;

end

end