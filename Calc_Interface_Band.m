function mask_band = Calc_Interface_Band(STATE,MODEL,NUM)
%CALC_INTERFACE_BAND Find geometric phase-interface band.
%
% Includes both mixed cells and neighbour jumps. This catches sharp or
% almost-pure interface rims where p is locally close to 0 or 1.

[ny,nx,Ng] = size(STATE.p);

p_cut = 1e-6;
band_thick = 2;

if isfield(NUM,'interface_p_cut')
    p_cut = NUM.interface_p_cut;
end

if isfield(NUM,'interface_band_thick')
    band_thick = NUM.interface_band_thick;
end

if isfield(MODEL,'phase_index') && ~isempty(MODEL.phase_index)
    phase_ids = unique(MODEL.phase_index,'stable');
else
    phase_ids = 1:Ng;
    MODEL.phase_index = 1:Ng;
end

Np = numel(phase_ids);
p_phase = zeros(ny,nx,Np);

for iph = 1:Np
    grains = find(MODEL.phase_index == phase_ids(iph));
    p_phase(:,:,iph) = sum(STATE.p(:,:,grains),3);
end

mask_core = false(ny,nx);

for ip = 1:Np
    p = p_phase(:,:,ip);

    mixed = p > p_cut & p < 1-p_cut;

    if nx == 1
        p_L = p;
        p_R = p;
    else
        p_L = p(:,[2,1:nx-1]);
        p_R = p(:,[2:nx,nx-1]);
    end

    if ny == 1
        p_U = p;
        p_D = p;
    else
        p_U = p([2,1:ny-1],:);
        p_D = p([2:ny,ny-1],:);
    end

    jump = abs(p-p_L) > p_cut | abs(p-p_R) > p_cut | ...
           abs(p-p_U) > p_cut | abs(p-p_D) > p_cut;

    mask_core = mask_core | mixed | jump;
end

mask_band = mask_core;

for k = 1:band_thick
    if nx == 1
        L = mask_band;
        R = mask_band;
    else
        L = mask_band(:,[2,1:nx-1]);
        R = mask_band(:,[2:nx,nx-1]);
    end

    if ny == 1
        U = mask_band;
        D = mask_band;
    else
        U = mask_band([2,1:ny-1],:);
        D = mask_band([2:ny,ny-1],:);
    end

    mask_band = mask_band | L | R | U | D;
end

end