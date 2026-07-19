function PARAM = Update_PF_SolverMasks(PARAM,STATE,MODEL,GRID,PHYS,NUM,mode)
%UPDATE_PF_SOLVERMASKS Build/reuse solver masks stored in PARAM.
%
% This minimal version only builds the ACCH phi mask used by
% PF_Coupled_ACCH_LETangent_CS through PARAM.maskPhi_ACCH.

if nargin < 7 || isempty(mode)
    mode = 'ACCH';
end

if ~strcmpi(mode,'ACCH')
    error('Update_PF_SolverMasks: only ACCH mode is implemented in this minimal version.')
end

phi = STATE.phi;
[ny,nx,Ngrain] = size(phi);

if isfield(NUM,'ACCH_mask_update') && ~isempty(NUM.ACCH_mask_update)
    update_every = NUM.ACCH_mask_update;
elseif isfield(NUM,'ilu_reuse_steps_ACCH') && ~isempty(NUM.ilu_reuse_steps_ACCH)
    update_every = NUM.ilu_reuse_steps_ACCH;
elseif isfield(NUM,'ilu_reuse_steps') && ~isempty(NUM.ilu_reuse_steps)
    update_every = NUM.ilu_reuse_steps;
else
    update_every = 1;
end

update_every = max(1,update_every);

rebuild = 0;

if ~isfield(PARAM,'maskPhi_ACCH') || isempty(PARAM.maskPhi_ACCH)
    rebuild = 1;
elseif ~isequal(size(PARAM.maskPhi_ACCH),[ny,nx,Ngrain])
    rebuild = 1;
end

if isfield(NUM,'istep') && ~isempty(NUM.istep)

    istep = NUM.istep;

    if ~isfield(PARAM,'maskPhi_ACCH_step') || isempty(PARAM.maskPhi_ACCH_step)
        rebuild = 1;
    elseif istep - PARAM.maskPhi_ACCH_step >= update_every
        rebuild = 1;
    end

else

    if ~isfield(PARAM,'maskPhi_ACCH_count') || isempty(PARAM.maskPhi_ACCH_count)
        PARAM.maskPhi_ACCH_count = update_every;
    end

    if PARAM.maskPhi_ACCH_count >= update_every
        rebuild = 1;
    end

end

if rebuild == 1

    if isfield(NUM,'phi_mask_cut')
        phi_cut = NUM.phi_mask_cut;
    else
        phi_cut = 1e-8;
    end

    if isfield(NUM,'phi_mask_pure_cut')
        pure_cut = NUM.phi_mask_pure_cut;
    else
        pure_cut = phi_cut;
    end

    if isfield(NUM,'phi_mask_thick')
        mask_thick = NUM.phi_mask_thick;
    else
        if (isfield(PHYS,'kap') && PHYS.kap ~= 0) || ...
                (isfield(PHYS,'kappa') && PHYS.kappa ~= 0)
            mask_thick = 2;
        else
            mask_thick = 1;
        end
    end

    maskPhi = Local_Calc_Interface_Mask(phi,phi_cut,pure_cut,mask_thick);

    if isfield(NUM,'phi_mask_source_tol') && ~isempty(NUM.phi_mask_source_tol)

        PARAM_MASK = PARAM;
        STATE_SRC  = STATE;

        if isfield(NUM,'use_Aac') && NUM.use_Aac == 1
            if isfield(NUM,'Aac_fac')
                Aac_fac = NUM.Aac_fac;
            else
                Aac_fac = 3;
            end
            PARAM_MASK.A_ac = Calc_Aac_FrozenOmega(STATE_SRC,PARAM_MASK,MODEL,Aac_fac,1e-6,0,[]);
        else
            PARAM_MASK.A_ac = zeros(ny,nx);
        end

        STATE_SRC = Calc_S_AllenCahn(STATE_SRC,PARAM_MASK,MODEL);
        S_AC      = STATE_SRC.S_AC;

        source_tol = NUM.phi_mask_source_tol;

        for alpha = 1:Ngrain
            core_source = abs(S_AC{alpha}) > source_tol;
            mask_source = Local_Dilate_Mask(core_source,mask_thick);
            maskPhi(:,:,alpha) = maskPhi(:,:,alpha) | mask_source;
        end

    end

    PARAM.maskPhi_ACCH = logical(maskPhi);

    if isfield(NUM,'istep') && ~isempty(NUM.istep)
        PARAM.maskPhi_ACCH_step = NUM.istep;
    else
        PARAM.maskPhi_ACCH_count = 0;
    end

else

    if ~isfield(NUM,'istep') || isempty(NUM.istep)
        PARAM.maskPhi_ACCH_count = PARAM.maskPhi_ACCH_count + 1;
    end

end

end


function mask = Local_Calc_Interface_Mask(phi,low_cut,pure_cut,thickness)

[ny,nx,Ngrain] = size(phi);

mask = false(ny,nx,Ngrain);

den = sum(phi.^2,3) + eps;

for alpha = 1:Ngrain

    q = phi(:,:,alpha).^2 ./ den;

    core = q > low_cut & q < 1 - pure_cut;

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

    jump = abs(q - qL) > low_cut | ...
           abs(q - qR) > low_cut | ...
           abs(q - qU) > low_cut | ...
           abs(q - qD) > low_cut;

    core = core | jump;

    mask(:,:,alpha) = Local_Dilate_Mask(core,thickness);

end

end


function mask = Local_Dilate_Mask(core,thickness)

if thickness <= 0
    mask = core;
    return
end

ker  = ones(2*thickness+1,2*thickness+1);
mask = conv2(double(core),ker,'same') > 0;

end
