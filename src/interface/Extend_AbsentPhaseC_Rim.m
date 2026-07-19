function STATE = Extend_AbsentPhaseC_Rim(STATE,PARAM)
%EXTEND_ABSENTPHASEC_RIM Fast nearest-neighbour extension of phase c.
%
% Present-phase region is kept fixed.
% Only absent rim / optional bad points are modified.
% STATE.c_ext is written for kappa(c) solvers.

[ny,nx,Ng] = size(STATE.p);

rim_width   = 3;
rim_smooth  = 0;
p_core      = 1e-5;
fix_bad_all = 0;

if isfield(PARAM,'c_rim_width')
    rim_width = PARAM.c_rim_width;
end
if isfield(PARAM,'c_rim_smooth')
    rim_smooth = PARAM.c_rim_smooth;
end
if isfield(PARAM,'c_rim_p_core')
    p_core = PARAM.c_rim_p_core;
end
if isfield(PARAM,'c_rim_fix_bad_all')
    fix_bad_all = PARAM.c_rim_fix_bad_all;
end

rim_width  = max(0,round(rim_width));
rim_smooth = max(0,round(rim_smooth));

if exist('bwdist','file') ~= 2
    error('Extend_AbsentPhaseC_Rim: bwdist is required for this fast version.')
end

c_ext = STATE.c;
ker   = ones(2*rim_width+1,2*rim_width+1);

for ig = 1:Ng

    present = STATE.p(:,:,ig) > p_core;

    if ~any(present(:))
        continue
    end

    if rim_width > 0
        band = conv2(double(present),ker,'same') > 0;
    else
        band = present;
    end

    rim = band & ~present;

    if ~any(rim(:)) && fix_bad_all == 0
        continue
    end

    % ------------------------------------------------------------
    % Stack all composition variables
    % ------------------------------------------------------------
    Nc = numel(STATE.c{ig});
    C0 = cat(3,STATE.c{ig}{:});
    C  = C0;

    % Valid source pixels must be inside the phase and finite for all c
    valid = present & all(isfinite(C0),3);

    if ~any(valid(:))
        continue
    end

    % ------------------------------------------------------------
    % Nearest present-phase source for all target points
    % ------------------------------------------------------------
    [~,idx] = bwdist(valid);

    target = rim;

    if fix_bad_all == 1
        target = target | any(~isfinite(C0),3);
    end

    ids = find(target);

    if ~isempty(ids)

        ids = double(ids(:));
        src = double(idx(ids));

        Nxy  = ny*nx;
        offs = double((0:Nc-1)*Nxy);

        linT = ids + offs;
        linS = src + offs;

        C(linT) = C(linS);

    end

    % ------------------------------------------------------------
    % Optional smoothing only in the rim
    % ------------------------------------------------------------
    if rim_smooth > 0 && any(rim(:))

        for it = 1:rim_smooth
            C = Smooth_Rim_Local(C,C0,present,rim,band);
        end

    end

    % ------------------------------------------------------------
    % Assign back
    % ------------------------------------------------------------
    for ic = 1:Nc
        c_ext{ig}{ic} = C(:,:,ic);
    end

end

STATE.c     = c_ext;
STATE.c_ext = c_ext;

end


function C = Smooth_Rim_Local(C,C0,present,rim,band)

[ny,nx,Nc] = size(C);

if nx == 1
    CL = C; CR = C;
    BL = band; BR = band;
else
    CL = C(:,[2,1:nx-1],:);
    CR = C(:,[2:nx,nx-1],:);
    BL = band(:,[2,1:nx-1]);
    BR = band(:,[2:nx,nx-1]);
end

if ny == 1
    CU = C; CD = C;
    BU = band; BD = band;
else
    CU = C([2,1:ny-1],:,:);
    CD = C([2:ny,ny-1],:,:);
    BU = band([2,1:ny-1],:);
    BD = band([2:ny,ny-1],:);
end

BL3 = repmat(BL,1,1,Nc);
BR3 = repmat(BR,1,1,Nc);
BU3 = repmat(BU,1,1,Nc);
BD3 = repmat(BD,1,1,Nc);

CL(~BL3) = C(~BL3);
CR(~BR3) = C(~BR3);
CU(~BU3) = C(~BU3);
CD(~BD3) = C(~BD3);

Avg = 0.25*(CL + CR + CU + CD);

R3 = repmat(rim,1,1,Nc);
P3 = repmat(present,1,1,Nc);

C(R3) = Avg(R3);
C(P3) = C0(P3);

end