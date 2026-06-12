function STATE = Extend_AbsentPhaseC_Rim(STATE,PARAM)
%EXTEND_ABSENTPHASEC_RIM Extend phase c a few grids outside each grain.
%
% This function is meant to be called after LE_Run_Mode_New.
% It does not solve local equilibrium. It only fills the absent-phase rim
% from neighbouring present-phase c so that kappa(c) does not see a sharp
% artificial jump at the phase boundary.
%
% The present-phase region is kept fixed. Only the outer rim is modified.
% STATE.c_ext is also written and can be used by kappa(c) solvers.

[ny,nx,Ng] = size(STATE.p);

rim_width   = 3;
rim_smooth  = 1;
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

c_ext = STATE.c;

for ig = 1:Ng

    p_ig    = STATE.p(:,:,ig);
    present = p_ig > p_core;

    if ~any(present(:))
        continue
    end

    band = Dilate_Mask_Local(present,rim_width);
    rim  = band & ~present;

    if ~any(rim(:)) && fix_bad_all == 0
        continue
    end

    for ic = 1:numel(STATE.c{ig})

        A0 = STATE.c{ig}{ic};
        A  = A0;

        good = isfinite(A);
        valid = present & good;

        if ~any(valid(:))
            continue
        end

        %Fill the outer rim from neighbouring present/filled values.
        todo = rim & ~valid;

        if fix_bad_all == 1
            todo = todo | ~good;
        end

        for it = 1:(rim_width+2)

            [sumv,cnt] = Neighbour_Sum_Valid(A,valid);
            take = todo & cnt > 0;

            if ~any(take(:))
                break
            end

            A(take) = sumv(take)./cnt(take);
            valid(take) = true;
            todo(take)  = false;

            if ~any(todo(:))
                break
            end
        end

        %If requested, repair isolated bad points using the present-phase mean.
        if fix_bad_all == 1 && any(todo(:))
            A(todo) = mean(A(valid),'all');
            valid(todo) = true;
        end

        %Optional smoothing only in the rim. The grain interior is fixed.
        for it = 1:rim_smooth
            Avg = Neighbour_Avg_Domain(A,band);
            A(rim)     = Avg(rim);
            A(present) = A0(present);
        end

        c_ext{ig}{ic} = A;

    end
end

%Use the rim-extended c as both auxiliary output and safe phase-c memory.
%Only rim/optional bad nodes are modified; present-phase c is unchanged.
STATE.c     = c_ext;
STATE.c_ext = c_ext;

end


function mask = Dilate_Mask_Local(core,thickness)

if thickness <= 0
    mask = logical(core);
    return
end

ker  = ones(2*thickness+1,2*thickness+1);
mask = conv2(double(core),ker,'same') > 0;

end


function [sumv,cnt] = Neighbour_Sum_Valid(A,valid)

[ny,nx] = size(A);

if nx == 1
    AL = A; AR = A;
    VL = valid; VR = valid;
else
    AL = A(:,[2,1:nx-1]);
    AR = A(:,[2:nx,nx-1]);
    VL = valid(:,[2,1:nx-1]);
    VR = valid(:,[2:nx,nx-1]);
end

if ny == 1
    AU = A; AD = A;
    VU = valid; VD = valid;
else
    AU = A([2,1:ny-1],:);
    AD = A([2:ny,ny-1],:);
    VU = valid([2,1:ny-1],:);
    VD = valid([2:ny,ny-1],:);
end

sumv = zeros(ny,nx);
cnt  = zeros(ny,nx);

sumv = sumv + AL.*VL; cnt = cnt + VL;
sumv = sumv + AR.*VR; cnt = cnt + VR;
sumv = sumv + AU.*VU; cnt = cnt + VU;
sumv = sumv + AD.*VD; cnt = cnt + VD;

end


function Avg = Neighbour_Avg_Domain(A,domain)

if size(A,2) == 1
    AL = A; AR = A;
    DL = domain; DR = domain;
else
    AL = A(:,[2,1:size(A,2)-1]);
    AR = A(:,[2:size(A,2),size(A,2)-1]);
    DL = domain(:,[2,1:size(A,2)-1]);
    DR = domain(:,[2:size(A,2),size(A,2)-1]);
end

if size(A,1) == 1
    AU = A; AD = A;
    DU = domain; DD = domain;
else
    AU = A([2,1:size(A,1)-1],:);
    AD = A([2:size(A,1),size(A,1)-1],:);
    DU = domain([2,1:size(A,1)-1],:);
    DD = domain([2:size(A,1),size(A,1)-1],:);
end

AL(~DL) = A(~DL);
AR(~DR) = A(~DR);
AU(~DU) = A(~DU);
AD(~DD) = A(~DD);

Avg = 0.25*(AL+AR+AU+AD);

end
