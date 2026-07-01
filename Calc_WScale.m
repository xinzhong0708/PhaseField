function [PARAM] = Calc_WScale(STATE,PARAM,MODEL,PHYS,NUM)
%CALC_WSCALE Prepare spatial W scaling for LE_Run_Mode.
%
% Output required by LE_Run_Mode:
%   PARAM.w_scale_phase(:,:,iph)
%
% Here iph is the collapsed thermodynamic phase order:
%   phase_id = unique(MODEL.phase_index,'stable')
%
% This function does not scale omega. It only prepares pars.w_scale for
% PhaseThermo/PhaseG inside LE_Calculator.

[ny,nx,Ng]  = size(STATE.p);
phase_index = MODEL.phase_index(:).';
phase_id    = unique(phase_index,'stable');
Nphase      = numel(phase_id);

% ------------------------------------------------------------
% Parameters
% ------------------------------------------------------------
p_cut   = 1e-8;
p_jump  = p_cut;
nint    = 2;
nsmooth = 4;

% pseudo-1D switch:
%   0 = full 2D
%   1 = force x-only WScale and copy to all rows
pseudo1D = 0;

if isfield(PARAM,'LE_p_tail')
    p_cut = min(p_cut,0.01*PARAM.LE_p_tail);
end
if isfield(PARAM,'WScale_p_cut'),    p_cut   = PARAM.WScale_p_cut;    end
if isfield(PARAM,'WScale_p_jump'),   p_jump  = PARAM.WScale_p_jump;   end
if isfield(PARAM,'WScale_nint'),     nint    = PARAM.WScale_nint;     end
if isfield(PARAM,'WScale_smooth'),   nsmooth = PARAM.WScale_smooth;   end
if isfield(PARAM,'WScale_pseudo1D'), pseudo1D = PARAM.WScale_pseudo1D; end

p_cut   = max(p_cut,0);
p_jump  = max(p_jump,0);
nint    = max(round(nint),0);
nsmooth = max(round(nsmooth),0);

% Auto-detect pseudo-1D only for very thin domains
if ~isfield(PARAM,'WScale_pseudo1D')
    if ny <= 4
        pseudo1D = 1;
    end
end

% ------------------------------------------------------------
% Collapse grain p to thermodynamic phase p
% ------------------------------------------------------------
p_phase = zeros(ny,nx,Nphase);

for iph = 1:Nphase

    grains = find(phase_index == phase_id(iph));

    for ii = 1:numel(grains)
        ig = grains(ii);
        p_phase(:,:,iph) = p_phase(:,:,iph) + STATE.p(:,:,ig);
    end
end

% Normalize only for interface detection
psum    = sum(p_phase,3);
q_phase = p_phase ./ max(psum,eps);

% For pseudo-1D, force the interface detector to see one y-averaged profile.
% This prevents row 1/2/3/4 from getting different WScale.
if pseudo1D == 1
    qbar    = mean(q_phase,1);
    q_phase = repmat(qbar,ny,1,1);
end

% ------------------------------------------------------------
% Read phase-specific W scaling factor
% ------------------------------------------------------------
Wfac_phase = ones(1,Nphase);

if isfield(PARAM,'WScale_phase_factor') && ~isempty(PARAM.WScale_phase_factor)

    Wfac0 = PARAM.WScale_phase_factor(:).';
    Wfac0 = max(min(Wfac0,1),0);

    for iph = 1:Nphase

        pid = phase_id(iph);

        if pid <= numel(Wfac0)
            Wfac_phase(iph) = Wfac0(pid);
        elseif numel(Wfac0) == Nphase
            Wfac_phase(iph) = Wfac0(iph);
        else
            Wfac_phase(iph) = 1;
        end
    end
end

% ------------------------------------------------------------
% Detect phase-specific interface regions
% ------------------------------------------------------------
interface_phase = false(ny,nx,Nphase);
interface_all   = false(ny,nx);

for iph = 1:Nphase

    q = q_phase(:,:,iph);

    if pseudo1D == 1

        q1 = q(1,:);

        mixed1 = q1 > p_cut & q1 < 1 - p_cut;

        jump1 = false(1,nx);

        if nx > 1
            dq = abs(q1(2:end) - q1(1:end-1));
            jj = dq > p_jump;

            jump1(2:end)   = jump1(2:end)   | jj;
            jump1(1:end-1) = jump1(1:end-1) | jj;
        end

        core1 = mixed1 | jump1;
        core  = repmat(core1,ny,1);

    else

        % Diffuse coexistence/interface tail.
        mixed = q > p_cut & q < 1 - p_cut;

        % Sharp jump, marking both sides of a phase boundary.
        jump = false(ny,nx);

        if nx > 1
            dq = abs(q(:,2:end) - q(:,1:end-1));
            jj = dq > p_jump;

            jump(:,2:end)   = jump(:,2:end)   | jj;
            jump(:,1:end-1) = jump(:,1:end-1) | jj;
        end

        if ny > 1
            dq = abs(q(2:end,:) - q(1:end-1,:));
            jj = dq > p_jump;

            jump(2:end,:)   = jump(2:end,:)   | jj;
            jump(1:end-1,:) = jump(1:end-1,:) | jj;
        end

        core = mixed | jump;

    end

    interface_phase(:,:,iph) = core;
    interface_all            = interface_all | core;
end

% ------------------------------------------------------------
% Extend and smooth each phase interface region
% ------------------------------------------------------------
halo_phase = zeros(ny,nx,Nphase);

for iph = 1:Nphase

    core = interface_phase(:,:,iph);

    if pseudo1D == 1

        core1 = core(1,:);

        if nint > 0
            kernel1 = ones(1,2*nint+1);
            halo1   = conv2(double(core1),kernel1,'same') > 0;
            halo1   = double(halo1);
        else
            halo1   = double(core1);
        end

        for is = 1:nsmooth
            halo1 = Smooth1_Replicate(halo1);
        end

        halo = repmat(halo1,ny,1);

    else

        if nint > 0
            halo = Dilate2_Replicate(core,nint);
            halo = double(halo);
        else
            halo = double(core);
        end

        for is = 1:nsmooth
            halo = Smooth2_Replicate(halo);
        end

    end

    halo = max(min(halo,1),0);

    halo_phase(:,:,iph) = halo;
end

% ------------------------------------------------------------
% Build WScale in collapsed thermodynamic phase order
% ------------------------------------------------------------
w_scale_phase = ones(ny,nx,Nphase);

for iph = 1:Nphase

    fac  = Wfac_phase(iph);
    halo = halo_phase(:,:,iph);

    % Interior: 1
    % Interface halo: approaches fac
    w_scale_phase(:,:,iph) = 1 - (1 - fac).*halo;
    w_scale_phase(:,:,iph) = max(min(w_scale_phase(:,:,iph),1),fac);
end

% Final safety for pseudo-1D: force every row identical
if pseudo1D == 1
    w_scale_phase = repmat(mean(w_scale_phase,1),ny,1,1);
    halo_phase    = repmat(mean(halo_phase,1),ny,1,1);
end

% ------------------------------------------------------------
% Store variables used by LE_Run_Mode
% ------------------------------------------------------------
PARAM.w_scale_phase = w_scale_phase;
PARAM.use_WScale    = any(abs(w_scale_phase(:) - 1) > 1e-14);

% Useful diagnostics/plotting fields
PARAM.WScale_phase       = w_scale_phase;
PARAM.WScale_interface   = interface_all;
PARAM.WScale_halo_phase  = halo_phase;
PARAM.WScale_factor_used = Wfac_phase;
PARAM.WScale_pseudo1D    = pseudo1D;

PARAM.WScale_diag = struct();
PARAM.WScale_diag.phase_id      = phase_id;
PARAM.WScale_diag.factor        = Wfac_phase;
PARAM.WScale_diag.min_phase     = squeeze(min(min(w_scale_phase,[],1),[],2)).';
PARAM.WScale_diag.max_phase     = squeeze(max(max(w_scale_phase,[],1),[],2)).';
PARAM.WScale_diag.halo_fraction = squeeze(mean(mean(halo_phase > 1e-12,1),2)).';
PARAM.WScale_diag.p_cut         = p_cut;
PARAM.WScale_diag.p_jump        = p_jump;
PARAM.WScale_diag.nint          = nint;
PARAM.WScale_diag.nsmooth       = nsmooth;
PARAM.WScale_diag.pseudo1D      = pseudo1D;
PARAM.WScale_diag.use_WScale    = PARAM.use_WScale;

if isfield(PARAM,'WScale_debug') && PARAM.WScale_debug == 1

    fprintf('\nWScale diagnostic:\n')

    for iph = 1:Nphase

        pid = phase_id(iph);

        if isfield(MODEL,'phs_name') && pid <= numel(MODEL.phs_name)
            pname = MODEL.phs_name{pid};
        else
            pname = sprintf('phase%d',pid);
        end

        fprintf('  iph %d phase_id %d %-18s factor %.3e W %.3e to %.3e halo %.3e\n', ...
            iph,pid,pname,Wfac_phase(iph), ...
            PARAM.WScale_diag.min_phase(iph), ...
            PARAM.WScale_diag.max_phase(iph), ...
            PARAM.WScale_diag.halo_fraction(iph))
    end

    fprintf('  pseudo1D   = %d\n',pseudo1D)
    fprintf('  use_WScale = %d\n',PARAM.use_WScale)
end

end


% =========================================================================
% Local helpers
% =========================================================================
function B = Smooth1_Replicate(A)

nx = numel(A);
B  = zeros(size(A));

for dx = -1:1

    ix = (1:nx) + dx;
    ix = max(min(ix,nx),1);

    B = B + A(ix);
end

B = B/3;

end


function B = Smooth2_Replicate(A)

[ny,nx] = size(A);
B       = zeros(ny,nx);

for dy = -1:1

    iy = (1:ny) + dy;
    iy = max(min(iy,ny),1);

    for dx = -1:1

        ix = (1:nx) + dx;
        ix = max(min(ix,nx),1);

        B = B + A(iy,ix);
    end
end

B = B/9;

end


function B = Dilate2_Replicate(A,nint)

[ny,nx] = size(A);
B       = false(ny,nx);

for dy = -nint:nint

    iy = (1:ny) + dy;
    iy = max(min(iy,ny),1);

    for dx = -nint:nint

        ix = (1:nx) + dx;
        ix = max(min(ix,nx),1);

        B = B | A(iy,ix);
    end
end

end