function PARAM = Compute_M_And_L(STATE,PARAM,MODEL,PHYS)
%COMPUTE_M_AND_L Build CH mobility and Allen-Cahn mobility.
%
% This function builds:
%
%   PARAM.M{ie,ie}   true CH diffusion mobility
%   PARAM.M{ie,je}   zero for ie ~= je
%
%   PARAM.L          Allen-Cahn mobility
%   PARAM.Lm         PARAM.L * PHYS.m
%   PARAM.LK         PARAM.L * PHYS.kap
%
% Mobility rule:
%
%   M = sum_i p_i M_i / sum_i p_i
%
% No log average.
% No harmonic average.
% No interface CH mobility floor.
% No artificial widening of low-M regions.
%
% PARAM.M is used by the CH solver.
% The mobility used to compute L can optionally be floored through
% PARAM.M_L_floor_fac, but by default it is zero.
PARAM.M_weight_mode = 'smoothstep';
M_weight_mode       = 'smoothstep';
PARAM.M_weight_mode = 'quintic';
if isfield(PARAM,'M_weight_mode')
    M_weight_mode = PARAM.M_weight_mode;
end


[ny,nx,~] = size(STATE.p);

if isfield(PARAM,'Ne')
    Ne = PARAM.Ne;
else
    Ne = numel(STATE.E);
end

if isfield(MODEL,'phs_name')
    Nphase = numel(MODEL.phs_name);
else
    Nphase = max(MODEL.phase_index);
end

if ~isfield(PHYS,'M_phs')
    error('Compute_M_And_L: PHYS.M_phs is missing.')
end

% ------------------------------------------------------------
% Defaults
% ------------------------------------------------------------
L_fac              = 1.0;
dceq               = 0.2;
M_min              = 1e-30;
z_min              = 1e-30;
M_L_floor_fac      = 0;
M_L_interface_only = 1;
M_L_p_cut          = 1e-8;

if isfield(PARAM,'L_fac'),              L_fac              = PARAM.L_fac;              end
if isfield(PHYS,'dceq'),                dceq               = PHYS.dceq;                end
if isfield(PARAM,'dceq'),               dceq               = PARAM.dceq;               end
if isfield(PARAM,'M_min'),              M_min              = PARAM.M_min;              end
if isfield(PARAM,'L_z_min'),            z_min              = PARAM.L_z_min;            end
if isfield(PARAM,'M_L_floor_fac'),      M_L_floor_fac      = PARAM.M_L_floor_fac;      end
if isfield(PHYS,'M_L_floor_fac'),       M_L_floor_fac      = PHYS.M_L_floor_fac;       end
if isfield(PARAM,'M_L_interface_only'), M_L_interface_only = PARAM.M_L_interface_only; end
if isfield(PARAM,'M_L_p_cut'),          M_L_p_cut          = PARAM.M_L_p_cut;          end

if isfield(PHYS,'kap') && ~isempty(PHYS.kap)
    kap = PHYS.kap;
else
    kap = PHYS.kappa;
end

if isscalar(dceq)
    dceq = dceq*ones(1,Ne);
end

if numel(dceq) ~= Ne
    error('Compute_M_And_L: dceq must be scalar or length Ne.')
end

if isscalar(M_L_floor_fac)
    M_L_floor_fac = M_L_floor_fac*ones(1,Ne);
end

if numel(M_L_floor_fac) ~= Ne
    error('Compute_M_And_L: M_L_floor_fac must be scalar or length Ne.')
end

% ------------------------------------------------------------
% Interpret PHYS.M_phs
% ------------------------------------------------------------
Mraw = PHYS.M_phs;

if isscalar(Mraw)

    M_phase_elem = Mraw*ones(Nphase,Ne);

elseif isvector(Mraw)

    Mraw = Mraw(:);

    if numel(Mraw) == Nphase

        M_phase_elem = repmat(Mraw,1,Ne);

    elseif numel(Mraw) == Ne

        M_phase_elem = repmat(Mraw(:).',Nphase,1);

    else

        error('Compute_M_And_L: vector PHYS.M_phs must have length Nphase or Ne.')

    end

else

    [n1,n2] = size(Mraw);

    if n1 == Nphase && n2 == Ne

        M_phase_elem = Mraw;

    elseif n1 == Ne && n2 == Nphase

        M_phase_elem = Mraw.';

    else

        error('Compute_M_And_L: PHYS.M_phs must be scalar, Nphase x Ne, Ne x Nphase, length Nphase, or length Ne.')

    end
end

M_phase_elem = max(M_phase_elem,M_min);

% ------------------------------------------------------------
% Collapse grain p to phase p
% ------------------------------------------------------------
p_phase = zeros(ny,nx,Nphase);

for iph = 1:Nphase

    grains = find(MODEL.phase_index == iph);

    if ~isempty(grains)
        p_phase(:,:,iph) = sum(STATE.p(:,:,grains),3);
    end
end

% ------------------------------------------------------------
% Mobility interpolation weight
%
% linear:
%   w = p
%
% smoothstep:
%   w = p^2*(3 - 2p)
% This keeps pure phases insensitive to tiny numerical phase tails.
% ------------------------------------------------------------
switch lower(M_weight_mode)

    case 'linear'

        p_mob = p_phase;

    case 'smoothstep'

        p_mob = p_phase.^2 .* (3 - 2*p_phase);

    case 'quintic'

        p_mob = p_phase.^3 .* (10 - 15*p_phase + 6*p_phase.^2);

    otherwise

        error('Compute_M_And_L: unknown PARAM.M_weight_mode "%s".',M_weight_mode)

end

psum     = sum(p_phase,3);
psum_mob = sum(p_mob,3);

% Interface mask from phase mixing.
% This affects only the optional L floor, not PARAM.M.
interface_mask = (1 - sum(p_phase.^2,3)) > M_L_p_cut;

% ------------------------------------------------------------
% Allocate true CH mobility
% ------------------------------------------------------------
PARAM.M = cell(Ne,Ne);

for ie = 1:Ne
    for je = 1:Ne
        PARAM.M{ie,je} = zeros(ny,nx);
    end
end

PARAM.M_diag = cell(1,Ne);

% ------------------------------------------------------------
% Build PARAM.M and compute zeta for L
% ------------------------------------------------------------
zeta = zeros(ny,nx);

for ie = 1:Ne

    % --------------------------------------------------------
    % Honest linear p-weighted mobility
    % --------------------------------------------------------
    Mgrid = zeros(ny,nx);

    for iph = 1:Nphase
        Mgrid = Mgrid + M_phase_elem(iph,ie).*p_mob(:,:,iph);
    end

    mask = psum_mob > eps;

    Mtmp = Mgrid;
    Mtmp(mask)  = Mgrid(mask)./psum_mob(mask);
    Mtmp(~mask) = mean(M_phase_elem(:,ie));

    Mtmp = max(Mtmp,M_min);

    % True CH diffusion mobility
    PARAM.M{ie,ie}   = Mtmp;
    PARAM.M_diag{ie} = Mtmp;

    % --------------------------------------------------------
    % Mobility used only to compute L
    % --------------------------------------------------------
    Mtmp_L = Mtmp;

    if isfield(PHYS,'M0') && ~isempty(PHYS.M0)

        if isscalar(PHYS.M0)
            Mref = PHYS.M0;
        elseif numel(PHYS.M0) == Ne
            Mref = PHYS.M0(ie);
        else
            Mref = mean(M_phase_elem(:,ie));
        end

    else

        Mref = mean(M_phase_elem(:,ie));

    end

    M_floor = Mref*M_L_floor_fac(ie);

    if M_floor > 0

        if M_L_interface_only == 1
            Mtmp_L(interface_mask) = max(Mtmp_L(interface_mask),M_floor);
        else
            Mtmp_L = max(Mtmp_L,M_floor);
        end

    end

    Mtmp_L = max(Mtmp_L,M_min);

    zeta = zeta + dceq(ie)^2 ./ Mtmp_L;

end

% ------------------------------------------------------------
% Compute Allen-Cahn mobility
% ------------------------------------------------------------
zeta = max(zeta,z_min);

PARAM.L  = L_fac * 4*PHYS.m ./ (3*kap.*zeta);
PARAM.Lm = PARAM.L .* PHYS.m;
PARAM.LK = PARAM.L .* kap;

PARAM.L(:)  = mean(PARAM.L( :));
PARAM.Lm(:) = mean(PARAM.Lm(:));
PARAM.LK(:) = mean(PARAM.LK(:));

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
PARAM.M_phase              = p_phase;
PARAM.M_phase_elem         = M_phase_elem;
PARAM.M_L_floor_fac        = M_L_floor_fac;
PARAM.M_L_interface_only   = M_L_interface_only;
PARAM.M_L_p_cut            = M_L_p_cut;

end