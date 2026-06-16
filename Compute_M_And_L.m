function PARAM = Compute_M_And_L(STATE,PARAM,MODEL,PHYS)
%COMPUTE_M_AND_L Build CH mobility, L-mobility, and AC mobility.
%
% This function does three things:
%
%   1. Build PARAM.M from PHYS.M_phs
%      This is the true CH diffusion mobility.
%
%   2. Build PARAM.M_L from PARAM.M
%      This is used only for computing Allen-Cahn L.
%      It can have an interface floor so that slow phase diffusion does not
%      artificially freeze the phase boundary.
%
%   3. Compute PARAM.L, PARAM.Lm, PARAM.LK
%
% PHYS.M_phs can be:
%
%   scalar:
%       same mobility for all phases and elements
%
%   vector, length = Nphase:
%       same phase mobility for all elements
%
%   matrix, Nphase x Ne:
%       rows    = phases
%       columns = elements
%
%   matrix, Ne x Nphase:
%       rows    = elements
%       columns = phases
%
% Output:
%   PARAM.M{ie,ie}   = true diagonal CH mobility
%   PARAM.M{ie,je}   = zero for ie ~= je
%
%   PARAM.M_L{ie,ie} = mobility used only for L
%   PARAM.M_L{ie,je} = zero for ie ~= je
%
%   PARAM.L, PARAM.Lm, PARAM.LK

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

% ------------------------------------------------------------
% Defaults
% ------------------------------------------------------------
M_L_floor_fac      = 1e-2;
M_L_interface_only = 1;
M_L_p_cut          = 1e-8;

L_fac = 1.0;
dceq  = 0.5;
M_min = 1e-30;
z_min = 1e-30;

if isfield(PARAM,'M_L_floor_fac')
    M_L_floor_fac = PARAM.M_L_floor_fac;
end
if isfield(PHYS,'M_L_floor_fac')
    M_L_floor_fac = PHYS.M_L_floor_fac;
end
if isfield(PARAM,'M_L_interface_only')
    M_L_interface_only = PARAM.M_L_interface_only;
end
if isfield(PARAM,'M_L_p_cut')
    M_L_p_cut = PARAM.M_L_p_cut;
end

if isfield(PARAM,'L_fac')
    L_fac = PARAM.L_fac;
end
if isfield(PHYS,'dceq')
    dceq = PHYS.dceq;
end
if isfield(PARAM,'dceq')
    dceq = PARAM.dceq;
end
if isfield(PARAM,'M_min')
    M_min = PARAM.M_min;
end
if isfield(PARAM,'L_z_min')
    z_min = PARAM.L_z_min;
end

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
if ~isfield(PHYS,'M_phs')
    error('Compute_M_And_L: PHYS.M_phs is missing.')
end

Mraw = PHYS.M_phs;

if isscalar(Mraw)

    M_phase_elem = Mraw*ones(Nphase,Ne);

elseif isvector(Mraw)

    Mraw = Mraw(:);

    if numel(Mraw) ~= Nphase
        error('Compute_M_And_L: vector PHYS.M_phs must have length Nphase.')
    end

    M_phase_elem = repmat(Mraw,1,Ne);

else

    [n1,n2] = size(Mraw);

    if n1 == Nphase && n2 == Ne

        M_phase_elem = Mraw;

    elseif n1 == Ne && n2 == Nphase

        M_phase_elem = Mraw.';

    else

        error('Compute_M_And_L: PHYS.M_phs must be scalar, Nphase x Ne, or Ne x Nphase.')

    end
end

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

psum = sum(p_phase,3);

% Interface mask from phase mixture.
% Pure phase: sum(p^2) close to 1.
% Interface: sum(p^2) smaller than 1.
interface_mask = (1 - sum(p_phase.^2,3)) > M_L_p_cut;

% ------------------------------------------------------------
% Allocate mobility matrices
% ------------------------------------------------------------
PARAM.M   = cell(Ne,Ne);
PARAM.M_L = cell(Ne,Ne);

for ie = 1:Ne
    for je = 1:Ne
        PARAM.M{ie,je}   = zeros(ny,nx);
        PARAM.M_L{ie,je} = zeros(ny,nx);
    end
end

PARAM.M_diag   = cell(1,Ne);
PARAM.M_L_diag = cell(1,Ne);

% ------------------------------------------------------------
% Build true CH mobility PARAM.M and L mobility PARAM.M_L
% ------------------------------------------------------------
for ie = 1:Ne

    Mgrid = zeros(ny,nx);

    for iph = 1:Nphase
        Mgrid = Mgrid + M_phase_elem(iph,ie).*p_phase(:,:,iph);
    end

    mask = psum > eps;

    Mtmp = Mgrid;
    Mtmp(mask)  = Mgrid(mask)./psum(mask);
    Mtmp(~mask) = mean(M_phase_elem(:,ie));

    % True diffusion mobility
    PARAM.M{ie,ie}    = Mtmp;
    PARAM.M_diag{ie}  = Mtmp;

    % Mobility used only for L
    Mtmp_L  = Mtmp;
    M_floor = PHYS.M0*M_L_floor_fac(ie);

    if M_floor > 0

        if M_L_interface_only == 1

            Mtmp_L(interface_mask) = max(Mtmp_L(interface_mask),M_floor);

        else

            Mtmp_L = max(Mtmp_L,M_floor);

        end
    end

    PARAM.M_L{ie,ie}    = Mtmp_L;
    PARAM.M_L_diag{ie}  = Mtmp_L;

end

% ------------------------------------------------------------
% Compute Allen-Cahn mobility from PARAM.M_L
% ------------------------------------------------------------
zeta = zeros(ny,nx);

for ie = 1:Ne

    Mie = PARAM.M_L{ie,ie};
    Mie = max(Mie,M_min);

    zeta = zeta + dceq(ie)^2 ./ Mie;

end

zeta = max(zeta,z_min);

PARAM.L  = L_fac * 4*PHYS.m ./ (3*kap*zeta);
PARAM.Lm = PARAM.L * PHYS.m;
PARAM.LK = PARAM.L * kap;

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
PARAM.M_phase_elem     = M_phase_elem;
PARAM.p_phase          = p_phase;
PARAM.M_L_floor_fac    = M_L_floor_fac;
PARAM.M_L_interface    = interface_mask;
PARAM.M_L_interface_only = M_L_interface_only;

end