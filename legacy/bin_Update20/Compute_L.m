function PARAM = Compute_L(STATE,PARAM,PHYS)

%COMPUTE_L Stable diffusion-controlled Allen-Cahn mobility.
%
% Uses:
%
%   L = L_fac * 4*m/(3*kap*zeta)
%
% with
%
%   zeta = sum_i dceq_i^2 / M_i
%
% This allows different elemental mobilities.
%
% Important:
%   PARAM.M   is still used by CH diffusion.
%   PARAM.M_L is optional and is used only for computing L.
%
% If PARAM.M_L is not given, PARAM.M is used.
%
% Recommended order in the main loop:
%
%   PARAM = Compute_L(STATE_OLD,PARAM,PHYS);
%   PARAM = Calc_AC_Anisotropy_FacetedStiffness(STATE_OLD,PARAM,MODEL,GRID);
%
% Then:
%
%   PARAM.LK_AC = PARAM.LK .* anisotropy_factor
%   PARAM.L_AC  = PARAM.L  .* anisotropy_factor
%   PARAM.Lm_AC = PARAM.Lm .* anisotropy_factor

Ne = numel(STATE.E);
[ny,nx] = size(STATE.E{1});

% ------------------------------------------------------------
% Parameters
% ------------------------------------------------------------
L_fac = 1.0;
dceq  = 0.1;
M_min = 1e-30;
z_min = 1e-30;

if isfield(PARAM,'L_fac'), PARAM.L_fac = PARAM.L_fac; end

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

% dceq can be scalar or one value per element
if isscalar(dceq)
    dceq = dceq*ones(1,Ne);
end

if numel(dceq) ~= Ne
    error('Compute_L: PARAM.dceq must be scalar or length Ne.')
end

% ------------------------------------------------------------
% Mobility used for L
% ------------------------------------------------------------
% PARAM.M is the real CH mobility.
% PARAM.M_L, if given, is the effective mobility used only for AC mobility.
if isfield(PARAM,'M_L') && ~isempty(PARAM.M_L)
    M_for_L = PARAM.M_L;
else
    M_for_L = PARAM.M;
end

% ------------------------------------------------------------
% Diffusion resistance
% ------------------------------------------------------------
zeta = zeros(ny,nx);

for ie = 1:Ne

    Mie = Local_Read_Element_Mobility(M_for_L,ie,ny,nx);
    Mie = max(Mie,M_min);

    zeta = zeta + dceq(ie)^2 ./ Mie;

end

zeta = max(zeta,z_min);

% ------------------------------------------------------------
% Allen-Cahn mobility
% ------------------------------------------------------------
PARAM.L  = L_fac * 4*PHYS.m ./ (3*kap*zeta);
PARAM.Lm = PARAM.L * PHYS.m;
PARAM.LK = PARAM.L * kap;

end


%% ========================================================================
%  Local helper
% ========================================================================
function Mie = Local_Read_Element_Mobility(Mcell,ie,ny,nx)

%LOCAL_READ_ELEMENT_MOBILITY Read elemental mobility for Compute_L.
%
% Supported formats:
%
% Old diagonal format:
%   PARAM.M{ie}
%
% Full cell matrix format:
%   PARAM.M{ie,je}
%
% Numeric format:
%   scalar
%   ny-by-nx
%   ny-by-nx-by-Ne
%
% For full mobility, only the diagonal M{ie,ie} is used here.

if iscell(Mcell)

    if size(Mcell,1) >= ie && size(Mcell,2) >= ie && size(Mcell,1) > 1
        Mie = Mcell{ie,ie};
    else
        Mie = Mcell{ie};
    end

    if isempty(Mie)
        Mie = zeros(ny,nx);
    end

elseif isnumeric(Mcell)

    if isscalar(Mcell)

        Mie = Mcell*ones(ny,nx);

    elseif isequal(size(Mcell),[ny,nx])

        Mie = Mcell;

    elseif ndims(Mcell) == 3 && size(Mcell,3) >= ie

        Mie = Mcell(:,:,ie);

    else

        error('Compute_L: unsupported numeric mobility size.')

    end

else

    error('Compute_L: unsupported mobility format.')

end

if isscalar(Mie)
    Mie = Mie*ones(ny,nx);
end

end