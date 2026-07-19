function PARAM = Compute_M_And_L(STATE,PARAM,MODEL,PHYS)
%COMPUTE_M_AND_L Build CH mobility and pairwise AC coefficients.
%
% Pairwise interface information is read by Read_PFM_Metadata:
%
%   PHYS.sigma_pair(iph,jph)
%   PHYS.m_pair(iph,jph)
%   PHYS.kap_pair(iph,jph)
%
% The phase-pair matrix is expanded to a grain-pair matrix using
% MODEL.phase_index. Therefore:
%
%   - two different grains of the same phase use the diagonal phase-pair
%     entry sigma_pair(iph,iph);
%   - one grain never interacts with itself, so the grain-pair diagonal is 0.
%
% For each grain alpha and grid point, the local effective pair coefficient
% is weighted by the other locally present grains:
%
%   sigma_eff_alpha =
%       sum_{beta~=alpha} p_beta*sigma_{alpha,beta}
%       / sum_{beta~=alpha} p_beta
%
% In a pure two-phase interface this selects exactly that pair. At a triple
% junction it blends all locally present pairs. Using p_beta rather than
% p_alpha*p_beta gives the same normalized result when p_alpha>0 and also
% assigns the correct pair coefficient to a very small or newly activated
% alpha tail.
%
% The actual AC coefficients are:
%
%   PARAM.L_AC(:,:,alpha)
%   PARAM.Lm_AC(:,:,alpha) = L*m_eff_alpha
%   PARAM.LK_AC(:,:,alpha) = L*kap_eff_alpha
%
% PARAM.Lm and PARAM.LK are only pair-weighted 2-D diagnostic/backward-
% compatibility maps. The AC solver should use Lm_AC and LK_AC.
%
% MODEL.dgdphi contains the scalar reference PHYS.m. PARAM.m_AC_base is kept
% equal to that reference so Calc_S_AllenCahn divides it out exactly before
% multiplying by the pairwise PARAM.Lm_AC.

[ny,nx,Ngrain] = size(STATE.p);

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

if ~isfield(MODEL,'phase_index') || numel(MODEL.phase_index) ~= Ngrain
    error(['Compute_M_And_L: MODEL.phase_index must have one entry ', ...
           'per grain.'])
end

phase_index = MODEL.phase_index(:).';

if any(phase_index < 1) || any(phase_index > Nphase) || ...
        any(phase_index ~= round(phase_index))
    error('Compute_M_And_L: MODEL.phase_index contains invalid phase IDs.')
end

if ~isfield(PHYS,'M_phs')
    error('Compute_M_And_L: PHYS.M_phs is missing.')
end

required_pair = {'sigma_pair','m_pair','kap_pair'};

for i = 1:numel(required_pair)

    name = required_pair{i};

    if ~isfield(PHYS,name) || isempty(PHYS.(name))
        error('Compute_M_And_L: PHYS.%s is missing.',name)
    end

    if ~isequal(size(PHYS.(name)),[Nphase Nphase])
        error('Compute_M_And_L: PHYS.%s must be Nphase x Nphase.',name)
    end

end

if ~isfield(PHYS,'m') || ~isscalar(PHYS.m) || ...
        ~isfinite(PHYS.m) || PHYS.m <= 0
    error('Compute_M_And_L: PHYS.m reference must be positive.')
end

if ~isfield(PHYS,'kap') || ~isscalar(PHYS.kap) || ...
        ~isfinite(PHYS.kap) || PHYS.kap <= 0
    error('Compute_M_And_L: PHYS.kap reference must be positive.')
end

% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------
M_weight_mode       = 'quintic';
L_fac               = 1.0;
dceq                = 0.2;
M_min               = 1e-30;
z_min               = 1e-30;
M_L_floor_fac       = 0;
M_L_interface_only  = 1;
M_L_p_cut           = 1e-8;
pair_p_cut          = 1e-12;

if isfield(PARAM,'M_weight_mode') && ~isempty(PARAM.M_weight_mode)
    M_weight_mode = PARAM.M_weight_mode;
else
    PARAM.M_weight_mode = M_weight_mode;
end

if isfield(PARAM,'L_fac'),              L_fac              = PARAM.L_fac;              end
if isfield(PHYS,'dceq'),                dceq               = PHYS.dceq;                end
if isfield(PARAM,'dceq'),               dceq               = PARAM.dceq;               end
if isfield(PARAM,'M_min'),              M_min              = PARAM.M_min;              end
if isfield(PARAM,'L_z_min'),            z_min              = PARAM.L_z_min;            end
if isfield(PARAM,'M_L_floor_fac'),      M_L_floor_fac      = PARAM.M_L_floor_fac;      end
if isfield(PHYS,'M_L_floor_fac'),       M_L_floor_fac      = PHYS.M_L_floor_fac;       end
if isfield(PARAM,'M_L_interface_only'), M_L_interface_only = PARAM.M_L_interface_only; end
if isfield(PARAM,'M_L_p_cut'),          M_L_p_cut          = PARAM.M_L_p_cut;          end
if isfield(PARAM,'interface_pair_p_cut')
    pair_p_cut = PARAM.interface_pair_p_cut;
end

pair_p_cut = max(pair_p_cut,0);

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

% -------------------------------------------------------------------------
% Interpret PHYS.M_phs
% -------------------------------------------------------------------------
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
        error(['Compute_M_And_L: vector PHYS.M_phs must have length ', ...
               'Nphase or Ne.'])
    end

else

    [n1,n2] = size(Mraw);

    if n1 == Nphase && n2 == Ne
        M_phase_elem = Mraw;
    elseif n1 == Ne && n2 == Nphase
        M_phase_elem = Mraw.';
    else
        error(['Compute_M_And_L: PHYS.M_phs must be scalar, ', ...
               'Nphase x Ne, Ne x Nphase, length Nphase, or length Ne.'])
    end

end

M_phase_elem = max(M_phase_elem,M_min);

% -------------------------------------------------------------------------
% Collapse grain p to thermodynamic-phase p
% -------------------------------------------------------------------------
p_phase = zeros(ny,nx,Nphase);

for iph = 1:Nphase

    grains = find(phase_index == iph);

    if ~isempty(grains)
        p_phase(:,:,iph) = sum(STATE.p(:,:,grains),3);
    end

end

% -------------------------------------------------------------------------
% Mobility interpolation weights
% -------------------------------------------------------------------------
switch lower(char(M_weight_mode))

    case 'linear'
        p_mob = p_phase;

    case 'smoothstep'
        p_mob = p_phase.^2 .* (3 - 2*p_phase);

    case 'quintic'
        p_mob = p_phase.^3 .* ...
            (10 - 15*p_phase + 6*p_phase.^2);

    otherwise
        error('Compute_M_And_L: unknown PARAM.M_weight_mode "%s".', ...
              char(M_weight_mode))
end

psum_mob = sum(p_mob,3);

% This mask affects only the optional mobility floor used to compute L.
interface_mask = (1 - sum(p_phase.^2,3)) > M_L_p_cut;

% -------------------------------------------------------------------------
% Allocate and build true CH mobility
% -------------------------------------------------------------------------
PARAM.M = cell(Ne,Ne);

for ie = 1:Ne
    for je = 1:Ne
        PARAM.M{ie,je} = zeros(ny,nx);
    end
end

PARAM.M_diag = cell(1,Ne);
zeta        = zeros(ny,nx);

for ie = 1:Ne

    Mgrid = zeros(ny,nx);

    for iph = 1:Nphase
        Mgrid = Mgrid + M_phase_elem(iph,ie).*p_mob(:,:,iph);
    end

    mask = psum_mob > eps;

    Mtmp = Mgrid;
    Mtmp(mask)  = Mgrid(mask)./psum_mob(mask);
    Mtmp(~mask) = mean(M_phase_elem(:,ie));
    Mtmp        = max(Mtmp,M_min);

    % True CH diffusion mobility.
    PARAM.M{ie,ie}   = Mtmp;
    PARAM.M_diag{ie} = Mtmp;

    % Mobility used only to calculate the AC kinetic L.
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
            Mtmp_L(interface_mask) = ...
                max(Mtmp_L(interface_mask),M_floor);
        else
            Mtmp_L = max(Mtmp_L,M_floor);
        end

    end

    Mtmp_L = max(Mtmp_L,M_min);
    zeta   = zeta + dceq(ie)^2 ./ Mtmp_L;

end

zeta = max(zeta,z_min);

% -------------------------------------------------------------------------
% Common Allen-Cahn kinetic mobility
%
% Since every pair uses the same interface thickness l:
%   m_pair/kap_pair = 8/l^2
% and L is independent of sigma_pair. The scalar reference ratio therefore
% gives the same L as every individual pair.
% -------------------------------------------------------------------------
PARAM.L = L_fac * 4*PHYS.m ./ (3*PHYS.kap.*zeta);

if any(~isfinite(PARAM.L(:))) || any(PARAM.L(:) < 0)
    error('Compute_M_And_L: invalid PARAM.L was generated.')
end

% -------------------------------------------------------------------------
% Expand phase-pair coefficients to grain-pair coefficients
% -------------------------------------------------------------------------
sigma_pair_grain = PHYS.sigma_pair(phase_index,phase_index);
m_pair_grain     = PHYS.m_pair(phase_index,phase_index);
kap_pair_grain   = PHYS.kap_pair(phase_index,phase_index);

% A grain has no interface with itself.
sigma_pair_grain(1:Ngrain+1:end) = 0;
m_pair_grain(1:Ngrain+1:end)     = 0;
kap_pair_grain(1:Ngrain+1:end)   = 0;

if isfield(PHYS,'sigma_pair_SI')
    PARAM.sigma_pair_phase_SI = PHYS.sigma_pair_SI;
end

PARAM.sigma_pair_phase = PHYS.sigma_pair;
PARAM.m_pair_phase     = PHYS.m_pair;
PARAM.kap_pair_phase   = PHYS.kap_pair;
PARAM.L_pair           = PARAM.L;
PARAM.L_pair_is_common = true;

PARAM.sigma_pair_grain = sigma_pair_grain;
PARAM.m_pair_grain     = m_pair_grain;
PARAM.kap_pair_grain   = kap_pair_grain;

% -------------------------------------------------------------------------
% Local grain-resolved effective pair coefficients
% -------------------------------------------------------------------------
PARAM.sigma_eff_AC = zeros(ny,nx,Ngrain);
PARAM.m_eff_AC     = zeros(ny,nx,Ngrain);
PARAM.kap_eff_AC   = zeros(ny,nx,Ngrain);

PARAM.L_AC_base  = zeros(ny,nx,Ngrain);
PARAM.Lm_AC_base = zeros(ny,nx,Ngrain);
PARAM.LK_AC_base = zeros(ny,nx,Ngrain);

for alpha = 1:Ngrain

    sigma_num = zeros(ny,nx);
    m_num     = zeros(ny,nx);
    kap_num   = zeros(ny,nx);
    wsum      = zeros(ny,nx);

    for beta = 1:Ngrain

        if beta == alpha
            continue
        end

        % Selecting by p_beta gives the exact alpha-beta value in a pure
        % binary interface and blends partners at a triple junction.
        w = STATE.p(:,:,beta);

        sigma_num = sigma_num + w.*sigma_pair_grain(alpha,beta);
        m_num     = m_num     + w.*m_pair_grain(alpha,beta);
        kap_num   = kap_num   + w.*kap_pair_grain(alpha,beta);
        wsum      = wsum      + w;

    end

    % Safe fallback in a one-grain pure interior. Interface contributions are
    % zero there because W'(0/1)=0 and grad(phi)=0.
    sigma_eff = PHYS.sigma_ref*ones(ny,nx);
    m_eff     = PHYS.m        *ones(ny,nx);
    kap_eff   = PHYS.kap      *ones(ny,nx);

    mask_pair = wsum > pair_p_cut;

    sigma_eff(mask_pair) = sigma_num(mask_pair)./wsum(mask_pair);
    m_eff(mask_pair)     = m_num(mask_pair)    ./wsum(mask_pair);
    kap_eff(mask_pair)   = kap_num(mask_pair)  ./wsum(mask_pair);

    if any(~isfinite(sigma_eff(:))) || any(sigma_eff(:) <= 0) || ...
       any(~isfinite(m_eff(:)))     || any(m_eff(:) <= 0)     || ...
       any(~isfinite(kap_eff(:)))   || any(kap_eff(:) <= 0)

        error(['Compute_M_And_L: invalid local pair coefficient for ', ...
               'grain %d.'],alpha)
    end

    PARAM.sigma_eff_AC(:,:,alpha) = sigma_eff;
    PARAM.m_eff_AC(:,:,alpha)     = m_eff;
    PARAM.kap_eff_AC(:,:,alpha)   = kap_eff;

    PARAM.L_AC_base(:,:,alpha)  = PARAM.L;
    PARAM.Lm_AC_base(:,:,alpha) = PARAM.L.*m_eff;
    PARAM.LK_AC_base(:,:,alpha) = PARAM.L.*kap_eff;

end

% Active coefficients consumed by Calc_S_AllenCahn and the ACCH solver.
% An anisotropy routine may overwrite these from the immutable *_AC_base
% arrays after Compute_M_And_L.
PARAM.L_AC  = PARAM.L_AC_base;
PARAM.Lm_AC = PARAM.Lm_AC_base;
PARAM.LK_AC = PARAM.LK_AC_base;

% Reference used only to remove m from MODEL.dgdphi.
PARAM.m_AC_base    = PHYS.m;
PARAM.m_dgdphi_ref = PHYS.m;

% -------------------------------------------------------------------------
% One common 2-D pair map for diagnostics and backward compatibility
%
% This map is pair weighted over every alpha<beta. It cannot replace the
% grain-resolved Lm_AC/LK_AC at a triple junction.
% -------------------------------------------------------------------------
sigma_num_global = zeros(ny,nx);
m_num_global     = zeros(ny,nx);
kap_num_global   = zeros(ny,nx);
pair_wsum        = zeros(ny,nx);

for alpha = 1:Ngrain-1
    for beta = alpha+1:Ngrain

        w = STATE.p(:,:,alpha).*STATE.p(:,:,beta);

        sigma_num_global = sigma_num_global + ...
            w.*sigma_pair_grain(alpha,beta);
        m_num_global = m_num_global + ...
            w.*m_pair_grain(alpha,beta);
        kap_num_global = kap_num_global + ...
            w.*kap_pair_grain(alpha,beta);
        pair_wsum = pair_wsum + w;

    end
end

PARAM.sigma_eff = PHYS.sigma_ref*ones(ny,nx);
PARAM.m_eff     = PHYS.m        *ones(ny,nx);
PARAM.kap_eff   = PHYS.kap      *ones(ny,nx);

mask_global = pair_wsum > pair_p_cut;

PARAM.sigma_eff(mask_global) = ...
    sigma_num_global(mask_global)./pair_wsum(mask_global);
PARAM.m_eff(mask_global) = ...
    m_num_global(mask_global)./pair_wsum(mask_global);
PARAM.kap_eff(mask_global) = ...
    kap_num_global(mask_global)./pair_wsum(mask_global);

PARAM.Lm = PARAM.L.*PARAM.m_eff;
PARAM.LK = PARAM.L.*PARAM.kap_eff;

% -------------------------------------------------------------------------
% Diagnostics
% -------------------------------------------------------------------------
PARAM.M_phase              = p_phase;
PARAM.M_phase_elem         = M_phase_elem;
PARAM.M_L_floor_fac        = M_L_floor_fac;
PARAM.M_L_interface_only   = M_L_interface_only;
PARAM.M_L_p_cut            = M_L_p_cut;
PARAM.interface_pair_p_cut = pair_p_cut;

PARAM.zeta_L               = zeta;
PARAM.pair_weight_sum      = pair_wsum;
PARAM.pairwise_interface_active = true;

PARAM.sigma_eff_minmax = ...
    [min(PARAM.sigma_eff,[],'all'),max(PARAM.sigma_eff,[],'all')];
PARAM.m_eff_minmax = ...
    [min(PARAM.m_eff,[],'all'),max(PARAM.m_eff,[],'all')];
PARAM.L_minmax = ...
    [min(PARAM.L,[],'all'),max(PARAM.L,[],'all')];
PARAM.Lm_minmax = ...
    [min(PARAM.Lm,[],'all'),max(PARAM.Lm,[],'all')];
PARAM.LK_minmax = ...
    [min(PARAM.LK,[],'all'),max(PARAM.LK,[],'all')];

end
