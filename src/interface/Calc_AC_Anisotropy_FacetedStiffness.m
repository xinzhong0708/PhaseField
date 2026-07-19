function PARAM = Calc_AC_Anisotropy_FacetedStiffness(STATE,PARAM,MODEL,GRID)
%CALC_AC_ANISOTROPY_FACETEDSTIFFNESS Kundin-style faceted AC anisotropy.
%
% This version is for the current scalar AC solver.
%
% It uses Kundin-type surface energy as the scalar multiplier:
%
%   sigma/sigma_ref = (A_ref/A_facet)
%       * sqrt(sin(dtheta)^2 + q^2*cos(dtheta)^2) / q
%
% The true Kundin stiffness:
%
%   sigma* = sigma_ref * (A_facet/A_ref) * q^2
%       / (sin(dtheta)^2 + q^2*cos(dtheta)^2)^(3/2)
%
% is stored as PARAM.sigma_star_fac only for diagnostics.
%
% Larger A means larger target facet area / more persistent facet.
%
% Outputs:
%   PARAM.L_AC(:,:,ig)
%   PARAM.Lm_AC(:,:,ig)
%   PARAM.LK_AC(:,:,ig)

[ny,nx,Ngrain] = size(STATE.phi);

% ------------------------------------------------------------
% Defaults
% ------------------------------------------------------------
nfold       = 4;
q_default   = 0.25;
s_min       = 0.1;
s_max       = 20.0;
grad_min    = 1e-12;
p_cut       = 1e-8;
normalize_s = 0;
rsl         = 1.0;

if isfield(PARAM,'aniso_nfold'),       nfold       = PARAM.aniso_nfold;       end
if isfield(PARAM,'aniso_q'),           q_default   = PARAM.aniso_q;           end
if isfield(PARAM,'aniso_min'),         s_min       = PARAM.aniso_min;         end
if isfield(PARAM,'aniso_max'),         s_max       = PARAM.aniso_max;         end
if isfield(PARAM,'aniso_grad_min'),    grad_min    = PARAM.aniso_grad_min;    end
if isfield(PARAM,'aniso_p_cut'),       p_cut       = PARAM.aniso_p_cut;       end
if isfield(PARAM,'aniso_normalize'),   normalize_s = PARAM.aniso_normalize;   end
if isfield(PARAM,'aniso_rsl'),         rsl         = PARAM.aniso_rsl;         end

if ~isfield(PARAM,'L') || ~isfield(PARAM,'Lm') || ~isfield(PARAM,'LK')
    error('Calc_AC_Anisotropy_FacetedStiffness: PARAM.L, PARAM.Lm, PARAM.LK must exist.')
end

theta_grain = zeros(1,Ngrain);

if isfield(PARAM,'theta_grain')
    theta_grain(1:min(numel(PARAM.theta_grain),Ngrain)) = ...
        PARAM.theta_grain(1:min(numel(PARAM.theta_grain),Ngrain));
end

% ------------------------------------------------------------
% Which grains are faceted?
% ------------------------------------------------------------
faceted_grain = false(1,Ngrain);

if isfield(PARAM,'aniso_phase') && ~isempty(PARAM.aniso_phase)

    for ig = 1:Ngrain
        iph = MODEL.phase_index(ig);

        if any(iph == PARAM.aniso_phase)
            faceted_grain(ig) = true;
        end
    end

elseif isfield(PARAM,'facet') && ~isempty(PARAM.facet)

    for ig = 1:Ngrain

        iph = MODEL.phase_index(ig);

        if iph <= numel(PARAM.facet)

            if isfield(PARAM.facet(iph),'mode') && ...
                    ~isempty(PARAM.facet(iph).mode) && ...
                    ~strcmpi(PARAM.facet(iph).mode,'iso')

                faceted_grain(ig) = true;

            end
        end
    end
end

if isfield(PARAM,'aniso_grains') && ~isempty(PARAM.aniso_grains)
    faceted_grain(:) = false;
    faceted_grain(PARAM.aniso_grains) = true;
end

% ------------------------------------------------------------
% Allocate output
% ------------------------------------------------------------
PARAM.L_AC  = zeros(ny,nx,Ngrain);
PARAM.Lm_AC = zeros(ny,nx,Ngrain);
PARAM.LK_AC = zeros(ny,nx,Ngrain);

PARAM.sigma_energy_fac = ones(ny,nx,Ngrain);
PARAM.sigma_star_fac   = ones(ny,nx,Ngrain);
PARAM.aniso_theta      = zeros(ny,nx,Ngrain);
PARAM.aniso_factor     = ones(ny,nx,Ngrain);
PARAM.aniso_ifacet     = zeros(ny,nx,Ngrain);

sigma_fac_raw  = ones(ny,nx,Ngrain);
sigma_star_raw = ones(ny,nx,Ngrain);
interface_mask = false(ny,nx,Ngrain);

for ig = 1:Ngrain
    PARAM.L_AC(:,:,ig)  = PARAM.L;
    PARAM.Lm_AC(:,:,ig) = PARAM.Lm;
    PARAM.LK_AC(:,:,ig) = PARAM.LK;
end

% ------------------------------------------------------------
% First pass: grain-wise surface-energy factor
% ------------------------------------------------------------
for alpha = 1:Ngrain

    phi = STATE.phi(:,:,alpha);

    [gx,gy] = Local_Grad_Reflex(phi,GRID.dx,GRID.dy);

    gnorm = sqrt(gx.^2 + gy.^2);
    mask  = gnorm > grad_min;

    theta_lab = atan2(gy,gx);

    if faceted_grain(alpha)

        iph = MODEL.phase_index(alpha);

        [facet_theta,facet_A,sigma_ref,A_ref,q] = ...
            Local_Get_Facet_Data(PARAM,iph,nfold,q_default);

        % Convert lab normal angle into crystal frame.
        theta_crystal = Local_WrapHalfPi(theta_lab - theta_grain(alpha));

        [dtheta,ifacet] = Local_Nearest_Facet_Angle(theta_crystal,facet_theta);

        A_use = facet_A(ifacet);

        den = sin(dtheta).^2 + q^2*cos(dtheta).^2;

        % Kundin surface-energy factor.
        % At facet normal:
        %   sigma_fac = A_ref/A_use
        %
        % Larger A_use gives lower surface energy and longer facet.
        sigma_energy = sigma_ref .* (A_ref./A_use) .* sqrt(den);
        sigma_fac    = sigma_energy ./ max(sigma_ref*q,eps);

        % True Kundin stiffness, diagnostic only.
        sigma_star = sigma_ref .* (A_use./A_ref) .* q^2 ./ den.^(3/2);

        sigma_fac(~mask)  = 1;
        sigma_star(~mask) = 1;

    else

        dtheta      = zeros(ny,nx);
        ifacet      = zeros(ny,nx);
        sigma_fac   = ones(ny,nx);
        sigma_star  = ones(ny,nx);

    end

    sigma_fac_raw(:,:,alpha)  = sigma_fac;
    sigma_star_raw(:,:,alpha) = sigma_star;
    interface_mask(:,:,alpha) = mask;

    PARAM.sigma_energy_fac(:,:,alpha) = sigma_fac;
    PARAM.sigma_star_fac(:,:,alpha)   = sigma_star;
    PARAM.aniso_theta(:,:,alpha)      = dtheta;
    PARAM.aniso_ifacet(:,:,alpha)     = ifacet;

end

% ------------------------------------------------------------
% Second pass: pair effective surface-energy factor
% ------------------------------------------------------------
for alpha = 1:Ngrain

    p_a   = STATE.p(:,:,alpha);
    s_sum = zeros(ny,nx);
    w_sum = zeros(ny,nx);

    for beta = 1:Ngrain

        if beta == alpha
            continue
        end

        p_b = STATE.p(:,:,beta);

        % Use p_a*p_b in diffuse interfaces.
        w = p_a .* p_b;

        % For sharp 0/1 initial maps, p_a*p_b may be zero.
        if max(w(:)) <= p_cut
            w = double(interface_mask(:,:,alpha) & interface_mask(:,:,beta));
        end

        if max(w(:)) <= p_cut
            continue
        end

        iph_a = MODEL.phase_index(alpha);
        iph_b = MODEL.phase_index(beta);

        sig_a = sigma_fac_raw(:,:,alpha);
        sig_b = sigma_fac_raw(:,:,beta);

        if iph_a == iph_b

            s_ab = ones(ny,nx);

            if isfield(PARAM,'aniso_same_phase') && PARAM.aniso_same_phase == 1
                s_ab = 0.5*rsl*(sig_a + sig_b);
            end

        else

            if faceted_grain(alpha) && ~faceted_grain(beta)

                s_ab = sig_a;

            elseif ~faceted_grain(alpha) && faceted_grain(beta)

                s_ab = sig_b;

            elseif faceted_grain(alpha) && faceted_grain(beta)

                s_ab = 0.5*rsl*(sig_a + sig_b);

            else

                s_ab = ones(ny,nx);

            end
        end

        s_sum = s_sum + w .* s_ab;
        w_sum = w_sum + w;

    end

    s_eff = ones(ny,nx);
    mask_pair = w_sum > p_cut;

    s_eff(mask_pair) = s_sum(mask_pair)./w_sum(mask_pair);

    if normalize_s == 1

        mask_norm = mask_pair & isfinite(s_eff);

        if any(mask_norm(:))
            s_mean = mean(s_eff(mask_norm));

            if isfinite(s_mean) && s_mean > 0
                s_eff = s_eff./s_mean;
            end
        end
    end

    s_eff = min(max(s_eff,s_min),s_max);
    s_eff(~mask_pair) = 1;

    PARAM.aniso_factor(:,:,alpha) = s_eff;

    % Apply anisotropic surface-energy factor to interface terms only.
    % Do not scale L_AC here.
    PARAM.L_AC(:,:,alpha)  = PARAM.L;
    PARAM.Lm_AC(:,:,alpha) = PARAM.Lm .* s_eff;
    PARAM.LK_AC(:,:,alpha) = PARAM.LK .* s_eff;

end

end


%% ========================================================================
%  Local helper functions
% ========================================================================

function [gx,gy] = Local_Grad_Reflex(A,dx,dy)

[ny,nx] = size(A);

if nx == 1
    AL = A;
    AR = A;
else
    AL = A(:,[2,1:nx-1]);
    AR = A(:,[2:nx,nx-1]);
end

if ny == 1
    AU = A;
    AD = A;
else
    AU = A([2,1:ny-1],:);
    AD = A([2:ny,ny-1],:);
end

gx = (AR - AL)/(2*dx);
gy = (AD - AU)/(2*dy);

end


function [facet_theta,facet_A,sigma_ref,A_ref,q] = ...
    Local_Get_Facet_Data(PARAM,iph,nfold,q_default)

facet_theta = pi*(0:nfold-1)/nfold;
facet_A     = ones(size(facet_theta));
sigma_ref   = 1.0;
A_ref       = 1.0;
q           = q_default;

if ~isfield(PARAM,'facet') || numel(PARAM.facet) < iph
    return
end

facet = PARAM.facet(iph);

if isfield(facet,'theta') && ~isempty(facet.theta)

    facet_theta = facet.theta(:).';

elseif isfield(facet,'theta_deg') && ~isempty(facet.theta_deg)

    facet_theta = facet.theta_deg(:).'*pi/180;

end

if isfield(facet,'A') && ~isempty(facet.A)

    facet_A = facet.A(:).';

elseif isfield(facet,'A_weight') && ~isempty(facet.A_weight)

    facet_A = facet.A_weight(:).';

else

    facet_A = ones(size(facet_theta));

end

if isfield(facet,'sigma_ref') && ~isempty(facet.sigma_ref)
    sigma_ref = facet.sigma_ref;
end

if isfield(facet,'A_ref') && ~isempty(facet.A_ref)

    A_ref = facet.A_ref;

elseif isfield(facet,'hkl') && ~isempty(facet.hkl)

    i001 = Local_Find_HKL(facet.hkl,'001');

    if ~isempty(i001)
        A_ref = facet_A(i001);
    else
        A_ref = facet_A(1);
    end

elseif ~isempty(facet_A)

    A_ref = facet_A(1);

end

if isfield(facet,'q') && ~isempty(facet.q)
    q = facet.q;
end

if numel(facet_A) ~= numel(facet_theta)
    error('facet_A and facet_theta must have the same length.')
end

if any(facet_A <= 0)
    error('All facet A values must be positive.')
end

if A_ref <= 0
    error('A_ref must be positive.')
end

end


function idx = Local_Find_HKL(hkl,target)

idx = [];

for i = 1:numel(hkl)

    s = hkl{i};

    if isstring(s)
        s = char(s);
    end

    s = regexprep(s,'[{}()\[\]\s]','');

    if strcmpi(s,target)
        idx = i;
        return
    end

end

end


function [dtheta_best,ifacet_best] = Local_Nearest_Facet_Angle(theta,facet_theta)

[ny,nx] = size(theta);
Nfacet  = numel(facet_theta);

abs_best    = inf(ny,nx);
dtheta_best = zeros(ny,nx);
ifacet_best = ones(ny,nx);

for k = 1:Nfacet

    dtheta = Local_WrapHalfPi(theta - facet_theta(k));
    take   = abs(dtheta) < abs_best;

    dtheta_best(take) = dtheta(take);
    abs_best(take)    = abs(dtheta(take));
    ifacet_best(take) = k;

end

end


function a = Local_WrapHalfPi(a)

a = mod(a + pi/2,pi) - pi/2;

end