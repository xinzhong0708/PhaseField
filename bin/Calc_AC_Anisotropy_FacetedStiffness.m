function PARAM = Calc_AC_Anisotropy_FacetedStiffness(STATE,PARAM,MODEL,GRID)
%CALC_AC_ANISOTROPY_FACETEDSTIFFNESS Kundin-style faceted AC anisotropy.
%
% Metadata fields used:
%
%   PARAM.facet(iph).mode
%   PARAM.facet(iph).theta
%   PARAM.facet(iph).theta_deg
%   PARAM.facet(iph).A
%   PARAM.facet(iph).A_weight
%   PARAM.facet(iph).q
%   PARAM.facet(iph).amin
%   PARAM.facet(iph).amax
%
% Kundin-style stiffness:
%
%   sigma* = sigma_ref * (A/A_ref) * q^2 ...
%          / (sin(dtheta)^2 + q^2*cos(dtheta)^2)^(3/2)
%
% Larger A means larger facet area / stronger persistence.
%
% This version keeps the legacy behavior:
%   L_AC, Lm_AC, and LK_AC are all scaled by the anisotropic factor.

[ny,nx,Ngrain] = size(STATE.phi);

% ------------------------------------------------------------
% Defaults
% ------------------------------------------------------------
nfold       = 6;
q_default   = 0.4;
s_min       = 0.3;
s_max       = 3.0;
grad_min    = 1e-12;
p_cut       = 1e-8;
normalize_s = 1;
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

PARAM.sigma_star_fac = ones(ny,nx,Ngrain);
PARAM.aniso_theta    = zeros(ny,nx,Ngrain);
PARAM.aniso_factor   = ones(ny,nx,Ngrain);
PARAM.aniso_ifacet   = zeros(ny,nx,Ngrain);

sigma_star_raw = ones(ny,nx,Ngrain);
interface_mask = false(ny,nx,Ngrain);

for ig = 1:Ngrain
    PARAM.L_AC(:,:,ig)  = PARAM.L;
    PARAM.Lm_AC(:,:,ig) = PARAM.Lm;
    PARAM.LK_AC(:,:,ig) = PARAM.LK;
end

% ------------------------------------------------------------
% First pass: raw Kundin stiffness for each grain
% ------------------------------------------------------------
for alpha = 1:Ngrain

    phi = STATE.phi(:,:,alpha);

    [gx,gy] = Local_Grad_Reflex(phi,GRID.dx,GRID.dy);

    gnorm = sqrt(gx.^2 + gy.^2);
    mask  = gnorm > grad_min;

    theta_n = atan2(gy,gx);

    if faceted_grain(alpha)

        iph = MODEL.phase_index(alpha);

        [facet_theta,facet_A,sigma_ref,A_ref,q,smin_phase,smax_phase] = ...
            Local_Get_Facet_Data(PARAM,iph,nfold,q_default,s_min,s_max);

        % Legacy behavior:
        % rotate crystal facet normals into lab frame.
        facet_theta = facet_theta + theta_grain(alpha);

        [dtheta,ifacet] = Local_Nearest_Facet_Angle(theta_n,facet_theta);

        A_use = facet_A(ifacet);

        den = sin(dtheta).^2 + q^2*cos(dtheta).^2;

        sigma_star = sigma_ref .* (A_use./A_ref) .* q^2 ./ den.^(3/2);

        sigma_star(~mask) = sigma_ref;
        sigma_star = min(max(sigma_star,smin_phase),smax_phase);

    else

        dtheta     = zeros(ny,nx);
        ifacet     = zeros(ny,nx);
        sigma_star = ones(ny,nx);

    end

    sigma_star_raw(:,:,alpha) = sigma_star;
    interface_mask(:,:,alpha) = mask;

    PARAM.sigma_star_fac(:,:,alpha) = sigma_star;
    PARAM.aniso_theta(:,:,alpha)    = dtheta;
    PARAM.aniso_ifacet(:,:,alpha)   = ifacet;

end

% ------------------------------------------------------------
% Second pass: pairwise effective stiffness
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

        % Diffuse-interface weight.
        w = p_a .* p_b;

        % Important for sharp 0/1 initial maps.
        % If p_a*p_b is zero, use overlapping gradient masks.
        if max(w(:)) <= p_cut
            w = double(interface_mask(:,:,alpha) & interface_mask(:,:,beta));
        end

        if max(w(:)) <= p_cut
            continue
        end

        iph_a = MODEL.phase_index(alpha);
        iph_b = MODEL.phase_index(beta);

        sig_a = sigma_star_raw(:,:,alpha);
        sig_b = sigma_star_raw(:,:,beta);

        if iph_a == iph_b

            % Same thermodynamic phase grain boundary.
            % Default: isotropic unless explicitly enabled.
            s_ab = ones(ny,nx);

            if isfield(PARAM,'aniso_same_phase') && PARAM.aniso_same_phase == 1
                s_ab = 0.5*rsl*(sig_a + sig_b);
            end

        else

            % Different phases.
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

    % Legacy behavior: normalize active interface.
    % This keeps the mean interface strength comparable to isotropic.
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

    % Legacy behavior: scale all three.
    PARAM.L_AC(:,:,alpha)  = PARAM.L  .* s_eff;
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


function [facet_theta,facet_A,sigma_ref,A_ref,q,s_min,s_max] = ...
    Local_Get_Facet_Data(PARAM,iph,nfold,q_default,s_min_default,s_max_default)

% Default equal facets over 2*pi, as in the legacy function.
facet_theta = 2*pi*(0:nfold-1)/nfold;
facet_A     = ones(size(facet_theta));
sigma_ref   = 1.0;
A_ref       = 1.0;
q           = q_default;
s_min       = s_min_default;
s_max       = s_max_default;

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

elseif ~isempty(facet_A)

    % Legacy behavior: use the first listed facet as reference.
    % If you want Kundin's exact {001} reference, put {001} first or set A_ref.
    A_ref = facet_A(1);

end

if isfield(facet,'q') && ~isempty(facet.q)
    q = facet.q;
end

if isfield(facet,'amin') && ~isempty(facet.amin)
    s_min = facet.amin;
end

if isfield(facet,'amax') && ~isempty(facet.amax)
    s_max = facet.amax;
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


function [dtheta_best,ifacet_best] = Local_Nearest_Facet_Angle(theta,facet_theta)

[ny,nx] = size(theta);
Nfacet  = numel(facet_theta);

abs_best    = inf(ny,nx);
dtheta_best = zeros(ny,nx);
ifacet_best = ones(ny,nx);

for k = 1:Nfacet

    dtheta = Local_WrapPi(theta - facet_theta(k));
    take   = abs(dtheta) < abs_best;

    dtheta_best(take) = dtheta(take);
    abs_best(take)    = abs(dtheta(take));
    ifacet_best(take) = k;

end

end


function a = Local_WrapPi(a)

a = mod(a + pi,2*pi) - pi;

end