function PARAM = Calc_AC_Anisotropy_FacetedStiffness(STATE,PARAM,MODEL,GRID)
%CALC_AC_ANISOTROPY_FACETEDSTIFFNESS Paper-style faceted AC anisotropy.
%
% This function computes frozen anisotropic interface coefficients from
% STATE.phi and writes:
%
%   PARAM.LK_AC(:,:,alpha)
%   PARAM.L_AC(:,:,alpha)
%   PARAM.Lm_AC(:,:,alpha)
%
% The anisotropy is based on the paper-style surface stiffness:
%
%   sigma* = sigma_ref * (A_facet/A_ref) * q^2 ...
%          / (sin(theta)^2 + q^2*cos(theta)^2)^(3/2)
%
% where theta is the angle between the local interface normal and the
% nearest facet normal.
%
% The coefficients are frozen during one ACCH solve.
%
% Required:
%   STATE.phi
%   STATE.p
%   PARAM.L
%   PARAM.LK
%   PARAM.Lm
%
% Optional simple symmetry:
%   PARAM.aniso_nfold       = 6;
%   PARAM.aniso_q           = 0.4;
%   PARAM.theta_grain       = zeros(1,Ngrain);
%
% Optional general facet table:
%   PARAM.facet(iph).theta      = [0 pi/3 2*pi/3 ...];  % facet normals
%   PARAM.facet(iph).A          = [1 1 1 ...];          % facet areas
%   PARAM.facet(iph).sigma_ref  = 1;
%   PARAM.facet(iph).A_ref      = 1;
%   PARAM.facet(iph).q          = 0.4;
%
% Notes:
%   - facet.theta is in crystal coordinates. theta_grain rotates it.
%   - hkl can be stored in PARAM.facet(iph).hkl for reference, but this
%     function uses facet normal angles in 2D.
%   - CH, LE and PARAM.M are not modified.

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

theta_grain = zeros(1,Ngrain);

if isfield(PARAM,'theta_grain')
    theta_grain(1:min(numel(PARAM.theta_grain),Ngrain)) = ...
        PARAM.theta_grain(1:min(numel(PARAM.theta_grain),Ngrain));
end

% Which grains are faceted?
faceted_grain = true(1,Ngrain);

if isfield(PARAM,'aniso_grains') && ~isempty(PARAM.aniso_grains)
    faceted_grain      = false(1,Ngrain);
    faceted_grain(PARAM.aniso_grains) = true;
end

if isfield(PARAM,'aniso_phase') && ~isempty(PARAM.aniso_phase)
    faceted_grain = false(1,Ngrain);

    for ig = 1:Ngrain
        iph = MODEL.phase_index(ig);

        if any(iph == PARAM.aniso_phase)
            faceted_grain(ig) = true;
        end
    end
end

% ------------------------------------------------------------
% Allocate output
% ------------------------------------------------------------
PARAM.LK_AC = zeros(ny,nx,Ngrain);
PARAM.L_AC  = zeros(ny,nx,Ngrain);
PARAM.Lm_AC = zeros(ny,nx,Ngrain);

PARAM.sigma_star_fac = zeros(ny,nx,Ngrain);
PARAM.aniso_theta    = zeros(ny,nx,Ngrain);
PARAM.aniso_factor   = zeros(ny,nx,Ngrain);

sigma_star_raw = ones(ny,nx,Ngrain);
interface_mask = false(ny,nx,Ngrain);

% ------------------------------------------------------------
% First pass: grain-wise solid-liquid stiffness
% ------------------------------------------------------------
for alpha = 1:Ngrain

    phi = STATE.phi(:,:,alpha);

    [gx,gy] = Local_Grad_Reflex(phi,GRID.dx,GRID.dy);

    gnorm = sqrt(gx.^2 + gy.^2);
    mask  = gnorm > grad_min;

    theta_n = atan2(gy,gx);

    if faceted_grain(alpha)

        iph = MODEL.phase_index(alpha);

        [facet_theta,facet_A,sigma_ref,A_ref,q] = ...
            Local_Get_Facet_Data(PARAM,iph,nfold,q_default);

        % Rotate crystal facet normals by grain orientation
        facet_theta = facet_theta + theta_grain(alpha);

        [dtheta,ifacet] = Local_Nearest_Facet_Angle(theta_n,facet_theta);

        A_use = facet_A(ifacet);

        % Paper-style stiffness:
        %   sigma* = sigma_ref * A/A_ref * q^2 ...
        %          / (sin^2(theta) + q^2*cos^2(theta))^(3/2)
        den = sin(dtheta).^2 + q^2*cos(dtheta).^2;

        sigma_star = sigma_ref .* (A_use./A_ref) .* q^2 ./ den.^(3/2);

        sigma_star(~mask) = sigma_ref;

    else

        dtheta     = zeros(ny,nx);
        sigma_star = ones(ny,nx);

    end

    sigma_star_raw(:,:,alpha) = sigma_star;
    interface_mask(:,:,alpha) = mask;

    PARAM.sigma_star_fac(:,:,alpha) = sigma_star;
    PARAM.aniso_theta(:,:,alpha)    = dtheta;

end

% ------------------------------------------------------------
% Second pass: pairwise effective stiffness for each grain
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
        w   = p_a .* p_b;

        if max(w(:)) <= p_cut
            continue
        end

        iph_a = MODEL.phase_index(alpha);
        iph_b = MODEL.phase_index(beta);

        sig_a = sigma_star_raw(:,:,alpha);
        sig_b = sigma_star_raw(:,:,beta);

        if iph_a == iph_b

            % Same thermodynamic phase grain boundary.
            % Default: chemically and capillarily simple/isotropic.
            s_ab = ones(ny,nx);

            if isfield(PARAM,'aniso_same_phase') && PARAM.aniso_same_phase == 1
                s_ab = 0.5*(sig_a + sig_b);
            end

        else

            % Different phases.
            % If one side is faceted and the other is melt/isotropic, the
            % faceted side dominates the pair stiffness. If both are
            % faceted, take the mean as a simple pair approximation.
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

    s_eff(mask_pair) = s_sum(mask_pair) ./ w_sum(mask_pair);

    % Normalize over active interface of this grain.
    % This keeps the mean interface strength close to PARAM.LK/PARAM.L.
    if normalize_s == 1
        mask_norm = mask_pair & isfinite(s_eff);

        if any(mask_norm(:))
            s_mean = mean(s_eff(mask_norm));
            if isfinite(s_mean) && s_mean > 0
                s_eff = s_eff ./ s_mean;
            end
        end
    end

    % Cap for numerical stability
    s_eff = min(max(s_eff,s_min),s_max);

    % No anisotropy in pure bulk
    s_eff(~mask_pair) = 1;

    PARAM.aniso_factor(:,:,alpha) = s_eff;

    PARAM.LK_AC(:,:,alpha) = PARAM.LK .* s_eff;
    PARAM.L_AC(:,:,alpha)  = PARAM.L  .* s_eff;
    PARAM.Lm_AC(:,:,alpha) = PARAM.Lm .* s_eff;

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

%LOCAL_GET_FACET_DATA Return facet normals and facet parameters.
%
% Simple mode:
%   nfold symmetry with equal facet areas.
%
% General mode:
%   PARAM.facet(iph).theta
%   PARAM.facet(iph).A
%   PARAM.facet(iph).sigma_ref
%   PARAM.facet(iph).A_ref
%   PARAM.facet(iph).q

if isfield(PARAM,'facet') && numel(PARAM.facet) >= iph && ...
        isfield(PARAM.facet(iph),'theta') && ~isempty(PARAM.facet(iph).theta)

    facet_theta = PARAM.facet(iph).theta(:).';

    if isfield(PARAM.facet(iph),'A') && ~isempty(PARAM.facet(iph).A)
        facet_A = PARAM.facet(iph).A(:).';
    else
        facet_A = ones(size(facet_theta));
    end

    if isfield(PARAM.facet(iph),'sigma_ref') && ~isempty(PARAM.facet(iph).sigma_ref)
        sigma_ref = PARAM.facet(iph).sigma_ref;
    else
        sigma_ref = 1.0;
    end

    if isfield(PARAM.facet(iph),'A_ref') && ~isempty(PARAM.facet(iph).A_ref)
        A_ref = PARAM.facet(iph).A_ref;
    else
        A_ref = facet_A(1);
    end

    if isfield(PARAM.facet(iph),'q') && ~isempty(PARAM.facet(iph).q)
        q = PARAM.facet(iph).q;
    else
        q = q_default;
    end

else

    facet_theta = 2*pi*(0:nfold-1)/nfold;
    facet_A     = ones(size(facet_theta));
    sigma_ref   = 1.0;
    A_ref       = 1.0;
    q           = q_default;

end

if numel(facet_A) ~= numel(facet_theta)
    error('facet_A and facet_theta must have the same length.')
end

end


function [dtheta_best,ifacet_best] = Local_Nearest_Facet_Angle(theta,facet_theta)

%LOCAL_NEAREST_FACET_ANGLE Nearest unoriented facet normal.
%
% A crystallographic plane has two equivalent normals:
%
%   n and -n
%
% Therefore theta and theta+pi represent the same facet plane.
% The angular difference is wrapped to [-pi/2,pi/2].

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

%LOCAL_WRAPHALFPI Wrap angle difference for unoriented plane normals.
%
% Result is in [-pi/2,pi/2].
%
% This makes theta and theta+pi equivalent.

a = mod(a + pi/2,pi) - pi/2;

end

function a = Local_WrapPi(a)

a = mod(a + pi,2*pi) - pi;

end