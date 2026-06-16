function PARAM = Calc_AC_Anisotropy_FacetedStiffness(STATE,PARAM,MODEL,GRID)
%CALC_AC_ANISOTROPY_FACETEDSTIFFNESS Build anisotropic AC coefficients.
%
% Uses metadata stored in:
%
%   PARAM.facet(iph).mode
%   PARAM.facet(iph).theta
%   PARAM.facet(iph).theta_deg
%   PARAM.facet(iph).A_weight
%   PARAM.facet(iph).q
%   PARAM.facet(iph).nfold
%
% Phase index:
%
%   iph = MODEL.phase_index(ig)
%
% Output:
%
%   PARAM.L_AC(:,:,ig)
%   PARAM.Lm_AC(:,:,ig)
%   PARAM.LK_AC(:,:,ig)

[ny,nx,Ngrain] = size(STATE.phi);

% ------------------------------------------------------------
% Defaults
% ------------------------------------------------------------
q_default    = 0.25;
amin_default = 0.2;
amax_default = 10.0;
eps_grad     = 1e-20;

if ~isfield(PARAM,'L') || ~isfield(PARAM,'Lm') || ~isfield(PARAM,'LK')
    error('Calc_AC_Anisotropy_FacetedStiffness: PARAM.L, PARAM.Lm, PARAM.LK must exist.')
end

if ~isfield(PARAM,'theta_grain') || numel(PARAM.theta_grain) ~= Ngrain
    PARAM.theta_grain = zeros(1,Ngrain);
end

% ------------------------------------------------------------
% Start from isotropic values
% ------------------------------------------------------------
PARAM.L_AC  = zeros(ny,nx,Ngrain);
PARAM.Lm_AC = zeros(ny,nx,Ngrain);
PARAM.LK_AC = zeros(ny,nx,Ngrain);

PARAM.sigma_star_fac = ones(ny,nx,Ngrain);
PARAM.aniso_theta    = zeros(ny,nx,Ngrain);
PARAM.aniso_factor   = ones(ny,nx,Ngrain);

for ig = 1:Ngrain
    PARAM.L_AC(:,:,ig)  = PARAM.L;
    PARAM.Lm_AC(:,:,ig) = PARAM.Lm;
    PARAM.LK_AC(:,:,ig) = PARAM.LK;
end

if ~isfield(PARAM,'facet') || isempty(PARAM.facet)
    return
end

% ------------------------------------------------------------
% Apply phase-specific anisotropy
% ------------------------------------------------------------
for ig = 1:Ngrain

    iph = MODEL.phase_index(ig);

    if iph > numel(PARAM.facet)
        continue
    end

    facet = PARAM.facet(iph);

    if ~isfield(facet,'mode') || isempty(facet.mode)
        mode = 'iso';
    else
        mode = lower(facet.mode);
    end

    if strcmpi(mode,'iso')
        continue
    end

    q    = q_default;
    amin = amin_default;
    amax = amax_default;

    if isfield(facet,'q')    && ~isempty(facet.q),    q    = facet.q;    end
    if isfield(facet,'amin') && ~isempty(facet.amin), amin = facet.amin; end
    if isfield(facet,'amax') && ~isempty(facet.amax), amax = facet.amax; end

    phi = STATE.phi(:,:,ig);

    [phix,phiy] = gradient(phi,GRID.dx,GRID.dy);

    gmag = sqrt(phix.^2 + phiy.^2);
    mask = gmag > eps_grad;

    theta_lab = atan2(phiy,phix);
    theta_gr  = PARAM.theta_grain(ig);

    % Interface normal angle in crystal frame.
    % pi-periodic because plane normals are unoriented.
    theta_rel = mod(theta_lab - theta_gr,pi);

    A = ones(ny,nx);

    if strcmpi(mode,'simple')

        nfold = 4;

        if isfield(facet,'nfold') && ~isempty(facet.nfold)
            nfold = facet.nfold;
        end

        A(mask) = 1 + q*cos(nfold*theta_rel(mask));

    elseif strcmpi(mode,'facet')

        if isfield(facet,'theta') && ~isempty(facet.theta)

            facet_theta = facet.theta(:).';

        elseif isfield(facet,'theta_deg') && ~isempty(facet.theta_deg)

            facet_theta = facet.theta_deg(:).'*pi/180;

        else

            warning('Facet mode requested for phase %d but no facet angles were found.',iph)
            continue

        end

        if isfield(facet,'A_weight') && ~isempty(facet.A_weight)
            facet_A = facet.A_weight(:).';
        else
            facet_A = ones(size(facet_theta));
        end

        if numel(facet_A) ~= numel(facet_theta)
            error('Facet angle and A_weight length mismatch for phase %d.',iph)
        end

        % Smooth nearest-facet interpolation.
        % Smaller angular distance gives stronger contribution.
        ksharp = 25;
        Wsum   = zeros(ny,nx);
        Asum   = zeros(ny,nx);

        for jf = 1:numel(facet_theta)

            dtheta = Angle_Diff_Pi(theta_rel,facet_theta(jf));
            W      = exp(-ksharp*dtheta.^2);

            Wsum = Wsum + W;
            Asum = Asum + W*facet_A(jf);

        end

        A(mask) = Asum(mask)./max(Wsum(mask),eps);

    else

        continue

    end

    A = max(min(A,amax),amin);

    % Only change interface cells. Keep bulk exactly isotropic.
    A(~mask) = 1;

    PARAM.sigma_star_fac(:,:,ig) = A;
    PARAM.aniso_theta(:,:,ig)    = theta_rel;
    PARAM.aniso_factor(:,:,ig)   = A;

    % Apply anisotropy.
    % Lm controls barrier force, LK controls gradient force.
    % L itself is also scaled so chemical driving is consistent.
    PARAM.L_AC(:,:,ig)  = PARAM.L  .* A;
    PARAM.Lm_AC(:,:,ig) = PARAM.Lm .* A;
    PARAM.LK_AC(:,:,ig) = PARAM.LK .* A;

end

end


function d = Angle_Diff_Pi(a,b)
%ANGLE_DIFF_PI Smallest distance between unoriented angles with period pi.

d = abs(mod(a-b+pi/2,pi)-pi/2);

end