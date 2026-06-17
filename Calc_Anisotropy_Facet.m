function [PARAM] = Calc_Anisotropy_Facet(STATE,PARAM,GRID,PHYS)

% ------------------------------------------------------------
%Step 1: Calculate the gradient of phi
% ------------------------------------------------------------
phi        = STATE.phi;
[ny,nx,Np] = size(phi);

% Grid spacing
dx = GRID.dx;
dy = GRID.dy;

% Allocate
PARAM.aniso_phix  = zeros(ny,nx,Np);
PARAM.aniso_phiy  = zeros(ny,nx,Np);
PARAM.aniso_gmag  = zeros(ny,nx,Np);
PARAM.aniso_nx    = zeros(ny,nx,Np);
PARAM.aniso_ny    = zeros(ny,nx,Np);
PARAM.aniso_theta = zeros(ny,nx,Np);

% Loop through phases
for ip = 1:Np

    % Current phi
    ph = phi(:,:,ip);

    % Compute gradient
    [phix,phiy] = Grad2D_Reflect(ph,dx,dy);

    gmag  = sqrt(phix.^2 + phiy.^2);
    nxloc = zeros(ny,nx);
    nyloc = zeros(ny,nx);

    % Mask it
    mask = gmag > 1e-20;

    nxloc(mask) = phix(mask)./gmag(mask);
    nyloc(mask) = phiy(mask)./gmag(mask);

    % Angle of interface normal
    theta = atan2d(nyloc,nxloc);
    theta = mod(theta,360);

    PARAM.aniso_phix(:,:,ip)  = phix;
    PARAM.aniso_phiy(:,:,ip)  = phiy;
    PARAM.aniso_gmag(:,:,ip)  = gmag;
    PARAM.aniso_nx(:,:,ip)    = nxloc;
    PARAM.aniso_ny(:,:,ip)    = nyloc;
    PARAM.aniso_theta(:,:,ip) = theta;

end


% ------------------------------------------------------------
%Step 2: Find closest facet normal
% ------------------------------------------------------------
PARAM.aniso_facet_id     = zeros(ny,nx,Np);
PARAM.aniso_facet_theta  = zeros(ny,nx,Np);
PARAM.aniso_facet_A      = ones(ny,nx,Np);
PARAM.aniso_theta_diff   = zeros(ny,nx,Np);

% Grain rotation
theta_grain = zeros(1,Np);

if isfield(PARAM,'theta_grain')
    theta_grain(1:min(Np,numel(PARAM.theta_grain))) = PARAM.theta_grain(1:min(Np,numel(PARAM.theta_grain)));
end

% Loop through phases
for ip = 1:Np

    % Skip isotropic phase
    if ~isfield(PARAM,'facet') || ip > numel(PARAM.facet)
        continue
    end

    if isempty(PARAM.facet(ip).theta_deg)
        continue
    end

    if strcmpi(PARAM.facet(ip).mode,'iso')
        continue
    end

    % Facet normal angle from metadata
    facet_theta = PARAM.facet(ip).theta_deg(:);

    % A weight from metadata
    facet_A = PARAM.facet(ip).A(:);

    % Add opposite normals
    facet_theta = [facet_theta; facet_theta + 180];
    facet_A     = [facet_A;     facet_A];

    % Rotate by grain angle
    facet_theta = mod(facet_theta + theta_grain(ip),360);

    Nfacet = numel(facet_theta);
    theta_n = PARAM.aniso_theta(:,:,ip);

    best_diff  = inf(ny,nx);
    best_id    = zeros(ny,nx);
    best_theta = zeros(ny,nx);
    best_A     = ones(ny,nx);

    for jf = 1:Nfacet

        dth  = AngleDiffDeg(theta_n,facet_theta(jf));
        take = dth < best_diff;

        best_diff(take)  = dth(take);
        best_id(take)    = jf;
        best_theta(take) = facet_theta(jf);
        best_A(take)     = facet_A(jf);

    end

    % Only meaningful at interface
    mask = PARAM.aniso_gmag(:,:,ip) > 1e-6;

    PARAM.aniso_facet_id(:,:,ip)     = best_id    .* mask;
    PARAM.aniso_facet_theta(:,:,ip)  = best_theta .* mask;
    PARAM.aniso_facet_A(:,:,ip)      = best_A;
    PARAM.aniso_theta_diff(:,:,ip)   = best_diff  .* mask;

end


% ------------------------------------------------------------
%Step 3: Calculate Kundin-type stiffness factor
% ------------------------------------------------------------
PARAM.aniso_S       = ones(ny,nx,Np);
PARAM.aniso_S_raw   = ones(ny,nx,Np);
PARAM.aniso_mask0   = false(ny,nx,Np);
PARAM.aniso_mask    = false(ny,nx,Np);


% Gradient cutoff
g_rel = 0.02;
phi_cut = 1e-3;

if isfield(PARAM,'aniso_g_rel')
    g_rel = PARAM.aniso_g_rel;
end

if isfield(PARAM,'aniso_phi_cut')
    phi_cut = PARAM.aniso_phi_cut;
end

% Expand interface mask by this many grids
ngrow = 0;

% Loop through phases
for ip = 1:Np

    % Skip isotropic phase
    if ~isfield(PARAM,'facet') || ip > numel(PARAM.facet)
        continue
    end

    if isempty(PARAM.facet(ip).theta_deg)
        continue
    end

    if strcmpi(PARAM.facet(ip).mode,'iso')
        continue
    end

    % Kundin q parameter
    q = PARAM.facet(ip).q;

    % Safety limits
    amin = PARAM.facet(ip).amin;
    amax = PARAM.facet(ip).amax;

    % Facet area weight
    A_use = PARAM.aniso_facet_A(:,:,ip);
    A_ref = PARAM.facet(ip).A_ref;

    if A_ref <= 0
        A_ref = 1.0;
    end

    % Angle to closest facet normal
    dth = PARAM.aniso_theta_diff(:,:,ip);

    % Kundin surface-energy factor, not stiffness.
    % sigma/sigma_ref = (A_ref/A_use) * sqrt(sin^2 + q^2 cos^2) / q
    den   = sind(dth).^2 + q^2*cosd(dth).^2;
    S_raw = (A_ref ./ max(A_use,eps)) .* sqrt(den) ./ max(q,eps);


    % Cap for numerical safety
    S_raw = max(S_raw,amin);
    S_raw = min(S_raw,amax);

    % Interface mask from relative gradient and diffuse phi band
    gmax = max(PARAM.aniso_gmag(:,:,ip),[],'all');

    if gmax > 0
        mask_g = PARAM.aniso_gmag(:,:,ip) > g_rel*gmax;
    else
        mask_g = false(ny,nx);
    end

    mask_phi = phi(:,:,ip) > phi_cut & phi(:,:,ip) < 1 - phi_cut;

    mask0 = mask_g | mask_phi;

    % Expand mask without periodic boundary
    % [S_ext,mask] = Expand_Field_Mask(S_raw,mask0,ngrow);
    S_ext = S_raw;
    mask  = mask0;

    % Apply anisotropy only in expanded interface band
    S = ones(ny,nx);
    S(mask) = S_ext(mask);

    % Store
    PARAM.aniso_S_raw(:,:,ip) = S_raw;
    PARAM.aniso_S(:,:,ip)     = S;
    PARAM.aniso_mask0(:,:,ip) = mask0;
    PARAM.aniso_mask(:,:,ip)  = mask;

end


% ------------------------------------------------------------
%Step 4: Calculate anisotropic LK and Lm
% ------------------------------------------------------------
PARAM.LK_aniso  = zeros(ny,nx,Np);
PARAM.Lm_aniso  = zeros(ny,nx,Np);
PARAM.LK_iso    = zeros(ny,nx,Np);
PARAM.Lm_iso    = zeros(ny,nx,Np);
PARAM.L_AC      = zeros(ny,nx,Np);
PARAM.m_AC_base = zeros(ny,nx,Np);

for ip = 1:Np

    % Base L
    if isscalar(PARAM.L)
        Lloc = PARAM.L*ones(ny,nx);
    elseif ndims(PARAM.L) == 2
        Lloc = PARAM.L;
    else
        Lloc = PARAM.L(:,:,ip);
    end

    % Base kap
    if isscalar(PHYS.kap)
        kaploc = PHYS.kap*ones(ny,nx);
    elseif isvector(PHYS.kap)
        kaploc = PHYS.kap(ip)*ones(ny,nx);
    else
        kaploc = PHYS.kap(:,:,ip);
    end

    % Base m
    if isscalar(PHYS.m)
        mloc = PHYS.m*ones(ny,nx);
    elseif isvector(PHYS.m)
        mloc = PHYS.m(ip)*ones(ny,nx);
    else
        mloc = PHYS.m(:,:,ip);
    end

    % Apply anisotropic stiffness
    S = PARAM.aniso_S(:,:,ip);

    PARAM.LK_iso(:,:,ip)    = Lloc .* kaploc;
    PARAM.Lm_iso(:,:,ip)    = Lloc .* mloc;

    PARAM.LK_aniso(:,:,ip)  = PARAM.LK_iso(:,:,ip) .* S;
    PARAM.Lm_aniso(:,:,ip)  = PARAM.Lm_iso(:,:,ip) .* S;

    % Chemical driving force should use unscaled L
    PARAM.L_AC(:,:,ip)      = Lloc;

    % Base m used to remove m from MODEL.dgdphi if needed
    PARAM.m_AC_base(:,:,ip) = mloc;

end

% Assign 3-D fields for AC solver
PARAM.LK_AC = PARAM.LK_aniso;
PARAM.Lm_AC = PARAM.Lm_aniso;

end


%% HELPERS
% ------------------------------------------------------------

function [fx,fy] = Grad2D_Reflect(f,dx,dy)

[ny,nx] = size(f);

fx = zeros(ny,nx);
fy = zeros(ny,nx);

% Interior central difference
fx(:,2:nx-1)   = (f(:,3:nx) - f(:,1:nx-2))/(2*dx);
fy(2:ny-1,:)   = (f(3:ny,:) - f(1:ny-2,:))/(2*dy);

% Boundary one-sided difference
fx(:,1)        = (f(:,2) - f(:,1))/dx;
fx(:,nx)       = (f(:,nx) - f(:,nx-1))/dx;

fy(1,:)        = (f(2,:) - f(1,:))/dy;
fy(ny,:)       = (f(ny,:) - f(ny-1,:))/dy;

end


function dth = AngleDiffDeg(theta1,theta2)

dth = abs(atan2d(sind(theta1 - theta2),cosd(theta1 - theta2)));

end


function [field,mask] = Expand_Field_Mask(field,mask,ngrow)

ker = ones(3,3);

for ig = 1:ngrow

    cnt  = conv2(double(mask),ker,'same');
    sumv = conv2(field.*double(mask),ker,'same');

    newmask = cnt > 0;
    fill    = newmask & ~mask;

    field(fill) = sumv(fill)./cnt(fill);
    mask        = newmask;

end

end