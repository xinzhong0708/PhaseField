function PARAM = Calc_AC_Anisotropy_Simple(STATE,PARAM,GRID)
%CALC_AC_ANISOTROPY_SIMPLE Simple orientation-dependent AC interface energy.
%
% Adds:
%   PARAM.LK_AC(:,:,alpha)
%   PARAM.Lm_AC(:,:,alpha)
%
% Only affects Allen-Cahn interface terms.
% CH diffusion and CH fourth-order terms are unchanged.

[ny,nx,Ngrain] = size(STATE.phi);

delta = 0.03;
nfold = 4;
grad_min = 1e-12;

if isfield(PARAM,'aniso_delta')
    delta = PARAM.aniso_delta;
end
if isfield(PARAM,'aniso_nfold')
    nfold = PARAM.aniso_nfold;
end
if isfield(PARAM,'aniso_grad_min')
    grad_min = PARAM.aniso_grad_min;
end

theta0 = zeros(1,Ngrain);

if isfield(PARAM,'theta_grain')
    theta0 = PARAM.theta_grain;
end

PARAM.LK_AC = zeros(ny,nx,Ngrain);
PARAM.Lm_AC = zeros(ny,nx,Ngrain);

for alpha = 1:Ngrain

    phi = STATE.phi(:,:,alpha);

    % Reflective central differences
    if nx == 1
        phi_L = phi;
        phi_R = phi;
    else
        phi_L = phi(:,[2,1:nx-1]);
        phi_R = phi(:,[2:nx,nx-1]);
    end

    if ny == 1
        phi_U = phi;
        phi_D = phi;
    else
        phi_U = phi([2,1:ny-1],:);
        phi_D = phi([2:ny,ny-1],:);
    end

    gx = (phi_R - phi_L)/(2*GRID.dx);
    gy = (phi_D - phi_U)/(2*GRID.dy);

    gnorm = sqrt(gx.^2 + gy.^2);

    theta = atan2(gy,gx) - theta0(alpha);

    a = 1 + delta*cos(nfold*theta);

    % Safety bound
    a = max(a,0.2);

    % In phase interior, orientation is meaningless
    a(gnorm < grad_min) = 1;

    PARAM.LK_AC(:,:,alpha) = PARAM.LK .* a;
    PARAM.Lm_AC(:,:,alpha) = PARAM.Lm .* a;

end

end