clear; clf; colormap(jet(256))
load temp

alpha        = 3;
dt           = NUM.dt_phy;
dx           = GRID.dx;
dy           = GRID.dy;

NUM.norm_phi = 0;
NUM.cut_phi  = 0;
NUM.norm_E   = 0;

NUM.use_Aac  = 0;
NUM.use_Jphi = 1;
NUM.Jphi_eps = 1e-6;
NUM.phi_mask_thick =  500;


% ------------------------------------------------------------
% Old-state AC source
% ------------------------------------------------------------
PARAM.A_ac   = zeros(GRID.ny,GRID.nx);

STATE_SRC    = Calc_S_AllenCahn(STATE_OLD,PARAM,MODEL);
S_old        = STATE_SRC.S_AC;

% Two outputs are required to get Jphi_diag and maskPhi
[STATE,DIAG] = PF_Coupled_ACCH_LETangent_complete( ...
                    STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);

dphi = STATE.phi(:,:,alpha) - STATE_OLD.phi(:,:,alpha);

% ------------------------------------------------------------
% Chemical dmu tangent
% This direct dpdphi form is appropriate for checking.
% ------------------------------------------------------------
chem = zeros(size(dphi));

for ip = 1:size(STATE_OLD.p,3)

    dpdphi = MODEL.dpdphi(alpha,ip,STATE_OLD.phi);

    for ie = 1:numel(STATE_OLD.mu_e)

        dmu  = STATE.mu_e{ie} - STATE_OLD.mu_e{ie};
        chem = chem - PARAM.L .* dpdphi .* STATE_OLD.e{ip}{ie} .* dmu;

    end

end

% ------------------------------------------------------------
% Diagonal real phi Jacobian used in the new AC matrix
% ------------------------------------------------------------
Jaa = DIAG.Jphi_diag{alpha};

% ------------------------------------------------------------
% New AC equation checked:
%
% dphi/dt - LK*Lap(dphi) - Jaa*dphi + chem
%       = LK*Lap(phi_old) + S_old
% ------------------------------------------------------------
lhs = dphi/dt ...
    - PARAM.LK .* Laplacian_Reflex(dphi,dx,dy) ...
    - Jaa .* dphi ...
    + chem;

rhs = PARAM.LK .* Laplacian_Reflex(STATE_OLD.phi(:,:,alpha),dx,dy) ...
    + S_old{alpha};

res = lhs - rhs;

% ------------------------------------------------------------
% Only active AC rows were solved
% ------------------------------------------------------------
mask = DIAG.maskPhi(:,:,alpha);

res_plot = res;
res_plot(~mask) = NaN;

figure(1); clf
pcolor(res_plot)
colorbar
shading interp
title('AC LHS - RHS on solved mask')

fprintf('AC max residual, active rows = %.3e\n',max(abs(res(mask))));
fprintf('AC relative residual         = %.3e\n', ...
    max(abs(res(mask))) / max(max(abs(rhs(mask))),eps));





% CH equation checker
ieq    = 1;
Ne     = numel(STATE_OLD.E);
Np     = size(STATE_OLD.p,3);
Nop    = size(STATE_OLD.phi,3);
dE_tan = zeros(size(STATE_OLD.E{ieq}));
for je = 1:Ne
    dmu_j  = STATE.mu_e{je} - STATE_OLD.mu_e{je};
    dE_tan = dE_tan + STATE_OLD.chi{ieq,je} .* dmu_j;
end
for alpha = 1:Nop
    dphi_a = STATE.phi(:,:,alpha) - STATE_OLD.phi(:,:,alpha);
    B_a = zeros(size(dphi_a));
    for ip = 1:Np
        dpdphi = MODEL.dpdphi(alpha,ip,STATE_OLD.phi);
        B_a    = B_a + STATE_OLD.e{ip}{ieq} .* dpdphi;
    end
    dE_tan = dE_tan + B_a .* dphi_a;
end

% ------------------------------------------------------------
M        =  PARAM.M{ieq};
dmu      =  STATE.mu_e{ieq} - STATE_OLD.mu_e{ieq};
mu_old   =  STATE_OLD.mu_e{ieq};
E_old    =  STATE_OLD.E{ieq};

D_dmu    = -M .* Laplacian_Reflex(dmu,dx,dy);
D_mu_old = -M .* Laplacian_Reflex(mu_old,dx,dy);

Q_dE     =  M .* PARAM.kappa_eff .* Laplacian_Reflex(Laplacian_Reflex(dE_tan,dx,dy),dx,dy);
Q_E_old  =  M .* PARAM.kappa_eff .* Laplacian_Reflex(Laplacian_Reflex(E_old,dx,dy),dx,dy);

% ------------------------------------------------------------
lhs_ch   =  dE_tan/dt + Q_dE + D_dmu;
rhs_ch   = -Q_E_old - D_mu_old;
res_ch   =  lhs_ch - rhs_ch;
  
figure(2); clf
pcolor(res_ch); colorbar; shading interp; title('CH LHS - RHS')






function Lap = Laplacian_Reflex(A,dx,dy)
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
Lap = (AL - 2*A + AR)/dx^2 + (AU - 2*A + AD)/dy^2;
end