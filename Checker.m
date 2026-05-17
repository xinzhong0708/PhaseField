clear;clf;colormap(jet(256))
load 100

alpha        = 4;
dt           = NUM.dt_phy;
dx           = GRID.dx;
dy           = GRID.dy;
NUM.norm_phi = 0;
NUM.norm_E   = 0;

% Recompute source and A_ac from old state
PARAM.A_ac   = Calc_Aac_FrozenOmega(STATE_OLD,PARAM,MODEL,3,1e-6,0,[]);
STATE_SRC    = Calc_S_AllenCahn(STATE_OLD,PARAM,MODEL);
S_old        = STATE_SRC.S_AC;
STATE        = PF_Coupled_ACCH_LETangent(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);

dphi         = STATE.phi(:,:,alpha) - STATE_OLD.phi(:,:,alpha);

chem = zeros(size(dphi));
for ip = 1:size(STATE_OLD.p,3)
    dpdphi = MODEL.dpdphi(alpha,ip,STATE_OLD.phi);
    for ie = 1:numel(STATE_OLD.mu_e)
        dmu  = STATE.mu_e{ie} - STATE_OLD.mu_e{ie};
        chem = chem - PARAM.L.*dpdphi.*STATE_OLD.e{ip}{ie}.*dmu;
    end
end

lhs = dphi/dt - PARAM.LK .* Laplacian(dphi,dx,dy) + PARAM.A_ac.*dphi + chem;

rhs = PARAM.LK.*Laplacian(STATE_OLD.phi(:,:,alpha),dx,dy) + S_old{alpha};

res = lhs - rhs;

pcolor(res);colorbar;shading interp




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

D_dmu    = -M .* Laplacian(dmu,dx,dy);
D_mu_old = -M .* Laplacian(mu_old,dx,dy);

Q_dE     =  M .* PHYS.kappa .* Laplacian(Laplacian(dE_tan,dx,dy),dx,dy);
Q_E_old  =  M .* PHYS.kappa .* Laplacian(Laplacian(E_old,dx,dy),dx,dy);

% ------------------------------------------------------------
lhs_ch   =  dE_tan/dt + Q_dE + D_dmu;
rhs_ch   = -Q_E_old - D_mu_old;
res_ch   =  lhs_ch - rhs_ch;
  
figure(2); clf
pcolor(res_ch); colorbar; shading interp; title('CH LHS - RHS')
