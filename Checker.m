clear; clf; colormap(jet(256))

%Load data
load temp

%% Solve it with the solver
NUM.norm_phi       =  0;
NUM.cut_phi        =  0;
NUM.norm_E         =  0;
NUM.use_Jphi       =  1;
STATE_OLD          =  STATE;
Ne                 =  numel(STATE_OLD.E);
PARAM.kappa_eff(:) =  mean(PARAM.kappa_eff(:),'all');

NUM.use_mu_band     = 1;
NUM.mu_band_thick   = 8;
NUM.mu_p_cut        = 1e-8;
NUM.use_order_cache = 1;
tic
[STATE,DIAG]        = PF_Coupled_ACCH_LETangent_CS(STATE_OLD,PARAM,MODEL,GRID,PHYS,NUM);
toc



%% 1) AC CHECK
alpha        = 3;
dt           = NUM.dt_phy;
dx           = GRID.dx;
dy           = GRID.dy;
dphi         = STATE.phi(:,:,alpha) - STATE_OLD.phi(:,:,alpha);

%Chemical part of the source term
chem         = zeros(size(dphi));
for ip = 1:size(STATE_OLD.p,3)
    dpdphi = MODEL.dpdphi(alpha,ip,STATE_OLD.phi);
    for ie = 1:numel(STATE_OLD.mu_e)
        dmu  = STATE.mu_e{ie} - STATE_OLD.mu_e{ie};
        chem = chem - PARAM.L .* dpdphi .* STATE_OLD.e{ip}{ie} .* dmu;
    end
end
Jaa          = DIAG.Jphi_diag{alpha};
STATE_SRC    = Calc_S_AllenCahn(STATE_OLD,PARAM,MODEL);
S_old        = STATE_SRC.S_AC;
lhs          = dphi/dt - PARAM.LK .* Laplacian_Reflex(dphi,dx,dy) - Jaa .* dphi + chem;
rhs          = PARAM.LK .* Laplacian_Reflex(STATE_OLD.phi(:,:,alpha),dx,dy) + S_old{alpha};

subplot(121);pcolor(lhs-rhs);shading interp;colorbar;title('AC residual')



%% 2) CH CHECK
ieq    = 2;
dmu    = cell(1,Ne);
for ie = 1:Ne
    dmu{ie} = STATE.mu_e{ie} - STATE_OLD.mu_e{ie};
end

% Rebuild dE exactly as used in CH matrix
dE     = zeros(size(STATE_OLD.E{ieq}));

% chi part
for je = 1:Ne
    dmu_j = STATE.mu_e{je} - STATE_OLD.mu_e{je};
    dE    = dE + STATE_OLD.chi{ieq,je} .* dmu_j;
end

% if using convex split chi, this should use chi_imp, not raw chi
if isfield(PARAM,'use_CS_chi') && PARAM.use_CS_chi == 1
    chi_imp = ConvexifyChi_ForImplicit(STATE_OLD.chi,PARAM);
    dE      = zeros(size(STATE_OLD.E{ieq}));
    for je = 1:Ne
        dmu_j = STATE.mu_e{je} - STATE_OLD.mu_e{je};
        dE    = dE + chi_imp{ieq,je} .* dmu_j;
    end
end

% phase-fraction tangent part: same formula as solver
eps_phi = 1e-14;
Dphi    = sum(STATE_OLD.phi.^2,3) + eps_phi;
p_tan   = STATE_OLD.phi.^2 ./ Dphi;
e_bar   = zeros(size(dE));
for ig = 1:size(STATE_OLD.p,3)
    e_bar = e_bar + p_tan(:,:,ig) .* STATE_OLD.e{ig}{ieq};
end

for alpha = 1:size(STATE_OLD.phi,3)
    dphi_a = STATE.phi(:,:,alpha) - STATE_OLD.phi(:,:,alpha);
    fac    = 2 .* STATE_OLD.phi(:,:,alpha) ./ Dphi;
    B_a    = fac .* (STATE_OLD.e{alpha}{ieq} - e_bar);
    dE     = dE + B_a .* dphi_a;
end


M      = PARAM.M{ieq};
Dmu    = Diff_Reflex(STATE.mu_e{ieq},M,dx,dy);
Kc     = KappaC_FromKAPC(ieq,dmu,M,PARAM.kappa_eff,DIAG.KAPC,dx,dy);
lhs_ch = dE/dt + Dmu + Kc;
rhs_ch = zeros(size(lhs_ch));
res_ch = lhs_ch - rhs_ch;

subplot(122)
pcolor(res_ch); shading interp; colorbar
title('CH residual')




%% ========================================================================
%% HELPERS
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



function Df = Diff_Reflex(f,M,dx,dy)

[ny,nx] = size(f);

if nx == 1
    fL = f; fR = f; ML = M; MR = M;
else
    fL = f(:,[2,1:nx-1]);
    fR = f(:,[2:nx,nx-1]);
    ML = M(:,[2,1:nx-1]);
    MR = M(:,[2:nx,nx-1]);
end

if ny == 1
    fU = f; fD = f; MU = M; MD = M;
else
    fU = f([2,1:ny-1],:);
    fD = f([2:ny,ny-1],:);
    MU = M([2,1:ny-1],:);
    MD = M([2:ny,ny-1],:);
end

dL = -(ML+M)/2/dx^2;
dR = -(MR+M)/2/dx^2;
dU = -(MU+M)/2/dy^2;
dD = -(MD+M)/2/dy^2;
dC = -(dL+dR+dU+dD);

Df = dC.*f + dL.*fL + dR.*fR + dU.*fU + dD.*fD;

end



function chi_imp = ConvexifyChi_ForImplicit(chi_raw,PARAM)
%CONVEXIFYCHI_FORIMPLICIT Positive bounded capacity for implicit CH terms.
%
% chi_raw is the physical LE susceptibility dE/dmu. In spinodal regions it
% can be indefinite or too large. chi_imp is used only in the implicit CH
% storage term and E recovery of the tangent/corrector. The raw mu_e field
% is still used in the CH residual, so the spinodal driving force is kept.

Ne = size(chi_raw,1);
[ny,nx] = size(chi_raw{1,1});
N = ny*nx;

chi_floor = 1e-8;
chi_cap   = inf;

if isfield(PARAM,'CS_chi_floor')
    chi_floor = PARAM.CS_chi_floor;
end
if isfield(PARAM,'CS_chi_cap')
    chi_cap = PARAM.CS_chi_cap;
end

chi_imp = cell(Ne,Ne);

for i = 1:Ne
    for j = 1:Ne
        chi_imp{i,j} = zeros(ny,nx);
    end
end

for id = 1:N

    C = zeros(Ne,Ne);

    for i = 1:Ne
        for j = 1:Ne
            tmp = chi_raw{i,j};
            C(i,j) = tmp(id);
        end
    end

    C = 0.5*(C+C.');

    [V,D] = eig(C);
    lam = diag(D);

    scale = max(1,max(abs(lam)));

    % Positive capacity for the implicit response.
    % Negative curvature is kept explicitly through raw mu_e in the residual.
    lam_imp = max(abs(lam),chi_floor*scale);

    if isfinite(chi_cap)
        lam_imp = min(lam_imp,chi_cap);
    end

    Cimp = V*diag(lam_imp)*V.';
    Cimp = 0.5*(Cimp+Cimp.');

    for i = 1:Ne
        for j = 1:Ne
            tmp = chi_imp{i,j};
            tmp(id) = Cimp(i,j);
            chi_imp{i,j} = tmp;
        end
    end
end

end


function Kc = KappaC_FromKAPC(ieq,dmu,M,kappa_eff,KAPC,dx,dy)

[ny,nx] = size(kappa_eff);
Ne      = numel(dmu);
Kc      = zeros(ny,nx);

for ik = 1:numel(KAPC)

    p_ig = KAPC(ik).p;
    K    = M .* kappa_eff .* p_ig;

    if max(abs(K(:))) == 0
        continue
    end

    for ic = 1:KAPC(ik).Nc
        cnew = reshape(KAPC(ik).c(ic,:),ny,nx);
        for je = 1:Ne
            Xij  = reshape(KAPC(ik).X(ic,je,:),ny,nx);
            cnew = cnew + Xij .* dmu{je};
        end
        Jij = reshape(KAPC(ik).J(ieq,ic,:),ny,nx);
        Kc = Kc + Jij .* Kappa4_Apply_13Point(K,cnew,dx,dy);
    end
end

end



function Q = Kappa4_Apply_13Point(K,A,dx,dy)
%KAPPA4_APPLY_13POINT Apply variable-coefficient lap(K*lap(A)).
%
% This matches the 13-point stencil used in PF_Coupled_ACCH_LETangent_CS.
%
% Q = lap(K * lap(A))

[ny,nx] = size(A);

[Igrid,Jgrid] = ndgrid(1:ny,1:nx);
ii = Igrid(:);
jj = Jgrid(:);

refI = @(i,sh) Reflect_Index(i+sh,ny);
refJ = @(j,sh) Reflect_Index(j+sh,nx);

jjL  = refJ(jj,-1);
jjR  = refJ(jj,+1);
iiU  = refI(ii,-1);
iiD  = refI(ii,+1);

jjL2 = refJ(jj,-2);
jjR2 = refJ(jj,+2);
iiU2 = refI(ii,-2);
iiD2 = refI(ii,+2);

iiUR = refI(ii,-1); jjUR = refJ(jj,+1);
iiDR = refI(ii,+1); jjDR = refJ(jj,+1);
iiUL = refI(ii,-1); jjUL = refJ(jj,-1);
iiDL = refI(ii,+1); jjDL = refJ(jj,-1);

idx_c  = sub2ind([ny,nx],ii,jj);
idx_L  = sub2ind([ny,nx],ii,jjL);
idx_R  = sub2ind([ny,nx],ii,jjR);
idx_U  = sub2ind([ny,nx],iiU,jj);
idx_D  = sub2ind([ny,nx],iiD,jj);

idx_L2 = sub2ind([ny,nx],ii,jjL2);
idx_R2 = sub2ind([ny,nx],ii,jjR2);
idx_U2 = sub2ind([ny,nx],iiU2,jj);
idx_D2 = sub2ind([ny,nx],iiD2,jj);

idx_UR = sub2ind([ny,nx],iiUR,jjUR);
idx_DR = sub2ind([ny,nx],iiDR,jjDR);
idx_UL = sub2ind([ny,nx],iiUL,jjUL);
idx_DL = sub2ind([ny,nx],iiDL,jjDL);

dx2 = dx^2;
dy2 = dy^2;
dx4 = dx2^2;
dy4 = dy2^2;

K_c = K(idx_c);
K_L = K(idx_L);
K_R = K(idx_R);
K_U = K(idx_U);
K_D = K(idx_D);

q_L  = -2*(K_L+K_c)/dx4 - 2*(K_L+K_c)/(dx2*dy2);
q_R  = -2*(K_R+K_c)/dx4 - 2*(K_R+K_c)/(dx2*dy2);
q_U  = -2*(K_U+K_c)/dy4 - 2*(K_U+K_c)/(dx2*dy2);
q_D  = -2*(K_D+K_c)/dy4 - 2*(K_D+K_c)/(dx2*dy2);

q_L2 = K_L/dx4;
q_R2 = K_R/dx4;
q_U2 = K_U/dy4;
q_D2 = K_D/dy4;

q_UR = (K_U+K_R)/(dx2*dy2);
q_DR = (K_D+K_R)/(dx2*dy2);
q_UL = (K_U+K_L)/(dx2*dy2);
q_DL = (K_D+K_L)/(dx2*dy2);

q_C  = (K_L+K_R)/dx4 ...
     + (K_U+K_D)/dy4 ...
     + 4*K_c/dx4 ...
     + 4*K_c/dy4 ...
     + 8*K_c/(dx2*dy2);

Qv = ...
    q_C .*A(idx_c)  + ...
    q_L .*A(idx_L)  + q_R .*A(idx_R)  + ...
    q_U .*A(idx_U)  + q_D .*A(idx_D)  + ...
    q_L2.*A(idx_L2) + q_R2.*A(idx_R2) + ...
    q_U2.*A(idx_U2) + q_D2.*A(idx_D2) + ...
    q_UR.*A(idx_UR) + q_DR.*A(idx_DR) + ...
    q_UL.*A(idx_UL) + q_DL.*A(idx_DL);

Q = reshape(Qv,ny,nx);

end



function idx = Reflect_Index(idx,n)

if n == 1
    idx = ones(size(idx));
    return
end

period = 2*n - 2;
r      = mod(idx - 1,period);
idx    = 1 + min(r,period - r);

end
