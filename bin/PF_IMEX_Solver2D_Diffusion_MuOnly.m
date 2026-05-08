function STATE = PF_IMEX_Solver2D_Diffusion_MuOnly(STATE,PARAM,GRID,PHYS,NUM)
%PF_IMEX_SOLVER2D_DIFFUSION_MUONLY Structured mu-only diffusion solver.

%Unpack
M       = PARAM.M;
kappa   = PHYS.kappa;
dx      = GRID.dx;
dy      = GRID.dy;
nx      = GRID.nx;
ny      = GRID.ny;
E       = STATE.E;
mu      = STATE.mu_e;
Chi     = STATE.chi;
s       = STATE.S_DIFF;
dt_phy  = NUM.dt_phy;
% x0      = STATE.sol_diff;

%Settings
chiRelFloor = 0;
chiAbsFloor = 0;

use_iter = false;
tol      = 1e-12;
maxit    = 500;

%Basic sizes
Nc     = length(Chi);
Nnodes = nx * ny;
Ntot   = Nc * Nnodes;

E_old  = E;
mu_old = mu;

dx2 = dx^2;
dy2 = dy^2;
dx4 = dx2 * dx2;
dy4 = dy2 * dy2;

% -------------------------------------------------------------------------
% Grid / neighbor indexing
% -------------------------------------------------------------------------
[Igrid,Jgrid] = ndgrid(1:ny,1:nx);
ii = Igrid(:);
jj = Jgrid(:);

refI = @(i,sh) reflect_index(i + sh, ny);
refJ = @(j,sh) reflect_index(j + sh, nx);

jjL  = refJ(jj,-1);   jjR  = refJ(jj,+1);
jjL2 = refJ(jj,-2);   jjR2 = refJ(jj,+2);

iiU  = refI(ii,-1);   iiD  = refI(ii,+1);
iiU2 = refI(ii,-2);   iiD2 = refI(ii,+2);

iiUR = refI(ii,-1);   jjUR = refJ(jj,+1);
iiDR = refI(ii,+1);   jjDR = refJ(jj,+1);
iiUL = refI(ii,-1);   jjUL = refJ(jj,-1);
iiDL = refI(ii,+1);   jjDL = refJ(jj,-1);

idx_c  = sub2ind([ny,nx], ii  , jj  );
idx_L  = sub2ind([ny,nx], ii  , jjL );
idx_R  = sub2ind([ny,nx], ii  , jjR );
idx_U  = sub2ind([ny,nx], iiU , jj  );
idx_D  = sub2ind([ny,nx], iiD , jj  );

idx_L2 = sub2ind([ny,nx], ii  , jjL2);
idx_R2 = sub2ind([ny,nx], ii  , jjR2);
idx_U2 = sub2ind([ny,nx], iiU2, jj  );
idx_D2 = sub2ind([ny,nx], iiD2, jj  );

idx_UR = sub2ind([ny,nx], iiUR, jjUR);
idx_DR = sub2ind([ny,nx], iiDR, jjDR);
idx_UL = sub2ind([ny,nx], iiUL, jjUL);
idx_DL = sub2ind([ny,nx], iiDL, jjDL);

%Global unknown ids
idNode = cell(1,Nc);
for l = 1:Nc
    idNode{l} = (l:Nc:Ntot).';
end

% -------------------------------------------------------------------------
% Chi regularization
% -------------------------------------------------------------------------
if chiRelFloor == 0 && chiAbsFloor == 0
    ChiReg = Chi;
else
    ChiReg = Regularize_Chi(Chi,Nc,ny,nx,idx_c,chiRelFloor,chiAbsFloor);
end

% -------------------------------------------------------------------------
% Sparse assembly
% -------------------------------------------------------------------------
max_nnz = Nnodes * (5*Nc + 13*Nc*Nc);

rows = zeros(max_nnz,1);
cols = zeros(max_nnz,1);
vals = zeros(max_nnz,1);
R    = zeros(Ntot,1);

k = 1;
n = Nnodes;

for l = 1:Nc

    idMu_c = idNode{l}(idx_c);
    idMu_L = idNode{l}(idx_L);
    idMu_R = idNode{l}(idx_R);
    idMu_U = idNode{l}(idx_U);
    idMu_D = idNode{l}(idx_D);

    Ml   = M{l};
    Ml_c = Ml(idx_c);
    Ml_L = Ml(idx_L);
    Ml_R = Ml(idx_R);
    Ml_U = Ml(idx_U);
    Ml_D = Ml(idx_D);

    %B_l = Eo_l - sum_m Chi_lm * muo_m
    B_l = E_old{l};

    for m = 1:Nc
        B_l = B_l - ChiReg{l,m} .* mu_old{m};
    end

    Bc  = B_l(idx_c);
    BL  = B_l(idx_L);
    BR  = B_l(idx_R);
    BU  = B_l(idx_U);
    BD  = B_l(idx_D);

    BL2 = B_l(idx_L2);
    BR2 = B_l(idx_R2);
    BU2 = B_l(idx_U2);
    BD2 = B_l(idx_D2);

    BUR = B_l(idx_UR);
    BDR = B_l(idx_DR);
    BUL = B_l(idx_UL);
    BDL = B_l(idx_DL);

    %Diffusion operator on mu_l: -div(M grad mu)
    d_L = -(Ml_L + Ml_c) / 2 / dx2;
    d_R = -(Ml_R + Ml_c) / 2 / dx2;
    d_U = -(Ml_U + Ml_c) / 2 / dy2;
    d_D = -(Ml_D + Ml_c) / 2 / dy2;
    d_C = -(d_L + d_R + d_U + d_D);

    %A_E operator on E_l, fourth-order CH term
    q_L  = Ml_c * kappa .* ( -4/dx4 - 4/(dx2*dy2) );
    q_R  = Ml_c * kappa .* ( -4/dx4 - 4/(dx2*dy2) );
    q_U  = Ml_c * kappa .* ( -4/dy4 - 4/(dx2*dy2) );
    q_D  = Ml_c * kappa .* ( -4/dy4 - 4/(dx2*dy2) );

    q_L2 = Ml_c * kappa .* (  1/dx4 );
    q_R2 = Ml_c * kappa .* (  1/dx4 );
    q_U2 = Ml_c * kappa .* (  1/dy4 );
    q_D2 = Ml_c * kappa .* (  1/dy4 );

    q_UR = Ml_c * kappa .* (  2/(dx2*dy2) );
    q_DR = Ml_c * kappa .* (  2/(dx2*dy2) );
    q_UL = Ml_c * kappa .* (  2/(dx2*dy2) );
    q_DL = Ml_c * kappa .* (  2/(dx2*dy2) );

    q_C  = Ml_c * kappa .* (  6/dx4 + 6/dy4 + 8/(dx2*dy2) );

    aC   = 1/dt_phy + q_C;
    aL   = q_L;     aR   = q_R;     aU   = q_U;     aD   = q_D;
    aL2  = q_L2;    aR2  = q_R2;    aU2  = q_U2;    aD2  = q_D2;
    aUR  = q_UR;    aDR  = q_DR;    aUL  = q_UL;    aDL  = q_DL;

    %RHS = Eo/dt - s - A_E * B_l
    Eo    = E_old{l}(idx_c);
    s_loc = s{l}(idx_c);

    R_mu = Eo/dt_phy - s_loc ...
         - aC  .* Bc ...
         - aL  .* BL  - aR  .* BR  - aU  .* BU  - aD  .* BD ...
         - aL2 .* BL2 - aR2 .* BR2 - aU2 .* BU2 - aD2 .* BD2 ...
         - aUR .* BUR - aDR .* BDR - aUL .* BUL - aDL .* BDL;

    R(idMu_c) = R_mu;

    %Diffusion block on mu_l
    rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idMu_L; vals(k:k+n-1) = d_L; k = k+n;
    rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idMu_R; vals(k:k+n-1) = d_R; k = k+n;
    rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idMu_U; vals(k:k+n-1) = d_U; k = k+n;
    rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idMu_D; vals(k:k+n-1) = d_D; k = k+n;
    rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idMu_c; vals(k:k+n-1) = d_C; k = k+n;

    %A_E * Chi contribution
    for m = 1:Nc

        idm_c  = idNode{m}(idx_c);
        idm_L  = idNode{m}(idx_L);
        idm_R  = idNode{m}(idx_R);
        idm_U  = idNode{m}(idx_U);
        idm_D  = idNode{m}(idx_D);

        idm_L2 = idNode{m}(idx_L2);
        idm_R2 = idNode{m}(idx_R2);
        idm_U2 = idNode{m}(idx_U2);
        idm_D2 = idNode{m}(idx_D2);

        idm_UR = idNode{m}(idx_UR);
        idm_DR = idNode{m}(idx_DR);
        idm_UL = idNode{m}(idx_UL);
        idm_DL = idNode{m}(idx_DL);

        Chi_lm = ChiReg{l,m};

        Chi_c  = Chi_lm(idx_c);
        Chi_L  = Chi_lm(idx_L);
        Chi_R  = Chi_lm(idx_R);
        Chi_U  = Chi_lm(idx_U);
        Chi_D  = Chi_lm(idx_D);

        Chi_L2 = Chi_lm(idx_L2);
        Chi_R2 = Chi_lm(idx_R2);
        Chi_U2 = Chi_lm(idx_U2);
        Chi_D2 = Chi_lm(idx_D2);

        Chi_UR = Chi_lm(idx_UR);
        Chi_DR = Chi_lm(idx_DR);
        Chi_UL = Chi_lm(idx_UL);
        Chi_DL = Chi_lm(idx_DL);

        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_c;  vals(k:k+n-1) = aC  .* Chi_c;  k = k+n;
        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_L;  vals(k:k+n-1) = aL  .* Chi_L;  k = k+n;
        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_R;  vals(k:k+n-1) = aR  .* Chi_R;  k = k+n;
        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_U;  vals(k:k+n-1) = aU  .* Chi_U;  k = k+n;
        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_D;  vals(k:k+n-1) = aD  .* Chi_D;  k = k+n;

        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_L2; vals(k:k+n-1) = aL2 .* Chi_L2; k = k+n;
        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_R2; vals(k:k+n-1) = aR2 .* Chi_R2; k = k+n;
        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_U2; vals(k:k+n-1) = aU2 .* Chi_U2; k = k+n;
        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_D2; vals(k:k+n-1) = aD2 .* Chi_D2; k = k+n;

        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_UR; vals(k:k+n-1) = aUR .* Chi_UR; k = k+n;
        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_DR; vals(k:k+n-1) = aDR .* Chi_DR; k = k+n;
        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_UL; vals(k:k+n-1) = aUL .* Chi_UL; k = k+n;
        rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idm_DL; vals(k:k+n-1) = aDL .* Chi_DL; k = k+n;

    end

end

rows = rows(1:k-1);
cols = cols(1:k-1);
vals = vals(1:k-1);

%Assemble and solve
Lmu = sparse(rows, cols, vals, Ntot, Ntot);

if use_iter

    if isempty(x0) || numel(x0) ~= Ntot
        x0 = zeros(size(R));

        for l = 1:Nc
            x0(l:Nc:end) = mu_old{l}(:);
        end
    else
        x0 = x0(:);
    end

    [sol,flag,relres,iter] = bicgstab(Lmu,R,tol,maxit,[],[],x0);

    if flag ~= 0 || any(~isfinite(sol))
        warning('bicgstab failed: flag=%d, relres=%g, iter=%d. Falling back to direct.', ...
            flag, relres, iter);
        sol = Lmu \ R;
    end

else

    sol = Lmu \ R;

end

%Unpack mu
for l = 1:Nc
    mu{l} = reshape(sol(l:Nc:end), ny, []);
end

%Reconstruct E from closure
for l = 1:Nc

    En = E_old{l};

    for m = 1:Nc
        En = En - ChiReg{l,m} .* mu_old{m} + ChiReg{l,m} .* mu{m};
    end

    E{l} = En;

end

%Enforce plain mean conservation with minimal change
E = EnforceMeanE(E, E_old);

%Write back
STATE.E        = E;
STATE.mu_e     = mu;
STATE.sol_diff = sol;

end


% =========================================================================
% Helper functions
% =========================================================================
function idx = reflect_index(idx, n)

if n == 1
    idx = ones(size(idx));
    return
end

period = 2*n - 2;
r      = mod(idx - 1, period);
idx    = 1 + min(r, period - r);

end

function ChiReg = Regularize_Chi(Chi,Nc,ny,nx,idx_all,chiRelFloor,chiAbsFloor)
Nnodes = numel(idx_all);
ChiReg = cell(Nc,Nc);
for a = 1:Nc
    for b = 1:Nc
        ChiReg{a,b} = zeros(ny,nx);
    end
end
for q = 1:Nnodes
    ChiMat = zeros(Nc,Nc);
    for a = 1:Nc
        for b = 1:Nc
            ChiMat(a,b) = Chi{a,b}(idx_all(q));
        end
    end
    ChiMat = 0.5 * (ChiMat + ChiMat.');
    [V,D]  = eig(ChiMat);
    lam    = diag(D);
    smax   = max(abs(lam));
    floorv = max(chiAbsFloor, chiRelFloor * max(smax, eps));
    lam2   = sign(lam) .* max(abs(lam), floorv);
    ChiLoc = V * diag(lam2) * V.';
    [iv,jv] = ind2sub([ny,nx], idx_all(q));
    for a = 1:Nc
        for b = 1:Nc
            ChiReg{a,b}(iv,jv) = ChiLoc(a,b);
        end
    end
end
end


function E = EnforceMeanE(E, E_old)
Nc = numel(E);
for a = 1:Nc
    target_mean = mean(E_old{a}(:));
    new_mean    = mean(E{a}(:));
    E{a} = E{a} + (target_mean - new_mean);
end
end