function [E, mu] = PF_IMEX_Solver2D_Diffusion_Periodic(M, kappa, dx, dy, nx, ny, E, mu, dt_phy, Chi, s, chiRelFloor, chiAbsFloor)
% Solver for mixed formulation (E, mu_E) for KKS model with CH term (4th order)
%
%   dE/dt = div(M*grad(mu_E)) - M*kappa*grad4(E) - s
%   Closure:
%       E^n - E^o = Chi*(mu^n - mu^o)
%
% Inputs:
%   M      : cell{Nc} of mobility fields, size [ny,nx]
%   kappa  : CH gradient coefficient (scalar)
%   dx,dy  : grid spacings
%   nx,ny  : number of points in x,y
%   E      : cell{Nc} of E at old time step (overwritten)
%   mu     : cell{Nc} of mu at old time step (overwritten)
%   dt_phy : physical time step
%   Chi    : cell{Nc,Nc} susceptibility matrix fields Chi_{l,m}(y,x)
%   s      : optional source term
%            - [] (default, no source)
%            - cell{Nc}, each [ny,nx]
%            - numeric [ny,nx,Nc]
%   chiRelFloor, chiAbsFloor : regularization floors for Chi
%
% Outputs:
%   E, mu  : updated fields at new time step

if nargin < 11 || isempty(s), s = []; end
if nargin < 12 || isempty(chiRelFloor), chiRelFloor = 0e-12; end
if nargin < 13 || isempty(chiAbsFloor), chiAbsFloor = 0e-14; end

Nc   = length(Chi);

Njp  = 2 * Nc;
Ntot = Njp * nx * ny;

[Igrid, Jgrid] = ndgrid(1:ny, 1:nx);
ii             = Igrid(:);
jj             = Jgrid(:);
Nnodes         = numel(ii);

gid   = @(loc, iv, jv) (((jv-1)*ny + (iv-1)) * Njp) + loc;

wrapI = @(i,sh) mod(i-1+sh, ny) + 1;
wrapJ = @(j,sh) mod(j-1+sh, nx) + 1;

jjL   = wrapJ(jj,-1);  jjR   = wrapJ(jj,+1);
jjL2  = wrapJ(jj,-2);  jjR2  = wrapJ(jj,+2);
iiU   = wrapI(ii,-1);  iiD   = wrapI(ii,+1);
iiU2  = wrapI(ii,-2);  iiD2  = wrapI(ii,+2);

iiUR  = wrapI(ii,-1);  jjUR  = wrapJ(jj,+1);
iiDR  = wrapI(ii,+1);  jjDR  = wrapJ(jj,+1);
iiUL  = wrapI(ii,-1);  jjUL  = wrapJ(jj,-1);
iiDL  = wrapI(ii,+1);  jjDL  = wrapJ(jj,-1);

max_nnz = 50 * Ntot;
rows    = zeros(max_nnz,1);
cols    = zeros(max_nnz,1);
vals    = zeros(max_nnz,1);
k       = 1;
n       = Nnodes;

R       = zeros(Ntot,1);
dx2     = dx^2;
dy2     = dy^2;
dx4     = dx2*dx2;
dy4     = dy2*dy2;

E_old   = E;
mu_old  = mu;

% ---- normalize/interpret s ----
% Accept:
%   - [] (no source)
%   - cell{1,Nc} or cell{Nc,1}, each [ny,nx]
%   - numeric [ny,nx,Nc]
if isempty(s)
    s = cell(1,Nc);
    for a = 1:Nc
        s{a} = zeros(ny,nx);
    end
elseif ~iscell(s)
    if isnumeric(s) && ndims(s) == 3 && size(s,1) == ny && size(s,2) == nx && size(s,3) == Nc
        sCell = cell(1,Nc);
        for a = 1:Nc
            sCell{a} = s(:,:,a);
        end
        s = sCell;
    else
        error('s must be [] or cell{Nc} with [ny,nx] fields or numeric [ny,nx,Nc].');
    end
else
    if numel(s) ~= Nc
        error('s must have Nc components.');
    end
    s = reshape(s,1,[]);
    for a = 1:Nc
        if ~isequal(size(s{a}), [ny,nx])
            error('s{%d} must have size [ny,nx].', a);
        end
    end
end

% ---- Chi regularization at every grid point (same spirit as 1D solver) ----
idx_all = sub2ind([ny,nx], ii, jj);
ChiReg  = cell(Nc,Nc);
for a = 1:Nc
    for b = 1:Nc
        ChiReg{a,b} = zeros(Nnodes,1);
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

    for a = 1:Nc
        for b = 1:Nc
            ChiReg{a,b}(q) = ChiLoc(a,b);
        end
    end
end

for l = 1:Nc
    locE  = l;
    locMu = Nc + l;

    idE_c   = gid(locE , ii  , jj  );
    idE_L   = gid(locE , ii  , jjL );
    idE_R   = gid(locE , ii  , jjR );
    idE_U   = gid(locE , iiU , jj  );
    idE_D   = gid(locE , iiD , jj  );
    idE_L2  = gid(locE , ii  , jjL2);
    idE_R2  = gid(locE , ii  , jjR2);
    idE_U2  = gid(locE , iiU2, jj  );
    idE_D2  = gid(locE , iiD2, jj  );
    idE_UR  = gid(locE , iiUR, jjUR);
    idE_DR  = gid(locE , iiDR, jjDR);
    idE_UL  = gid(locE , iiUL, jjUL);
    idE_DL  = gid(locE , iiDL, jjDL);

    idMu_c  = gid(locMu, ii  , jj  );
    idMu_L  = gid(locMu, ii  , jjL );
    idMu_R  = gid(locMu, ii  , jjR );
    idMu_U  = gid(locMu, iiU , jj  );
    idMu_D  = gid(locMu, iiD , jj  );

    idx_c  = sub2ind([ny,nx], ii  , jj  );

    Ml    = M{l};
    Ml_c  = Ml(idx_c);
    Ml_L  = Ml(sub2ind([ny,nx], ii  , jjL ));
    Ml_R  = Ml(sub2ind([ny,nx], ii  , jjR ));
    Ml_U  = Ml(sub2ind([ny,nx], iiU , jj  ));
    Ml_D  = Ml(sub2ind([ny,nx], iiD , jj  ));

    % ============================================================
    %  1) E-equation
    % ============================================================

    coef_time = 1/dt_phy;
    rows(k:k+n-1) = idE_c;
    cols(k:k+n-1) = idE_c;
    vals(k:k+n-1) = coef_time;
    k = k + n;

    d_L = -(Ml_L + Ml_c)/2/dx2;
    d_R = -(Ml_R + Ml_c)/2/dx2;
    d_U = -(Ml_U + Ml_c)/2/dy2;
    d_D = -(Ml_D + Ml_c)/2/dy2;
    d_C = -(d_L + d_R + d_U + d_D);

    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idMu_L; vals(k:k+n-1) = d_L; k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idMu_R; vals(k:k+n-1) = d_R; k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idMu_U; vals(k:k+n-1) = d_U; k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idMu_D; vals(k:k+n-1) = d_D; k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idMu_c; vals(k:k+n-1) = d_C; k = k+n;

    q_L  =  Ml_c*kappa .* ( -4/dx4 - 4/(dx2*dy2) );
    q_L2 =  Ml_c*kappa .* (  1/dx4 );
    q_R  =  Ml_c*kappa .* ( -4/dx4 - 4/(dx2*dy2) );
    q_R2 =  Ml_c*kappa .* (  1/dx4 );
    q_U  =  Ml_c*kappa .* ( -4/dy4 - 4/(dx2*dy2) );
    q_U2 =  Ml_c*kappa .* (  1/dy4 );
    q_D  =  Ml_c*kappa .* ( -4/dy4 - 4/(dx2*dy2) );
    q_D2 =  Ml_c*kappa .* (  1/dy4 );
    q_UR =  Ml_c*kappa .* (  2/(dx2*dy2) );
    q_DR =  Ml_c*kappa .* (  2/(dx2*dy2) );
    q_UL =  Ml_c*kappa .* (  2/(dx2*dy2) );
    q_DL =  Ml_c*kappa .* (  2/(dx2*dy2) );
    q_C  =  Ml_c*kappa .* (  6/dx4 + 6/dy4 + 8/(dx2*dy2) );

    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_L;   vals(k:k+n-1) = q_L;   k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_R;   vals(k:k+n-1) = q_R;   k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_U;   vals(k:k+n-1) = q_U;   k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_D;   vals(k:k+n-1) = q_D;   k = k+n;

    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_L2;  vals(k:k+n-1) = q_L2;  k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_R2;  vals(k:k+n-1) = q_R2;  k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_U2;  vals(k:k+n-1) = q_U2;  k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_D2;  vals(k:k+n-1) = q_D2;  k = k+n;

    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_UR;  vals(k:k+n-1) = q_UR;  k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_DR;  vals(k:k+n-1) = q_DR;  k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_UL;  vals(k:k+n-1) = q_UL;  k = k+n;
    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_DL;  vals(k:k+n-1) = q_DL;  k = k+n;

    rows(k:k+n-1) = idE_c; cols(k:k+n-1) = idE_c;   vals(k:k+n-1) = q_C;   k = k+n;

    Eo    = E_old{l}(idx_c);
    s_loc = s{l}(idx_c);
    R_El  = Eo/dt_phy - s_loc;
    R(idE_c) = R(idE_c) + R_El;

    % ============================================================
    %  2) Closure equation (Chi-form)
    % ============================================================

    Elo  = E_old{l}(idx_c);
    R_cl = Elo;

    % +1 * E_l^n
    idE_l = gid(locE, ii, jj);
    rows(k:k+n-1) = idMu_c;
    cols(k:k+n-1) = idE_l;
    vals(k:k+n-1) = 1;
    k = k + n;

    % -ChiReg_{l,m} * mu_m^n  ; RHS: -ChiReg_{l,m} * mu_m^o
    for m = 1:Nc
        Chi_lm = ChiReg{l,m};
        mu_mo  = mu_old{m}(idx_c);

        locMu_m = Nc + m;
        idMu_m  = gid(locMu_m, ii, jj);

        rows(k:k+n-1) = idMu_c;
        cols(k:k+n-1) = idMu_m;
        vals(k:k+n-1) = -Chi_lm;
        k = k + n;

        R_cl = R_cl - Chi_lm .* mu_mo;
    end

    R(idMu_c) = R(idMu_c) + R_cl;
end

rows = rows(1:k-1);
cols = cols(1:k-1);
vals = vals(1:k-1);

%Assemble L
L   = sparse(rows, cols, vals, Ntot, Ntot);

%Direct solver
% tic
% sol = L\R;
% toc

tic
p = colamd(L);
sol = zeros(size(R));
sol(p) = real((L(:,p)) \ R);
toc


for l = 1:Nc
    locE  = l;
    locMu = Nc + l;
    E{l}  = reshape(sol(locE:Njp:end),  ny, []);
    mu{l} = reshape(sol(locMu:Njp:end), ny, []);
end
end