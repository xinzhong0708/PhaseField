function [E, mu] = PF_IMEX_Solver2D_Diffusion_MuOnly(M, kappa, dx, dy, nx, ny, E, mu, dt_phy, Chi, s, chiRelFloor, chiAbsFloor)
% Solver for mixed formulation (E, mu_E) for KKS model with CH term (4th order)
% Solves only for mu, then reconstructs E from closure:
%
%   dE/dt = div(M*grad(mu_E)) - M*kappa*grad4(E) - s
%   Closure:
%       E^n - E^o = Chi*(mu^n - mu^o)
%
% so that
%       E^n = E^o - Chi*mu^o + Chi*mu^n
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

Nc    = length(Chi);
Ntot  = Nc * nx * ny;

[Igrid, Jgrid] = ndgrid(1:ny, 1:nx);
ii             = Igrid(:);
jj             = Jgrid(:);
Nnodes         = numel(ii);

gid   = @(loc, iv, jv) (((jv-1)*ny + (iv-1)) * Nc) + loc;

% Reflective indexing for zero-flux boundaries
refI = @(i,sh) reflect_index(i + sh, ny);
refJ = @(j,sh) reflect_index(j + sh, nx);

jjL   = refJ(jj,-1);  jjR   = refJ(jj,+1);
jjL2  = refJ(jj,-2);  jjR2  = refJ(jj,+2);
iiU   = refI(ii,-1);  iiD   = refI(ii,+1);
iiU2  = refI(ii,-2);  iiD2  = refI(ii,+2);

iiUR  = refI(ii,-1);  jjUR  = refJ(jj,+1);
iiDR  = refI(ii,+1);  jjDR  = refJ(jj,+1);
iiUL  = refI(ii,-1);  jjUL  = refJ(jj,-1);
iiDL  = refI(ii,+1);  jjDL  = refJ(jj,-1);

max_nnz = (20*Nc + 8) * Ntot;
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

% ---- Chi regularization at every grid point ----
idx_all = sub2ind([ny,nx], ii, jj);
ChiReg  = cell(Nc,Nc);
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

for l = 1:Nc
    idMu_c  = gid(l, ii  , jj  );
    idMu_L  = gid(l, ii  , jjL );
    idMu_R  = gid(l, ii  , jjR );
    idMu_U  = gid(l, iiU , jj  );
    idMu_D  = gid(l, iiD , jj  );

    idx_c   = sub2ind([ny,nx], ii  , jj  );
    idx_L   = sub2ind([ny,nx], ii  , jjL );
    idx_R   = sub2ind([ny,nx], ii  , jjR );
    idx_U   = sub2ind([ny,nx], iiU , jj  );
    idx_D   = sub2ind([ny,nx], iiD , jj  );
    idx_L2  = sub2ind([ny,nx], ii  , jjL2);
    idx_R2  = sub2ind([ny,nx], ii  , jjR2);
    idx_U2  = sub2ind([ny,nx], iiU2, jj  );
    idx_D2  = sub2ind([ny,nx], iiD2, jj  );
    idx_UR  = sub2ind([ny,nx], iiUR, jjUR);
    idx_DR  = sub2ind([ny,nx], iiDR, jjDR);
    idx_UL  = sub2ind([ny,nx], iiUL, jjUL);
    idx_DL  = sub2ind([ny,nx], iiDL, jjDL);

    Ml    = M{l};
    Ml_c  = Ml(idx_c);
    Ml_L  = Ml(idx_L);
    Ml_R  = Ml(idx_R);
    Ml_U  = Ml(idx_U);
    Ml_D  = Ml(idx_D);

    % ============================================================
    %  Precompute closure RHS field:
    %      B_l = Eo_l - sum_m ChiReg_{l,m} * muo_m
    % ============================================================
    B_l = E_old{l};
    for m = 1:Nc
        B_l = B_l - ChiReg{l,m} .* mu_old{m};
    end

    Bc   = B_l(idx_c);
    BL   = B_l(idx_L);
    BR   = B_l(idx_R);
    BU   = B_l(idx_U);
    BD   = B_l(idx_D);
    BL2  = B_l(idx_L2);
    BR2  = B_l(idx_R2);
    BU2  = B_l(idx_U2);
    BD2  = B_l(idx_D2);
    BUR  = B_l(idx_UR);
    BDR  = B_l(idx_DR);
    BUL  = B_l(idx_UL);
    BDL  = B_l(idx_DL);

    % ============================================================
    %  Mu-only equation obtained by eliminating E
    % ============================================================

    % D operator on mu_l
    d_L = -(Ml_L + Ml_c)/2/dx2;
    d_R = -(Ml_R + Ml_c)/2/dx2;
    d_U = -(Ml_U + Ml_c)/2/dy2;
    d_D = -(Ml_D + Ml_c)/2/dy2;
    d_C = -(d_L + d_R + d_U + d_D);

    % A_E operator on E_l
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

    aC   = 1/dt_phy + q_C;
    aL   = q_L;
    aR   = q_R;
    aU   = q_U;
    aD   = q_D;
    aL2  = q_L2;
    aR2  = q_R2;
    aU2  = q_U2;
    aD2  = q_D2;
    aUR  = q_UR;
    aDR  = q_DR;
    aUL  = q_UL;
    aDL  = q_DL;

    % RHS = Eo/dt - s - A_E * B_l
    Eo    = E_old{l}(idx_c);
    s_loc = s{l}(idx_c);

    R_mu  = Eo/dt_phy - s_loc ...
          - aC  .* Bc  ...
          - aL  .* BL  - aR  .* BR  - aU  .* BU  - aD  .* BD ...
          - aL2 .* BL2 - aR2 .* BR2 - aU2 .* BU2 - aD2 .* BD2 ...
          - aUR .* BUR - aDR .* BDR - aUL .* BUL - aDL .* BDL;

    R(idMu_c) = R(idMu_c) + R_mu;

    % Assemble mu matrix:
    % diffusion part on mu_l
    rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idMu_L; vals(k:k+n-1) = d_L; k = k+n;
    rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idMu_R; vals(k:k+n-1) = d_R; k = k+n;
    rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idMu_U; vals(k:k+n-1) = d_U; k = k+n;
    rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idMu_D; vals(k:k+n-1) = d_D; k = k+n;
    rows(k:k+n-1) = idMu_c; cols(k:k+n-1) = idMu_c; vals(k:k+n-1) = d_C; k = k+n;

    % A_E * Chi contribution on all mu_m
    for m = 1:Nc
        idm_c   = gid(m, ii  , jj  );
        idm_L   = gid(m, ii  , jjL );
        idm_R   = gid(m, ii  , jjR );
        idm_U   = gid(m, iiU , jj  );
        idm_D   = gid(m, iiD , jj  );
        idm_L2  = gid(m, ii  , jjL2);
        idm_R2  = gid(m, ii  , jjR2);
        idm_U2  = gid(m, iiU2, jj  );
        idm_D2  = gid(m, iiD2, jj  );
        idm_UR  = gid(m, iiUR, jjUR);
        idm_DR  = gid(m, iiDR, jjDR);
        idm_UL  = gid(m, iiUL, jjUL);
        idm_DL  = gid(m, iiDL, jjDL);

        Chi_c   = ChiReg{l,m}(idx_c);
        Chi_L   = ChiReg{l,m}(idx_L);
        Chi_R   = ChiReg{l,m}(idx_R);
        Chi_U   = ChiReg{l,m}(idx_U);
        Chi_D   = ChiReg{l,m}(idx_D);
        Chi_L2  = ChiReg{l,m}(idx_L2);
        Chi_R2  = ChiReg{l,m}(idx_R2);
        Chi_U2  = ChiReg{l,m}(idx_U2);
        Chi_D2  = ChiReg{l,m}(idx_D2);
        Chi_UR  = ChiReg{l,m}(idx_UR);
        Chi_DR  = ChiReg{l,m}(idx_DR);
        Chi_UL  = ChiReg{l,m}(idx_UL);
        Chi_DL  = ChiReg{l,m}(idx_DL);

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

% Assemble mu-only system
Lmu = sparse(rows, cols, vals, Ntot, Ntot);

% Solve for mu only
sol = Lmu \ R;

for l = 1:Nc
    mu{l} = reshape(sol(l:Nc:end), ny, []);
end

% Reconstruct E from closure:
%    E^n = E^o - ChiReg*mu^o + ChiReg*mu^n
for l = 1:Nc
    En = E_old{l};
    for m = 1:Nc
        En = En - ChiReg{l,m} .* mu_old{m} + ChiReg{l,m} .* mu{m};
    end
    E{l} = En;
end

% Reconstruct E from closure:
%    E^n = E^o - ChiReg*mu^o + ChiReg*mu^n
for l = 1:Nc
    En = E_old{l};
    for m = 1:Nc
        En = En - ChiReg{l,m} .* mu_old{m} + ChiReg{l,m} .* mu{m};
    end
    E{l} = En;
end

% Enforce plain mean conservation with minimal change
E = EnforceMeanE(E, E_old, dt_phy);

end


function idx = reflect_index(idx, n)
% Reflective indexing:
% ... 3 2 | 1 2 3 ... n-2 n-1 n | n-1 n-2 ...
if n == 1
    idx = ones(size(idx));
    return
end

period = 2*n - 2;
r      = mod(idx - 1, period);
idx    = 1 + min(r, period - r);
end



function E = EnforceMeanE(E, E_old, varargin)
%ENFORCEMEANE Enforce exact conservation of total bulk composition.
%
% For closed / periodic / no-flux systems:
%
%     sum(E_new{a}) = sum(E_old{a})
%
% for every conserved component a.
%
% This is appropriate when the source term s comes from phase-fraction
% redistribution and should not create or destroy total composition.

Nc = numel(E);

for a = 1:Nc
    target_mean = mean(E_old{a}(:));
    new_mean    = mean(E{a}(:));

    E{a} = E{a} + (target_mean - new_mean);
end

end