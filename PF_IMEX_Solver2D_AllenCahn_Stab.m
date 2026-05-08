function STATE = PF_IMEX_Solver2D_AllenCahn_Stab(STATE,MODEL,PARAM,GRID,NUM,Norm)
%PF_IMEX_SOLVER2D_ALLENCAHN_STAB
% Stabilized semi-implicit Allen-Cahn update.
%
% Original:
%   (1/dt - LK*Lap) phi_new = phi_old/dt + S_old
%
% Stabilized:
%   (1/dt + A_ac - LK*Lap) phi_new =
%       phi_old/dt + S_old + A_ac*phi_old
%
% A_ac is a positive stabilizer approximating local source stiffness.

LK     = PARAM.LK;
dx     = GRID.dx;
dy     = GRID.dy;
nx     = GRID.nx;
ny     = GRID.ny;
dt_phy = NUM.dt_phy;
phi    = STATE.phi;
s      = STATE.S_AC;
mask   = logical(STATE.mask);

Np         = size(phi,3);
LKvec_full = LK(:);

% -------------------------------------------------------------------------
% Stabilization field
% PARAM.A_ac can be:
%   scalar
%   ny x nx
%   ny x nx x Np
% PARAM.A_SC is accepted as alternative name
% -------------------------------------------------------------------------
if isfield(PARAM,'A_ac')
    Aac = PARAM.A_ac;
elseif isfield(PARAM,'A_SC')
    Aac = PARAM.A_SC;
else
    Aac = 0;
end

if isscalar(Aac)
    Aac = Aac*ones(ny,nx,Np);
elseif ndims(Aac) == 2
    Aac = repmat(Aac,1,1,Np);
end

tol_direct = 40000;
tol_iter   = 1e-10;
maxit      = 500;

for ip = 1:Np

    mask_p  = mask(:,:,ip);
    act_lin = find(mask_p);
    Nact    = numel(act_lin);

    if Nact == 0
        phi(:,:,ip) = 0;
        continue
    end

    phi_old = phi(:,:,ip);
    phi_vec = phi_old(:);
    s_vec   = s{ip}(:);

    % Stabilization vector for this phase
    Aac_p        = Aac(:,:,ip);
    Aac_vec_full = Aac_p(:);
    Aac_act      = Aac_vec_full(act_lin);

    [ia,ja] = ind2sub([ny,nx],act_lin);

    hasL = ja > 1;
    hasR = ja < nx;
    hasU = ia > 1;
    hasD = ia < ny;

    indL = zeros(Nact,1);
    indR = zeros(Nact,1);
    indU = zeros(Nact,1);
    indD = zeros(Nact,1);

    indL(hasL) = sub2ind([ny,nx],ia(hasL)  ,ja(hasL)-1);
    indR(hasR) = sub2ind([ny,nx],ia(hasR)  ,ja(hasR)+1);
    indU(hasU) = sub2ind([ny,nx],ia(hasU)-1,ja(hasU)  );
    indD(hasD) = sub2ind([ny,nx],ia(hasD)+1,ja(hasD)  );

    g2c = zeros(ny*nx,1,'int32');
    g2c(act_lin) = int32(1:Nact);

    comp = (1:Nact).';

    LKc = LKvec_full(act_lin);

    coefL = zeros(Nact,1);
    coefR = zeros(Nact,1);
    coefU = zeros(Nact,1);
    coefD = zeros(Nact,1);

    coefL(hasL) = -(LKvec_full(indL(hasL)) + LKc(hasL))/2/dx^2;
    coefR(hasR) = -(LKvec_full(indR(hasR)) + LKc(hasR))/2/dx^2;
    coefU(hasU) = -(LKvec_full(indU(hasU)) + LKc(hasU))/2/dy^2;
    coefD(hasD) = -(LKvec_full(indD(hasD)) + LKc(hasD))/2/dy^2;

    % ---------------------------------------------------------------------
    % Stabilized diagonal:
    %
    % old:
    %   coefC = -(coefL + coefR + coefU + coefD) + 1/dt
    %
    % new:
    %   coefC = -(coefL + coefR + coefU + coefD) + 1/dt + A_ac
    % ---------------------------------------------------------------------
    coefC = -(coefL + coefR + coefU + coefD) + 1/dt_phy + Aac_act;

    % ---------------------------------------------------------------------
    % Stabilized RHS:
    %
    % old:
    %   RHS = phi_old/dt + S_old
    %
    % new:
    %   RHS = phi_old/dt + S_old + A_ac*phi_old
    % ---------------------------------------------------------------------
    RHS = phi_vec(act_lin)/dt_phy + s_vec(act_lin) + Aac_act.*phi_vec(act_lin);

    est_nnz = 5*Nact;
    rows    = zeros(est_nnz,1);
    cols    = zeros(est_nnz,1);
    vals    = zeros(est_nnz,1);
    ptr     = 1;

    rows(ptr:ptr+Nact-1) = comp;
    cols(ptr:ptr+Nact-1) = comp;
    vals(ptr:ptr+Nact-1) = coefC;
    ptr = ptr + Nact;

    neigh = g2c(indL(hasL));
    rows0 = comp(hasL);
    vals0 = coefL(hasL);
    ok    = neigh > 0;
    nadd  = nnz(ok);
    if nadd > 0
        rows(ptr:ptr+nadd-1) = rows0(ok);
        cols(ptr:ptr+nadd-1) = double(neigh(ok));
        vals(ptr:ptr+nadd-1) = vals0(ok);
        ptr = ptr + nadd;
    end

    neigh = g2c(indR(hasR));
    rows0 = comp(hasR);
    vals0 = coefR(hasR);
    ok    = neigh > 0;
    nadd  = nnz(ok);
    if nadd > 0
        rows(ptr:ptr+nadd-1) = rows0(ok);
        cols(ptr:ptr+nadd-1) = double(neigh(ok));
        vals(ptr:ptr+nadd-1) = vals0(ok);
        ptr = ptr + nadd;
    end

    neigh = g2c(indU(hasU));
    rows0 = comp(hasU);
    vals0 = coefU(hasU);
    ok    = neigh > 0;
    nadd  = nnz(ok);
    if nadd > 0
        rows(ptr:ptr+nadd-1) = rows0(ok);
        cols(ptr:ptr+nadd-1) = double(neigh(ok));
        vals(ptr:ptr+nadd-1) = vals0(ok);
        ptr = ptr + nadd;
    end

    neigh = g2c(indD(hasD));
    rows0 = comp(hasD);
    vals0 = coefD(hasD);
    ok    = neigh > 0;
    nadd  = nnz(ok);
    if nadd > 0
        rows(ptr:ptr+nadd-1) = rows0(ok);
        cols(ptr:ptr+nadd-1) = double(neigh(ok));
        vals(ptr:ptr+nadd-1) = vals0(ok);
        ptr = ptr + nadd;
    end

    rows = rows(1:ptr-1);
    cols = cols(1:ptr-1);
    vals = vals(1:ptr-1);

    Lmat = sparse(rows,cols,vals,Nact,Nact);

    x0 = phi_vec(act_lin);

    if Nact <= tol_direct

        x = Lmat \ RHS;

    else

        try
            opts.type     = 'ict';
            opts.droptol  = 1e-3;
            opts.diagcomp = 1e-3;

            Lpre = ichol(Lmat,opts);

            [x,flag] = pcg(Lmat,RHS,tol_iter,maxit,Lpre,Lpre',x0);

            if flag ~= 0
                warning('pcg not converged, use direct solver');
                x = Lmat \ RHS;
            end

        catch

            x = Lmat \ RHS;

        end

    end

    phi_new           = zeros(ny*nx,1);
    phi_new(act_lin)  = x;

    tmp               = reshape(phi_new,ny,nx);
    tmp(~mask_p)      = 0;

    phi(:,:,ip)       = tmp;

end

% Normalization
if Norm == 1
    phi = Norm_Phi(phi);
end

STATE.phi = phi;

% Update p
STATE.p = Calc_p(MODEL,phi);

end