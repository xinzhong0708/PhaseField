function [phi] = PF_IMEX_Solver2D_AllenCahn_Periodic(LK,dx,dy,nx,ny,dt_phy,phi,s,mask)
% PF_IMEX_Solver2D_AllenCahn — periodic BC version, no nested functions
% (phi^{n+1}-phi^n)/dt = div( LK * grad phi^{n+1} ) + s
% Periodic in x and y. Only DOFs where mask(:,:,p)=true are solved.

mask = logical(mask);
Np   = size(phi,3);
LKv  = LK(:);

% periodic index wrappers
wrapI = @(i,sh) mod(i-1+sh, ny) + 1;
wrapJ = @(j,sh) mod(j-1+sh, nx) + 1;

tol_direct = 2000;

for p = 1:Np
    msk = mask(:,:,p);
    act = find(msk);
    Nact = numel(act);

    if Nact==0
        phi(:,:,p) = 0; continue;
    end

    % global->compact map
    g2c = zeros(ny*nx,1,'int32');
    g2c(act) = int32(1:Nact);

    % active (i,j)
    [ia,ja] = ind2sub([ny,nx], act);

    % periodic neighbors
    iL = ia;           jL = wrapJ(ja,-1);
    iR = ia;           jR = wrapJ(ja,+1);
    iU = wrapI(ia,-1); jU = ja;
    iD = wrapI(ia,+1); jD = ja;

    % linear neighbor indices
    idxC = act;
    idxL = sub2ind([ny,nx],iL,jL);
    idxR = sub2ind([ny,nx],iR,jR);
    idxU = sub2ind([ny,nx],iU,jU);
    idxD = sub2ind([ny,nx],iD,jD);

    % coefficients (face-averaged LK)
    LKc = LKv(idxC);
    cL  = -( LKv(idxL) + LKc )/2/dx^2;
    cR  = -( LKv(idxR) + LKc )/2/dx^2;
    cU  = -( LKv(idxU) + LKc )/2/dy^2;
    cD  = -( LKv(idxD) + LKc )/2/dy^2;
    cC  = -(cL + cR + cU + cD) + 1/dt_phy;

    % compact neighbor ids (0 if neighbor not active in mask)
    nL = g2c(idxL);
    nR = g2c(idxR);
    nU = g2c(idxU);
    nD = g2c(idxD);

    % estimate and allocate triplets
    est = 5*Nact;  % diag + up to 4 neighbors
    rows = zeros(est,1,'int32');
    cols = zeros(est,1,'int32');
    vals = zeros(est,1);
    ptr  = 1;

    comp = int32((1:Nact).');

    % diagonal
    rows(ptr:ptr+Nact-1) = comp;
    cols(ptr:ptr+Nact-1) = comp;
    vals(ptr:ptr+Nact-1) = cC;
    ptr = ptr + Nact;

    % add a neighbor block utility (inline)
    ok = nL>0; m = nnz(ok);
    if m>0
        rows(ptr:ptr+m-1) = comp(ok);
        cols(ptr:ptr+m-1) = double(nL(ok));
        vals(ptr:ptr+m-1) = cL(ok);
        ptr = ptr + m;
    end

    ok = nR>0; m = nnz(ok);
    if m>0
        rows(ptr:ptr+m-1) = comp(ok);
        cols(ptr:ptr+m-1) = double(nR(ok));
        vals(ptr:ptr+m-1) = cR(ok);
        ptr = ptr + m;
    end

    ok = nU>0; m = nnz(ok);
    if m>0
        rows(ptr:ptr+m-1) = comp(ok);
        cols(ptr:ptr+m-1) = double(nU(ok));
        vals(ptr:ptr+m-1) = cU(ok);
        ptr = ptr + m;
    end

    ok = nD>0; m = nnz(ok);
    if m>0
        rows(ptr:ptr+m-1) = comp(ok);
        cols(ptr:ptr+m-1) = double(nD(ok));
        vals(ptr:ptr+m-1) = cD(ok);
        ptr = ptr + m;
    end

    % trim triplets
    rows = rows(1:ptr-1); cols = cols(1:ptr-1); vals = vals(1:ptr-1);

    % RHS
    ph = phi(:,:,p); phv = ph(:);
    sv = s{p}(:);
    RHS = phv(act)/dt_phy + sv(act);

    % build & solve
    L = sparse(double(rows), double(cols), vals, Nact, Nact);
    if Nact <= tol_direct
        x = L \ RHS;
    else
        try
            setup_opts.type = 'ilutp'; setup_opts.droptol = 1e-3;
            [Lpre,Upre] = ilu(L, setup_opts);
            [x,flag] = bicgstab(L, RHS, 1e-10, 500, Lpre, Upre);
            if flag~=0, x = L\RHS; end
        catch
            x = L \ RHS;
        end
    end

    % write back (zero outside mask)
    out = zeros(ny*nx,1);
    out(act) = x;
    tmp = reshape(out,ny,nx);
    tmp(~msk) = 0;
    phi(:,:,p) = tmp;
end
end 