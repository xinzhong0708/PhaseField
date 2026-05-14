function STATE_NEW = PF_IMEX_Solver2D_AllenCahn_StabRef(STATE_OLD,STATE_REF,MODEL,PARAM,GRID,NUM,Norm)
%PF_IMEX_SOLVER2D_ALLENCAHN_STABREF
%
% Reference-linearized stabilized Allen-Cahn step.
%
% It solves one physical timestep from STATE_OLD, but evaluates the
% nonlinear AC source and stabilizer at STATE_REF:
%
%   (1/dt - div(LK grad) + A_ref) phi_new
%       = phi_old/dt + S_ref + A_ref*phi_ref
%
% This is different from the usual stabilized AC step, where the same
% state is used both as old-time value and source reference.

%Source and stabilizer from reference state
STATE_SRC       =  STATE_REF;
STATE_SRC       =  Calc_S_AllenCahn(STATE_SRC,PARAM,MODEL);

A_ac            =  Calc_Aac_FrozenOmega(STATE_SRC,PARAM,MODEL,3,1e-6,0,[]);

%Unpack
LK              =  PARAM.LK;
dx              =  GRID.dx;
dy              =  GRID.dy;
nx              =  GRID.nx;
ny              =  GRID.ny;
dt_phy          =  NUM.dt_phy;

phi_old_all     =  STATE_OLD.phi;
phi_ref_all     =  STATE_REF.phi;
s               =  STATE_SRC.S_AC;

if isfield(STATE_REF,'mask')
    mask        =  logical(STATE_REF.mask);
else
    mask        =  true(ny,nx,size(phi_ref_all,3));
end

Np              =  size(phi_ref_all,3);
LKvec_full      =  LK(:);

tol_direct      =  40000;
tol_iter        =  1e-10;
maxit           =  500;

phi_new_all     =  zeros(size(phi_ref_all));

for ip = 1:Np

    mask_p      =  mask(:,:,ip);
    act_lin     =  find(mask_p);
    Nact        =  numel(act_lin);

    if Nact == 0
        phi_new_all(:,:,ip) = 0;
        continue
    end

    phi_old     =  phi_old_all(:,:,ip);
    phi_ref     =  phi_ref_all(:,:,ip);

    phi_old_vec =  phi_old(:);
    phi_ref_vec =  phi_ref(:);
    s_vec       =  s{ip}(:);

    if ndims(A_ac) == 2
        A_ip    =  A_ac;
    else
        A_ip    =  A_ac(:,:,ip);
    end

    A_vec       =  A_ip(:);

    [ia,ja]     =  ind2sub([ny,nx],act_lin);

    hasL        =  ja > 1;
    hasR        =  ja < nx;
    hasU        =  ia > 1;
    hasD        =  ia < ny;

    indL        =  zeros(Nact,1);
    indR        =  zeros(Nact,1);
    indU        =  zeros(Nact,1);
    indD        =  zeros(Nact,1);

    indL(hasL)  =  sub2ind([ny,nx],ia(hasL)  ,ja(hasL)-1);
    indR(hasR)  =  sub2ind([ny,nx],ia(hasR)  ,ja(hasR)+1);
    indU(hasU)  =  sub2ind([ny,nx],ia(hasU)-1,ja(hasU)  );
    indD(hasD)  =  sub2ind([ny,nx],ia(hasD)+1,ja(hasD)  );

    g2c         =  zeros(ny*nx,1,'int32');
    g2c(act_lin)=  int32(1:Nact);

    comp        =  (1:Nact).';

    LKc         =  LKvec_full(act_lin);

    coefL       =  zeros(Nact,1);
    coefR       =  zeros(Nact,1);
    coefU       =  zeros(Nact,1);
    coefD       =  zeros(Nact,1);

    coefL(hasL) =  -(LKvec_full(indL(hasL)) + LKc(hasL))/2/dx^2;
    coefR(hasR) =  -(LKvec_full(indR(hasR)) + LKc(hasR))/2/dx^2;
    coefU(hasU) =  -(LKvec_full(indU(hasU)) + LKc(hasU))/2/dy^2;
    coefD(hasD) =  -(LKvec_full(indD(hasD)) + LKc(hasD))/2/dy^2;

    coefC       =  -(coefL + coefR + coefU + coefD) ...
                   + 1/dt_phy ...
                   + A_vec(act_lin);

    RHS         =  phi_old_vec(act_lin)/dt_phy ...
                   + s_vec(act_lin) ...
                   + A_vec(act_lin).*phi_ref_vec(act_lin);

    est_nnz     =  5*Nact;
    rows        =  zeros(est_nnz,1);
    cols        =  zeros(est_nnz,1);
    vals        =  zeros(est_nnz,1);
    ptr         =  1;

    rows(ptr:ptr+Nact-1) = comp;
    cols(ptr:ptr+Nact-1) = comp;
    vals(ptr:ptr+Nact-1) = coefC;
    ptr = ptr + Nact;

    neigh       =  g2c(indL(hasL));
    rows0       =  comp(hasL);
    vals0       =  coefL(hasL);
    ok          =  neigh > 0;
    nadd        =  nnz(ok);
    if nadd > 0
        rows(ptr:ptr+nadd-1) = rows0(ok);
        cols(ptr:ptr+nadd-1) = double(neigh(ok));
        vals(ptr:ptr+nadd-1) = vals0(ok);
        ptr = ptr + nadd;
    end

    neigh       =  g2c(indR(hasR));
    rows0       =  comp(hasR);
    vals0       =  coefR(hasR);
    ok          =  neigh > 0;
    nadd        =  nnz(ok);
    if nadd > 0
        rows(ptr:ptr+nadd-1) = rows0(ok);
        cols(ptr:ptr+nadd-1) = double(neigh(ok));
        vals(ptr:ptr+nadd-1) = vals0(ok);
        ptr = ptr + nadd;
    end

    neigh       =  g2c(indU(hasU));
    rows0       =  comp(hasU);
    vals0       =  coefU(hasU);
    ok          =  neigh > 0;
    nadd        =  nnz(ok);
    if nadd > 0
        rows(ptr:ptr+nadd-1) = rows0(ok);
        cols(ptr:ptr+nadd-1) = double(neigh(ok));
        vals(ptr:ptr+nadd-1) = vals0(ok);
        ptr = ptr + nadd;
    end

    neigh       =  g2c(indD(hasD));
    rows0       =  comp(hasD);
    vals0       =  coefD(hasD);
    ok          =  neigh > 0;
    nadd        =  nnz(ok);
    if nadd > 0
        rows(ptr:ptr+nadd-1) = rows0(ok);
        cols(ptr:ptr+nadd-1) = double(neigh(ok));
        vals(ptr:ptr+nadd-1) = vals0(ok);
        ptr = ptr + nadd;
    end

    rows        =  rows(1:ptr-1);
    cols        =  cols(1:ptr-1);
    vals        =  vals(1:ptr-1);

    Lmat        =  sparse(rows,cols,vals,Nact,Nact);

    x0          =  phi_ref_vec(act_lin);

    if Nact <= tol_direct

        x       =  Lmat \ RHS;

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

    phi_new          =  zeros(ny*nx,1);
    phi_new(act_lin) =  x;

    tmp             =  reshape(phi_new,ny,nx);
    tmp(~mask_p)    =  0;

    phi_new_all(:,:,ip) = tmp;

end

STATE_NEW       =  STATE_REF;
STATE_NEW.phi   =  phi_new_all;

if Norm == 1
    STATE_NEW.phi = Norm_Phi(max(STATE_NEW.phi,0));
end

STATE_NEW.p     =  Calc_p(MODEL,STATE_NEW.phi);

end