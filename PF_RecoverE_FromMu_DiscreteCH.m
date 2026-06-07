function STATE_NEW = PF_RecoverE_FromMu_DiscreteCH(STATE_TIME,STATE_NEW,PARAM,GRID,PHYS,NUM)
%PF_RECOVERE_FROMMU_DISCRETECH Recover conserved E from discrete CH balance.
%
% This replaces the local tangent reconstruction
%
%     E_new = E_old + chi*dmu + e*dp
%
% by the conservative discrete CH equation
%
%     A_E * E_new + D * mu_new = E_time / dt
%
% where A_E and D are the same operators used in PF_Coupled_ACCH_LETangent
% and PF_CH_LECorrector_FixedP_Band.
%
% This does not change phi, p, or mu_e. It only replaces STATE_NEW.E.

E_time = STATE_TIME.E;
mu_new = STATE_NEW.mu_e;

[ny,nx,~] = size(STATE_NEW.p);
Ne        = numel(E_time);
Nnode     = ny*nx;

dt        = NUM.dt_phy;
dx        = GRID.dx;
dy        = GRID.dy;

dx2       = dx^2;
dy2       = dy^2;
dx4       = dx2^2;
dy4       = dy2^2;

% Spatially varying kappa
if isfield(PARAM,'kappa_eff') && ~isempty(PARAM.kappa_eff)
    kappa_eff = PARAM.kappa_eff;
    if isscalar(kappa_eff)
        kappa_eff = kappa_eff*ones(ny,nx);
    end
else
    kappa_eff = PHYS.kappa*ones(ny,nx);
end

% Reflective neighbour indices
[Igrid,Jgrid] = ndgrid(1:ny,1:nx);

ii       = Igrid(:);
jj       = Jgrid(:);

refI     = @(i,sh) reflect_index_local(i+sh,ny);
refJ     = @(j,sh) reflect_index_local(j+sh,nx);

jjL      = refJ(jj,-1);
jjR      = refJ(jj,+1);
iiU      = refI(ii,-1);
iiD      = refI(ii,+1);

jjL2     = refJ(jj,-2);
jjR2     = refJ(jj,+2);
iiU2     = refI(ii,-2);
iiD2     = refI(ii,+2);

iiUR     = refI(ii,-1);
jjUR     = refJ(jj,+1);
iiDR     = refI(ii,+1);
jjDR     = refJ(jj,+1);
iiUL     = refI(ii,-1);
jjUL     = refJ(jj,-1);
iiDL     = refI(ii,+1);
jjDL     = refJ(jj,-1);

idx_c    = sub2ind([ny,nx],ii,jj);
idx_L    = sub2ind([ny,nx],ii,jjL);
idx_R    = sub2ind([ny,nx],ii,jjR);
idx_U    = sub2ind([ny,nx],iiU,jj);
idx_D    = sub2ind([ny,nx],iiD,jj);

idx_L2   = sub2ind([ny,nx],ii,jjL2);
idx_R2   = sub2ind([ny,nx],ii,jjR2);
idx_U2   = sub2ind([ny,nx],iiU2,jj);
idx_D2   = sub2ind([ny,nx],iiD2,jj);

idx_UR   = sub2ind([ny,nx],iiUR,jjUR);
idx_DR   = sub2ind([ny,nx],iiDR,jjDR);
idx_UL   = sub2ind([ny,nx],iiUL,jjUL);
idx_DL   = sub2ind([ny,nx],iiDL,jjDL);

% Recover each conserved elemental field independently
E_out = cell(size(E_time));

for ie = 1:Ne

    Ml = PARAM.M{ie};

    M_c = Ml(idx_c);
    M_L = Ml(idx_L);
    M_R = Ml(idx_R);
    M_U = Ml(idx_U);
    M_D = Ml(idx_D);

    % Diffusion operator D on mu
    d_L = -(M_L + M_c)/2/dx2;
    d_R = -(M_R + M_c)/2/dx2;
    d_U = -(M_U + M_c)/2/dy2;
    d_D = -(M_D + M_c)/2/dy2;
    d_C = -(d_L + d_R + d_U + d_D);

    % A_E operator on E: 1/dt plus optional fourth-order kappa term
    kappa_c = kappa_eff(idx_c);
    MK_c    = M_c .* kappa_c;

    q_L  = MK_c .* (-4/dx4 - 4/(dx2*dy2));
    q_R  = MK_c .* (-4/dx4 - 4/(dx2*dy2));
    q_U  = MK_c .* (-4/dy4 - 4/(dx2*dy2));
    q_D  = MK_c .* (-4/dy4 - 4/(dx2*dy2));

    q_L2 = MK_c .* (1/dx4);
    q_R2 = MK_c .* (1/dx4);
    q_U2 = MK_c .* (1/dy4);
    q_D2 = MK_c .* (1/dy4);

    q_UR = MK_c .* (2/(dx2*dy2));
    q_DR = MK_c .* (2/(dx2*dy2));
    q_UL = MK_c .* (2/(dx2*dy2));
    q_DL = MK_c .* (2/(dx2*dy2));

    q_C  = MK_c .* (6/dx4 + 6/dy4 + 8/(dx2*dy2));

    a_C  = 1/dt + q_C;
    a_L  = q_L;
    a_R  = q_R;
    a_U  = q_U;
    a_D  = q_D;

    a_L2 = q_L2;
    a_R2 = q_R2;
    a_U2 = q_U2;
    a_D2 = q_D2;

    a_UR = q_UR;
    a_DR = q_DR;
    a_UL = q_UL;
    a_DL = q_DL;

    % D * mu_new
    mu = mu_new{ie};

    mu_c = mu(idx_c);
    mu_L = mu(idx_L);
    mu_R = mu(idx_R);
    mu_U = mu(idx_U);
    mu_D = mu(idx_D);

    D_mu = d_C.*mu_c + d_L.*mu_L + d_R.*mu_R + d_U.*mu_U + d_D.*mu_D;

    % RHS: E_time/dt - D*mu_new
    Et = E_time{ie};
    rhs = Et(idx_c)/dt - D_mu;

    % Assemble A_E matrix
    max_nnz = Nnode*13 + 100;
    rows = zeros(max_nnz,1);
    cols = zeros(max_nnz,1);
    vals = zeros(max_nnz,1);
    k    = 1;

    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_c, a_C);
    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_L, a_L);
    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_R, a_R);
    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_U, a_U);
    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_D, a_D);

    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_L2,a_L2);
    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_R2,a_R2);
    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_U2,a_U2);
    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_D2,a_D2);

    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_UR,a_UR);
    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_DR,a_DR);
    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_UL,a_UL);
    [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,idx_c,idx_DL,a_DL);

    rows = rows(1:k-1);
    cols = cols(1:k-1);
    vals = vals(1:k-1);

    AE = sparse(rows,cols,vals,Nnode,Nnode);

    Evec = AE\rhs;

    E_out{ie} = reshape(Evec,ny,nx);

end

STATE_NEW.E = E_out;

end


function [rows,cols,vals,k] = add_block_local(rows,cols,vals,k,r,c,v)

n = numel(r);

if k+n-1 > numel(rows)
    grow = max(numel(rows),n+1000);
    rows = [rows; zeros(grow,1)];
    cols = [cols; zeros(grow,1)];
    vals = [vals; zeros(grow,1)];
end

rows(k:k+n-1) = r(:);
cols(k:k+n-1) = c(:);
vals(k:k+n-1) = v(:);

k = k+n;

end


function idx = reflect_index_local(idx,n)

idx(idx < 1) = 2 - idx(idx < 1);
idx(idx > n) = 2*n - idx(idx > n);

idx = max(min(idx,n),1);

end