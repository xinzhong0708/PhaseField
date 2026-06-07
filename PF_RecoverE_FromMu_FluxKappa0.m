function STATE_NEW = PF_RecoverE_FromMu_FluxKappa0(STATE_TIME,STATE_NEW,PARAM,GRID,NUM)
%PF_RECOVERE_FROMMU_FLUXKAPPA0 Recover E from conservative CH flux.
%
% This is only for kappa = 0 tests.
%
% The discrete CH residual used by the corrector is:
%
%     E_TIME/dt - A_E*E_NEW - D*mu_NEW = 0
%
% For kappa = 0:
%
%     A_E = I/dt
%
% therefore:
%
%     E_NEW = E_TIME - dt * D*mu_NEW
%
% This recovers E from the PDE flux, not from local tangent:
%
%     E_NEW ~= E_OLD + chi*dmu + e*dp
%
% Use this only as a test to see whether LE kinks come from tangent
% reconstruction of E.

E_time = STATE_TIME.E;
mu_new = STATE_NEW.mu_e;

[ny,nx,~] = size(STATE_NEW.p);
Ne        = numel(E_time);

dt = NUM.dt_phy;
dx = GRID.dx;
dy = GRID.dy;

dx2 = dx^2;
dy2 = dy^2;

% Reflective neighbor indices
[Igrid,Jgrid] = ndgrid(1:ny,1:nx);

ii = Igrid(:);
jj = Jgrid(:);

jjL = reflect_index_local(jj-1,nx);
jjR = reflect_index_local(jj+1,nx);
iiU = reflect_index_local(ii-1,ny);
iiD = reflect_index_local(ii+1,ny);

idx_c = sub2ind([ny,nx],ii,jj);
idx_L = sub2ind([ny,nx],ii,jjL);
idx_R = sub2ind([ny,nx],ii,jjR);
idx_U = sub2ind([ny,nx],iiU,jj);
idx_D = sub2ind([ny,nx],iiD,jj);

E_out = cell(size(E_time));

for ie = 1:Ne

    Ml = PARAM.M{ie};

    M_c = Ml(idx_c);
    M_L = Ml(idx_L);
    M_R = Ml(idx_R);
    M_U = Ml(idx_U);
    M_D = Ml(idx_D);

    % D operator copied from the CH/corrector sign convention.
    % For constant M:
    %   D*mu = -M*laplacian(mu)
    d_L = -(M_L + M_c)/2/dx2;
    d_R = -(M_R + M_c)/2/dx2;
    d_U = -(M_U + M_c)/2/dy2;
    d_D = -(M_D + M_c)/2/dy2;
    d_C = -(d_L + d_R + d_U + d_D);

    mu = mu_new{ie};

    Dmu = d_C.*mu(idx_c) + ...
          d_L.*mu(idx_L) + ...
          d_R.*mu(idx_R) + ...
          d_U.*mu(idx_U) + ...
          d_D.*mu(idx_D);

    Et = E_time{ie};

    Evec = Et(idx_c) - dt.*Dmu;

    E_out{ie} = reshape(Evec,ny,nx);

end

STATE_NEW.E = E_out;

end


function idx = reflect_index_local(idx,n)

idx(idx < 1) = 2 - idx(idx < 1);
idx(idx > n) = 2*n - idx(idx > n);

idx = max(min(idx,n),1);

end