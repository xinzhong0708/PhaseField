alpha = 1;
dx    = GRID.dx;
dy    = GRID.dy;

if isfield(PARAM,'LK_AC') && size(PARAM.LK_AC,3) >= alpha
    LK_a = PARAM.LK_AC(:,:,alpha);
else
    LK_a = PARAM.LK;
end

STATE_SRC = Calc_S_AllenCahn(STATE_OLD,PARAM,MODEL);
S_old     = STATE_SRC.S_AC;

G_cap  = LK_a .* Laplacian_Reflex(STATE_OLD.phi(:,:,alpha),dx,dy);
G_src  = S_old{alpha};

figure; clf
subplot(131); pcolor(G_cap); shading interp; colorbar; title('capillary term')
subplot(132); pcolor(G_src); shading interp; colorbar; title('source term')
subplot(133); pcolor(abs(G_src)./max(abs(G_cap),eps)); shading interp; colorbar; title('|source| / |capillary|')

fprintf('max capillary = %.6e\n',max(abs(G_cap(:))));
fprintf('max source    = %.6e\n',max(abs(G_src(:))));
fprintf('ratio source/capillary = %.6e\n',max(abs(G_src(:)))/max(max(abs(G_cap(:))),eps));







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
