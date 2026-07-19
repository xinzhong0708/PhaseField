%Laplacian operator
function [Lap] = Laplacian_Periodic(M,A,dx,dy)
Ax  = [A(:,end),A,A(:,1)];
Ay  = [A(end,:);A;A(1,:)];
qx  = M.*diff(diff(Ax,1,2),1,2)/dx/dx;
qy  = M.*diff(diff(Ay,1,1),1,1)/dy/dy;
Lap = qx+qy;
end