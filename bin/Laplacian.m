%Laplacian operator
function [Lap] = Laplacian(A,dx,dy)
qx  = diff(diff(A,1,2),1,2)/dx/dx;
qy  = diff(diff(A,1,1),1,1)/dy/dy;
Lap = expand2(qx(2:end-1,:)+qy(:,2:end-1),5);
end
