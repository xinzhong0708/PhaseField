function [phi] = Norm_Phi(phi)
phi(phi<0) = 0;
phi(phi>1) = 1;
phi        = phi ./ sum(phi,3);
end