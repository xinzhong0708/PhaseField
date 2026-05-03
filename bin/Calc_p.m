function [p] = Calc_p(F,phi)
%Interpolation function
for i = 1:size(phi,3)
    p(:,:,i)     =  real(F.p_fun(i,phi));
end
end