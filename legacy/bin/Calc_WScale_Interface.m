function lambda = Calc_WScale_Interface(p_phase,p0,p1)
%CALC_WSCALE_INTERFACE Smooth scaling factor for excess/Margules energy.
%
% lambda = 0 for p <= p0
% lambda = 1 for p >= p1
% smoothstep transition in between.
%
% This is used to turn off/soften spinodal excess energy in diffuse
% interface regions, while preserving full thermodynamics in bulk phase.

if p1 <= p0
    error('p1 must be larger than p0.');
end

x = (p_phase - p0) ./ (p1 - p0);
x = min(max(x,0),1);

lambda = x.^2 .* (3 - 2*x);

end