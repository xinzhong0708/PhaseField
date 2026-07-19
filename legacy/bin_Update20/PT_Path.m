function [T,P] = PT_Path(t,PT)

%PT_PATH Interpolate P-T path at physical time t.

T = interp1(PT.t_path,PT.T_path,t,'linear','extrap');
P = interp1(PT.t_path,PT.P_path,t,'linear','extrap');

end