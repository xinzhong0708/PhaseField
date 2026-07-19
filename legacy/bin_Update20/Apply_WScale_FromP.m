function pars_out = Apply_WScale_FromP(pars,p,p0,p1,gamma)
%APPLY_WSCALE_FROMP Add p-dependent w_scale to every pars{ip}.
% p is ny x nx x Np.
% pars is pars{1} ... pars{Np}.
if nargin < 5 || isempty(gamma)
    gamma = 1;
end
pars_out = pars;
Np       = numel(pars);
for ip = 1:Np
    lambda = Calc_WScale_Interface(p(:,:,ip),p0,p1);
    % Stronger suppression
    lambda = lambda.^gamma;
    pars_out{ip}.w_scale = lambda(:).';
end
end