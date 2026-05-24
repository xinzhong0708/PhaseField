function STATE = Calc_S_AllenCahn(STATE,PARAM,MODEL)
%CALC_S_ALLENCAHN Allen-Cahn source term using structured variables.
%
% Fast version for:
%   p_i = phi_i^2 / sum_j(phi_j^2)
%
% Avoids the double loop over alpha and beta.

phi        = STATE.phi;
omg        = STATE.omg;
LL         = PARAM.L;
[ny,nx,Np] = size(phi);

if isfield(PARAM,'eps_phi')
    eps_phi = PARAM.eps_phi;
else
    eps_phi = 1e-14;
end

D = sum(phi.^2,3) + eps_phi;

% weighted omega average:
%   omg_bar = sum_beta phi_beta^2 * omega_beta / D
omg_bar = sum(phi.^2 .* omg,3) ./ D;

s = cell(1,Np);

for alp = 1:Np
    phi_a  = phi(:,:,alp);
    % chemical part:
    Q_a    = 2 .* phi_a ./ D .* (omg(:,:,alp) - omg_bar);
    s{alp} = -LL .* MODEL.dgdphi(phi_a) -LL .* Q_a;
end
STATE.S_AC = s;

end