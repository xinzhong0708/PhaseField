function STATE = Calc_S_AllenCahn(STATE,PARAM,MODEL)
%CALC_S_ALLENCAHN Allen-Cahn source term using structured variables.
%
% Fast version for:
%   p_i = phi_i^2 / sum_j(phi_j^2)
%
% Optional anisotropy:
%   PARAM.L_AC(:,:,alpha)
%
% IMPORTANT:
% MODEL.dgdphi already contains PHYS.m in the current code.
% Therefore the barrier term must be multiplied by L, not by L*m.

phi        = STATE.phi;
omg        = STATE.omg;
Lchem      = PARAM.L;
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

    phi_a = phi(:,:,alp);

    if isfield(PARAM,'L_AC') && ~isempty(PARAM.L_AC) && size(PARAM.L_AC,3) >= alp
        L_a = PARAM.L_AC(:,:,alp);
    else
        L_a = PARAM.L;
    end

    % Chemical part:
    %   d/dphi_alpha sum_beta p_beta omega_beta
    Q_a = 2 .* phi_a ./ D .* (omg(:,:,alp) - omg_bar);

    % Barrier/source part uses anisotropic L_a.
    % Chemical omega driving still uses PARAM.L.
    s{alp} = -L_a .* MODEL.dgdphi(phi_a) - Lchem .* Q_a;

end

STATE.S_AC = s;

end