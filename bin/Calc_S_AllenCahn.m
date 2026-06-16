function STATE = Calc_S_AllenCahn(STATE,PARAM,MODEL)
%CALC_S_ALLENCAHN Allen-Cahn source term using structured variables.
%
% Backward compatible with the old isotropic form, but if anisotropic
% coefficients exist it uses
%
%   PARAM.L_AC(:,:,alpha)   for chemical driving force
%   PARAM.Lm_AC(:,:,alpha)  for double-well/barrier term
%
% This is needed for the Kundin faceted model, where stiffness scales the
% interface terms but should not scale the chemical driving force.

phi = STATE.phi;
omg = STATE.omg;

[ny,nx,Np] = size(phi); %#ok<ASGLU>

if isfield(PARAM,'eps_phi')
    eps_phi = PARAM.eps_phi;
else
    eps_phi = 1e-14;
end

D = sum(phi.^2,3) + eps_phi;

% weighted omega average: omg_bar = sum_beta phi_beta^2 omega_beta / D
omg_bar = sum(phi.^2 .* omg,3) ./ D;

s = cell(1,Np);

for alp = 1:Np

    phi_a = phi(:,:,alp);

    % chemical part: Q_a = d/dphi_a sum_i p_i omega_i
    Q_a = 2 .* phi_a ./ D .* (omg(:,:,alp) - omg_bar);

    if isfield(PARAM,'L_AC') && ndims(PARAM.L_AC) == 3 && size(PARAM.L_AC,3) >= alp
        LL = PARAM.L_AC(:,:,alp);
    else
        LL = PARAM.L;
    end

    if isfield(PARAM,'Lm_AC') && ndims(PARAM.Lm_AC) == 3 && size(PARAM.Lm_AC,3) >= alp

        Lm = PARAM.Lm_AC(:,:,alp);

        if isfield(PARAM,'Lm') && isfield(PARAM,'L')
            m0 = PARAM.Lm ./ max(PARAM.L,eps);
            dW = MODEL.dgdphi(phi_a) ./ max(m0,eps);
        else
            dW = MODEL.dgdphi(phi_a);
        end

        s{alp} = -Lm .* dW - LL .* Q_a;

    else

        s{alp} = -LL .* MODEL.dgdphi(phi_a) - LL .* Q_a;

    end

end

STATE.S_AC = s;

end
