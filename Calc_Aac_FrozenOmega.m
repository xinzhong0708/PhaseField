function A_ac = Calc_Aac_FrozenOmega(STATE,PARAM,MODEL,fac,eps_fd,Amin,Amax)
%CALC_AAC_FROZENOMEGA
% Estimate local Allen-Cahn source stiffness for stabilized semi-implicit AC.
%
% A_ac approximates/bounds:
%   | d/dphi [ L*(-m*dg/dphi - d/dphi(sum p*omega)) ] |
%
% omega is frozen at the old LE state.
%
% INPUT
%   fac    : safety factor, e.g. 2 to 10
%   eps_fd : finite-difference perturbation, e.g. 1e-6
%   Amin   : minimum A_ac
%   Amax   : maximum A_ac, use [] if not wanted
%
% OUTPUT
%   A_ac   : ny x nx stabilizer field

phi = STATE.phi;
omg = STATE.omg;

[ny,nx,Np] = size(phi);

% ------------------------------------------------------------
% Use omega contrast only, because common omega offset has no force
% ------------------------------------------------------------
omg_c = omg - mean(omg,3);

% ------------------------------------------------------------
% Double-well stiffness
% For g = phi^2*(1-phi)^2, max |g''| = 2 in [0,1]
% source stiffness = L*m*|g''|
% ------------------------------------------------------------
A_dw = 2 * PARAM.Lm;

% ------------------------------------------------------------
% Chemical stiffness from frozen omega and p(phi)
% Q_alpha = sum_beta omega_beta * dp_beta/dphi_alpha
% estimate dQ_alpha/dphi_gamma by finite difference
% ------------------------------------------------------------
A_chem = zeros(ny,nx);

Q0 = Local_Q(phi,omg_c,MODEL);

for ig = 1:Np

    phi_p = phi;
    phi_p(:,:,ig) = phi_p(:,:,ig) + eps_fd;

    Qp = Local_Q(phi_p,omg_c,MODEL);

    dQ = (Qp - Q0)/eps_fd;   % ny x nx x Np, derivative wrt phi_gamma

    % row-sum style bound:
    % for each alpha, add |dQ_alpha/dphi_gamma|
    A_chem = A_chem + max(abs(dQ),[],3);

end

% Chemical source has factor L
A_chem = PARAM.L .* A_chem;

% ------------------------------------------------------------
% Total stabilizer
% ------------------------------------------------------------
A_ac = fac * (A_dw + A_chem);

% Floor and cap
if nargin >= 6 && ~isempty(Amin)
    A_ac = max(A_ac,Amin);
end

if nargin >= 7 && ~isempty(Amax)
    A_ac = min(A_ac,Amax);
end

end

function Q = Local_Q(phi,omg_c,MODEL)
% Q_alpha = sum_beta omega_beta * dp_beta/dphi_alpha

[ny,nx,Np] = size(phi);
Q = zeros(ny,nx,Np);

for alpha = 1:Np
    tmp = zeros(ny,nx);

    for beta = 1:Np
        tmp = tmp + omg_c(:,:,beta) .* MODEL.dpdphi(alpha,beta,phi);
    end

    Q(:,:,alpha) = tmp;
end

end