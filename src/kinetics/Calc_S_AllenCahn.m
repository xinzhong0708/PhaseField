function STATE = Calc_S_AllenCahn(STATE,PARAM,MODEL)
%CALC_S_ALLENCAHN Allen-Cahn source term using structured variables.
%
% If anisotropic coefficients exist:
%
%   PARAM.L_AC(:,:,alpha)   for chemical driving force
%   PARAM.Lm_AC(:,:,alpha)  for double-well/barrier term
%
% The Kundin stiffness scales the interface term only. It should not scale
% the chemical driving force.

phi = STATE.phi;
omg = STATE.omg;

[ny,nx,Np] = size(phi);

if isfield(PARAM,'eps_phi')
    eps_phi = PARAM.eps_phi;
else
    eps_phi = 1e-14;
end

eps_div = 1e-30;

% ------------------------------------------------------------
% Interpolation denominator
% ------------------------------------------------------------
D = sum(phi.^2,3) + eps_phi;

% weighted omega average: omg_bar = sum_beta phi_beta^2 omega_beta / D
omg_bar = sum(phi.^2 .* omg,3) ./ D;

s = cell(1,Np);

% ------------------------------------------------------------
% Loop through phases
% ------------------------------------------------------------
for alp = 1:Np

    phi_a = phi(:,:,alp);

    % chemical part: Q_a = d/dphi_a sum_i p_i omega_i
    Q_a = 2 .* phi_a ./ D .* (omg(:,:,alp) - omg_bar);

    % --------------------------------------------------------
    % Chemical mobility: unscaled L
    % --------------------------------------------------------
    if isfield(PARAM,'L_AC') && ~isempty(PARAM.L_AC)
        LL = Get_Field(PARAM.L_AC,alp,ny,nx);
    else
        LL = Get_Field(PARAM.L,alp,ny,nx);
    end

    % --------------------------------------------------------
    % Barrier term
    % --------------------------------------------------------
    if isfield(PARAM,'Lm_AC') && ~isempty(PARAM.Lm_AC) && size(PARAM.Lm_AC,3) >= alp

        Lm = Get_Field(PARAM.Lm_AC,alp,ny,nx);

        % MODEL.dgdphi already contains m.
        % Therefore divide by base m before multiplying by Lm = L*m*S.
        if isfield(PARAM,'m_AC_base') && ~isempty(PARAM.m_AC_base)
            m0 = Get_Field(PARAM.m_AC_base,alp,ny,nx);

        elseif isfield(PARAM,'Lm') && isfield(PARAM,'L')
            Lm0 = Get_Field(PARAM.Lm,alp,ny,nx);
            L0  = Get_Field(PARAM.L,alp,ny,nx);
            m0  = Lm0 ./ max(L0,eps_div);

        else
            m0 = ones(ny,nx);
        end

        dW = MODEL.dgdphi(phi_a) ./ max(m0,eps_div);

        s{alp} = -Lm .* dW - LL .* Q_a;

    else

        % Old isotropic form
        s{alp} = -LL .* MODEL.dgdphi(phi_a) - LL .* Q_a;

    end

end

STATE.S_AC = s;

end


%% HELPERS
% ------------------------------------------------------------

function F = Get_Field(A,ip,ny,nx)

if isscalar(A)

    F = A*ones(ny,nx);

elseif ndims(A) == 2 && isequal(size(A),[ny nx])

    F = A;

elseif ndims(A) == 3

    F = A(:,:,ip);

elseif isvector(A)

    F = A(min(ip,numel(A))) * ones(ny,nx);

else

    F = A(1)*ones(ny,nx);

end

end