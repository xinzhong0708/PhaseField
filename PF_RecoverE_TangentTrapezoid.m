function STATE_OUT = PF_RecoverE_TangentTrapezoid(STATE_OLD,STATE_RAW,STATE_COEF,NUM)
%PF_RECOVERE_TANGENTTRAPEZOID Improve E recovery after ACCH predictor.
%
% Replaces first-order tangent recovery:
%
%   E_new = E_old + chi_old*dmu + e_old*dp
%
% by a local trapezoidal tangent recovery:
%
%   E_new = E_old
%         + 0.5*(chi_old + chi_new)*dmu
%         + 0.5*(e_old   + e_new)  *dp
%
% STATE_COEF is usually STATE_LE0, i.e. local thermodynamic response
% evaluated from STATE_RAW. This keeps the correction local and consistent
% with the same tangent form used by the ACCH solver.
%
% It does not solve an extra CH equation and does not overwrite E using
% a different global projection.

STATE_OUT = STATE_RAW;

Ne     = numel(STATE_OLD.E);
Ngrain = size(STATE_OLD.p,3);

% dmu from the ACCH predictor
dmu = cell(1,Ne);
for ie = 1:Ne
    dmu{ie} = STATE_RAW.mu_e{ie} - STATE_OLD.mu_e{ie};
end

% dp from the accepted raw phase update
dp = STATE_RAW.p - STATE_OLD.p;

% Trapezoidal recovery
E_new = STATE_OLD.E;

for ie = 1:Ne

    En = STATE_OLD.E{ie};

    % Chemical part: average old and predicted chi
    for je = 1:Ne
        chi_avg = 0.5*(STATE_OLD.chi{ie,je} + STATE_COEF.chi{ie,je});
        En = En + chi_avg.*dmu{je};
    end

    % Phase-change part: average old and predicted e
    for ig = 1:Ngrain
        e_avg = 0.5*(STATE_OLD.e{ig}{ie} + STATE_COEF.e{ig}{ie});
        En = En + e_avg.*dp(:,:,ig);
    end

    E_new{ie} = En;
end

STATE_OUT.E = E_new;

if nargin >= 4 && isfield(NUM,'norm_E') && NUM.norm_E == 1
    STATE_OUT.E = EnforceMeanE_Local(STATE_OUT.E,STATE_OLD.E);
end

end