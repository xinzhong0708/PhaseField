function STATE_OUT = PF_RecoverE_MatchMu_LE(STATE_RAW,STATE_LE,NUM)
%PF_RECOVERE_MATCHMU_LE Correct E so fixed-E LE matches smooth ACCH mu.
%
% ACCH gives a smooth PDE chemical potential STATE_RAW.mu_e.
% Fixed-E LE gives the chemical potential implied by STATE_RAW.E.
%
% If the LE mu has a kink, correct E by:
%
%   dE = chi_LE * (mu_RAW - mu_LE)
%
% Then rerun LE on the corrected E.
%
% This keeps the final closure fixed-E LE and avoids accepting GP.

STATE_OUT = STATE_RAW;

Ne = numel(STATE_RAW.E);

% Strength of correction
theta = 0.5;
if isfield(NUM,'E_mu_match_theta')
    theta = NUM.E_mu_match_theta;
end

% Optional cap on |dE| per correction
dE_cap = inf;
if isfield(NUM,'E_mu_match_cap')
    dE_cap = NUM.E_mu_match_cap;
end

% Preserve mean E from STATE_RAW
preserve_mean = 1;
if isfield(NUM,'E_mu_match_preserve_mean')
    preserve_mean = NUM.E_mu_match_preserve_mean;
end

% Chemical-potential mismatch
dmu = cell(1,Ne);

for ie = 1:Ne
    dmu{ie} = STATE_RAW.mu_e{ie} - STATE_LE.mu_e{ie};
end

% E correction
dE_all = cell(1,Ne);

for ie = 1:Ne

    dE = zeros(size(STATE_RAW.E{ie}));

    for je = 1:Ne
        dE = dE + STATE_LE.chi{ie,je}.*dmu{je};
    end

    dE = theta.*dE;

    if isfinite(dE_cap)
        dE = max(dE,-dE_cap);
        dE = min(dE, dE_cap);
    end

    if preserve_mean == 1
        dE = dE - mean(dE(:));
    end

    dE_all{ie} = dE;
end

for ie = 1:Ne
    STATE_OUT.E{ie} = STATE_RAW.E{ie} + dE_all{ie};
end

end