function DIAG = Diag_dE_Source(STATE,STATE_RAW,STATE_CORR,STATE_TRIAL)
%DIAG_DE_SOURCE Diagnose where E spikes are created in one timestep.
%
% Call after STATE_RAW, STATE_CORR, STATE_TRIAL have been computed.

DIAG.raw   = MaxCellDiff(STATE_RAW.E,STATE.E);
DIAG.corr  = MaxCellDiff(STATE_CORR.E,STATE_RAW.E);
DIAG.trial = MaxCellDiff(STATE_TRIAL.E,STATE_CORR.E);
DIAG.total = MaxCellDiff(STATE_TRIAL.E,STATE.E);

fprintf('dE source: raw %.3e, corr %.3e, trial %.3e, total %.3e\n', ...
    DIAG.raw,DIAG.corr,DIAG.trial,DIAG.total);

end

function d = MaxCellDiff(A,B)

d = 0;

for i = 1:numel(A)
    d = max(d,max(abs(A{i}(:)-B{i}(:))));
end

end
