function [STATE,NUM,DIAG] = Update_TimeStep(STATE,STATE_T,PARAM,MODEL,NUM)
%UPDATE_TIMESTEP Accept/reject trial step and update next timestep.
%
% Key behavior:
%   - reject if dE or dp exceeds target
%   - dmu is only a soft limiter, except for extreme jumps
%   - accepted steps do not shrink dt
%   - dt grows only after several very good steps

%======================================================================
% Final LE correction BEFORE diagnostics
%======================================================================
STATE_T = LE_Run(STATE_T,PARAM,MODEL);

%======================================================================
% Diagnostics on fully corrected trial state
%======================================================================
dE      = Max_Cell_Diff(STATE_T.E,STATE.E);
dmu     = Max_Cell_Diff(STATE_T.mu_e,STATE.mu_e);
dp      = max(abs(STATE_T.p(:) - STATE.p(:)));

% Error ratios
rE      = dE  / max(NUM.dE_target,eps);
rp      = dp  / max(NUM.dp_target,eps);
rmu     = dmu / max(NUM.dmu_target,eps);

% Hard and soft errors
err_hard = max([rE,rp]);

% dmu affects timestep growth, but weakly
err_eff  = max(err_hard,0.2*rmu);

%======================================================================
% Accept/reject
%======================================================================
reject = false;

% Hard reject only from primary fields
if rE > 1 || rp > 1
    reject = true;
end

% Optional extreme mu guard
if rmu > 20
    reject = true;
end

%======================================================================
% Timestep update with hysteresis
%======================================================================
scale = 1.0;
good_count_trial = NUM.dt_good_count;

if reject

    % Shrink only on rejection
    scale = NUM.dt_shrink_fac;
    dt_next = max(NUM.dt_phy * scale, NUM.dt_min);
    good_count_trial = 0;

else

    % Accepted step.
    % Do not shrink dt even if err_eff is close to 1.
    % Only grow after several very safe steps.
    if err_eff < NUM.err_grow
        good_count_trial = NUM.dt_good_count + 1;
    else
        good_count_trial = 0;
    end

    if good_count_trial >= NUM.dt_grow_after
        scale = NUM.dt_grow_fac;
        dt_next = min(NUM.dt_phy * scale, NUM.dt_max);
        good_count_trial = 0;
    else
        scale = 1.0;
        dt_next = min(NUM.dt_phy, NUM.dt_max);
    end

end

%======================================================================
% Diagnostics
%======================================================================
DIAG.dE       = dE;
DIAG.dmu      = dmu;
DIAG.dp       = dp;
DIAG.rE       = rE;
DIAG.rp       = rp;
DIAG.rmu      = rmu;
DIAG.err_hard = err_hard;
DIAG.err_eff  = err_eff;
DIAG.scale    = scale;
DIAG.dt_try   = NUM.dt_phy;
DIAG.dt_next  = dt_next;

[~,idx] = max([rE,rp,rmu]);
names = {'E','p','mu'};
DIAG.limiter = names{idx};

%======================================================================
% Apply accept/reject
%======================================================================
if reject

    % Reject: keep accepted STATE unchanged
    NUM.dt_phy        = dt_next;
    NUM.dt_good_count = 0;
    DIAG.accept       = false;

else

    % Accept trial state
    STATE             = STATE_T;
    NUM.time          = NUM.time + DIAG.dt_try;
    NUM.dt_phy        = dt_next;
    NUM.dt_good_count = good_count_trial;
    DIAG.accept       = true;

end

end

function d = Max_Cell_Diff(A,B)
d = 0;
for i = 1:numel(A)
    d = max(d,max(abs(A{i}(:) - B{i}(:))));
end
end