function [STATE,NUM,DIAG] = Update_TimeStep_Soft(STATE,STATE_T,PARAM,MODEL,NUM)
%UPDATE_TIMESTEP Smooth accept/reject timestep controller.
%
% Main rule:
%   dmu controls hard stability.
%   dp and dE control timestep growth, but do not automatically reject
%   unless they exceed optional hard caps.

dt_try = NUM.dt_phy;

% ------------------------------------------------------------
% Defaults
% ------------------------------------------------------------
if ~isfield(NUM,'dE_target') || isempty(NUM.dE_target)
    NUM.dE_target = 1e-4;
end

if ~isfield(NUM,'dp_target') || isempty(NUM.dp_target)
    NUM.dp_target = 1e-3;
end

if ~isfield(NUM,'dmu_target') || isempty(NUM.dmu_target)
    NUM.dmu_target = 0.05;
end

if ~isfield(NUM,'dmu_hard_cap') || isempty(NUM.dmu_hard_cap)
    NUM.dmu_hard_cap = 3*NUM.dmu_target;
end

if ~isfield(NUM,'dp_hard_cap') || isempty(NUM.dp_hard_cap)
    NUM.dp_hard_cap = inf;
end

if ~isfield(NUM,'dE_hard_cap') || isempty(NUM.dE_hard_cap)
    NUM.dE_hard_cap = inf;
end

if ~isfield(NUM,'dt_shrink_fac') || isempty(NUM.dt_shrink_fac)
    NUM.dt_shrink_fac = 0.5;
end

if ~isfield(NUM,'dt_shrink_min_fac') || isempty(NUM.dt_shrink_min_fac)
    NUM.dt_shrink_min_fac = 0.1;
end

if ~isfield(NUM,'dt_grow_fac') || isempty(NUM.dt_grow_fac)
    NUM.dt_grow_fac = 1.05;
end

if ~isfield(NUM,'dt_grow_after') || isempty(NUM.dt_grow_after)
    NUM.dt_grow_after = 10;
end

if ~isfield(NUM,'err_grow') || isempty(NUM.err_grow)
    NUM.err_grow = 0.3;
end

if ~isfield(NUM,'dt_good_count') || isempty(NUM.dt_good_count)
    NUM.dt_good_count = 0;
end

% ------------------------------------------------------------
% Current physical time
% ------------------------------------------------------------
if ~isfield(NUM,'t_phy') || isempty(NUM.t_phy)
    if isfield(NUM,'time') && ~isempty(NUM.time)
        NUM.t_phy = NUM.time;
    else
        NUM.t_phy = 0;
    end
end

if ~isfield(NUM,'time') || isempty(NUM.time)
    NUM.time = NUM.t_phy;
end

t_old   = NUM.t_phy;
t_trial = t_old + dt_try;

% ------------------------------------------------------------
% Trial state
% ------------------------------------------------------------
STATE_C = STATE_T;

% ------------------------------------------------------------
% Diagnostics on primary evolved fields
% ------------------------------------------------------------
dE = Max_Cell_Diff(STATE_C.E,STATE.E);
dp = max(abs(STATE_C.p(:) - STATE.p(:)));

% ------------------------------------------------------------
% Mu diagnostic
% ------------------------------------------------------------
[dmu_rob,dmu_hard,mu_bad_number] = Mu_Cell_Diff(STATE_C.mu_e,STATE.mu_e,0.99);

% ------------------------------------------------------------
% Relative errors
% ------------------------------------------------------------
err_E        = dE / max(NUM.dE_target,eps);
err_p        = dp / max(NUM.dp_target,eps);
err_mu       = dmu_rob / max(NUM.dmu_target,eps);
err_mu_hard  = dmu_hard / max(NUM.dmu_hard_cap,eps);

% This error is used only for growth decision.
err_growth = max([err_E,err_p,err_mu]);

% ------------------------------------------------------------
% Solver failure flags
% ------------------------------------------------------------
thermo_bad = false;
acch_bad   = false;
chle_bad   = false;

if isfield(STATE_C,'LE_failed') && STATE_C.LE_failed
    thermo_bad = true;
end

% if isfield(STATE_C,'LE_diag') && isfield(STATE_C.LE_diag,'failed')
%     if STATE_C.LE_diag.failed
%         thermo_bad = true;
%     end
% end

if isfield(STATE_C,'ACCH_diag')
    D = STATE_C.ACCH_diag;
    if isfield(D,'solve_flag') && D.solve_flag ~= 0
        acch_bad = true;
    end
end

if isfield(STATE_C,'CHLE_diag')
    D = STATE_C.CHLE_diag;

    if isfield(D,'solve_flag') && D.solve_flag ~= 0
        chle_bad = true;
    end

    % Do NOT hard reject on edge_dmu_ratio or outside_res_ratio here.
    % These are useful diagnostics, but too sensitive for timestep rejection.
end

% ------------------------------------------------------------
% Reject conditions
% ------------------------------------------------------------
bad_number   = mu_bad_number || ~isfinite(dE) || ~isfinite(dp);
bad_mu_rob   = dmu_rob  > NUM.dmu_target;
bad_mu_hard  = dmu_hard > NUM.dmu_hard_cap;

bad_p_hard   = dp > NUM.dp_hard_cap;
bad_E_hard   = dE > NUM.dE_hard_cap;

reject_step = bad_number || bad_mu_rob || bad_mu_hard || ...
              bad_p_hard || bad_E_hard || ...
              thermo_bad || acch_bad || chle_bad;

% ------------------------------------------------------------
% Timestep update
% ------------------------------------------------------------
if reject_step

    fac = NUM.dt_shrink_fac;

    % If mu is the reason, shrink proportional to overshoot.
    if bad_mu_rob || bad_mu_hard
        fac_mu = 0.8 / max([err_mu,err_mu_hard,1.0]);
        fac    = min(fac,fac_mu);
    end

    % Solver failure or bad numbers need stronger shrink.
    if bad_number || thermo_bad || acch_bad || chle_bad
        fac = min(fac,0.2);
    end

    fac = max(fac,NUM.dt_shrink_min_fac);

    dt_next = max(dt_try*fac,NUM.dt_min);
    NUM.dt_good_count = 0;

else

    % Only grow if all monitored errors are comfortably small.
    if err_growth < NUM.err_grow
        NUM.dt_good_count = NUM.dt_good_count + 1;
    else
        NUM.dt_good_count = 0;
    end

    if NUM.dt_good_count >= NUM.dt_grow_after
        dt_next = min(dt_try*NUM.dt_grow_fac,NUM.dt_max);
        NUM.dt_good_count = 0;
    else
        dt_next = dt_try;
    end

end

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
DIAG.dE           = dE;
DIAG.dp           = dp;
DIAG.dmu_rob      = dmu_rob;
DIAG.dmu_hard     = dmu_hard;

DIAG.err_E        = err_E;
DIAG.err_p        = err_p;
DIAG.err_mu       = err_mu;
DIAG.err_mu_hard  = err_mu_hard;
DIAG.err_growth   = err_growth;

DIAG.bad_number   = bad_number;
DIAG.bad_mu_rob   = bad_mu_rob;
DIAG.bad_mu_hard  = bad_mu_hard;
DIAG.bad_p_hard   = bad_p_hard;
DIAG.bad_E_hard   = bad_E_hard;

DIAG.thermo_bad   = thermo_bad;
DIAG.acch_bad     = acch_bad;
DIAG.chle_bad     = chle_bad;

DIAG.dt_try       = dt_try;
DIAG.dt_next      = dt_next;
DIAG.accept       = ~reject_step;

DIAG.t_old        = t_old;
DIAG.t_trial      = t_trial;

% Main reason string
if bad_number
    DIAG.reject_reason = 'bad_number';
elseif bad_mu_hard
    DIAG.reject_reason = 'bad_mu_hard';
elseif bad_mu_rob
    DIAG.reject_reason = 'bad_mu_rob';
elseif bad_p_hard
    DIAG.reject_reason = 'bad_p_hard';
elseif bad_E_hard
    DIAG.reject_reason = 'bad_E_hard';
elseif thermo_bad
    DIAG.reject_reason = 'thermo_bad';
elseif acch_bad
    DIAG.reject_reason = 'acch_bad';
elseif chle_bad
    DIAG.reject_reason = 'chle_bad';
else
    DIAG.reject_reason = 'accepted';
end

DIAG.t_old        = t_old;
DIAG.t_trial      = t_trial;

% ------------------------------------------------------------
% Accept or reject
% ------------------------------------------------------------
if reject_step

    NUM.t_phy  = t_old;
    NUM.time   = t_old;
    NUM.dt_phy = dt_next;

    DIAG.t_new = t_old;

else

    STATE      = STATE_C;
    NUM.t_phy  = t_old + dt_try;
    NUM.time   = NUM.t_phy;
    NUM.dt_phy = dt_next;

    DIAG.t_new = NUM.t_phy;

end

end


function d = Max_Cell_Diff(A,B)

d = 0;

for i = 1:numel(A)

    diff_i = abs(A{i}(:) - B{i}(:));

    if any(~isfinite(diff_i))
        d = inf;
        return
    end

    d = max(d,max(diff_i));

end

end


function [drob,dhard,bad_number] = Mu_Cell_Diff(A,B,q)

bad_number = false;
dhard      = 0;
vals       = [];

for i = 1:numel(A)

    diff_i = abs(A{i}(:) - B{i}(:));

    if any(~isfinite(diff_i))
        bad_number = true;
    end

    diff_i = diff_i(isfinite(diff_i));

    if isempty(diff_i)
        continue
    end

    dhard = max(dhard,max(diff_i));

    vals = [vals; diff_i];

end

if isempty(vals)
    drob = inf;
    dhard = inf;
    bad_number = true;
    return
end

vals = sort(vals);
id   = max(1,min(numel(vals),ceil(q*numel(vals))));

drob = vals(id);

end