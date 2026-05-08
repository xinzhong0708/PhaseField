function [STATE,NUM,DIAG] = Update_TimeStep_Soft(STATE,STATE_T,PARAM,MODEL,NUM)
%UPDATE_TIMESTEP Smooth accept/reject timestep controller.

dt_try = NUM.dt_phy;

%Final LE correction before timestep decision
STATE_C = LE_Run(STATE_T,PARAM,MODEL);

%Diagnostics on primary evolved fields
dE = Max_Cell_Diff(STATE_C.E,STATE.E);
dp = max(abs(STATE_C.p(:) - STATE.p(:)));

%Robust mu diagnostic
Dmu      = abs(cell2mat(STATE_C.mu_e) - cell2mat(STATE.mu_e));
dmu_rob  = Robust_Max(Dmu,0.99);
dmu_hard = max(Dmu,[],'all');

%Relative errors
err_E  = dE / max(NUM.dE_target,eps);
err_p  = dp / max(NUM.dp_target,eps);
err_mu = dmu_rob / max(NUM.dmu_target,eps);

err_primary = max(err_E,err_p);
err         = max([err_primary,err_mu]);

%Reject condition
hard_bad_mu = dmu_hard > 50*NUM.dmu_target && dmu_rob > 2*NUM.dmu_target;
bad_number  = any(~isfinite(Dmu(:))) || ~isfinite(dE) || ~isfinite(dp);

reject_step = err_primary > 1.0 || hard_bad_mu || bad_number;

%Smooth timestep update
if reject_step

    dt_next = max(dt_try*NUM.dt_shrink_fac,NUM.dt_min);
    NUM.dt_good_count = 0;
    % disp('Reject')
else
    % disp('Accept')
    if err < NUM.err_grow
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

%Diagnostics
DIAG.dE          = dE;
DIAG.dp          = dp;
DIAG.dmu_rob     = dmu_rob;
DIAG.dmu_hard    = dmu_hard;
DIAG.err_E       = err_E;
DIAG.err_p       = err_p;
DIAG.err_mu      = err_mu;
DIAG.err         = err;
DIAG.dt_try      = dt_try;
DIAG.dt_next     = dt_next;
DIAG.accept      = ~reject_step;
DIAG.hard_bad_mu = hard_bad_mu;

%Accept or reject
if reject_step

    %Reject: keep accepted STATE unchanged
    NUM.dt_phy = dt_next;

else

    %Accept final LE-corrected state
    STATE      = STATE_C;
    NUM.time   = NUM.time + dt_try;
    NUM.dt_phy = dt_next;

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


function val = Robust_Max(A,q)

a = abs(A(:));
a = a(isfinite(a));

if isempty(a)
    val = inf;
    return
end

a = sort(a);
id = max(1,min(numel(a),ceil(q*numel(a))));

val = a(id);

end