function [MODEL,PARAM] = Update_Model_PT(MODEL,PARAM,PHYS,T,P)
%UPDATE_MODEL_PT Rebuild thermodynamic pars for present phases only.
%
% Full-name-only version.
%
% PARAM.update_PT_every controls how often thermodynamic pars are rebuilt.
% If update_PT_every = 1, update every call.
% If update_PT_every = 10, rebuild every 10 calls.
%
% Required MODEL fields:
%   MODEL.phs_name
%   MODEL.phase_index
%   MODEL.Cname
%   MODEL.solmod
%
% Required PHYS fields:
%   PHYS.E_sc
%
% Optional PHYS fields:
%   PHYS.vref

if isfield(PHYS,'vref')
    vref = PHYS.vref;
else
    vref = 2e-5;
end

% ------------------------------------------------------------
% Update control
% ------------------------------------------------------------
if isfield(PARAM,'update_PT_every')
    update_every = PARAM.update_PT_every;
else
    update_every = 1;
end

update_every = round(update_every);
if update_every < 1
    update_every = 1;
end

if ~isfield(PARAM,'PT_update_count') || isempty(PARAM.PT_update_count)
    PARAM.PT_update_count = 0;
end

PARAM.PT_update_count = PARAM.PT_update_count + 1;

% Current requested P-T, for diagnostics/path tracking
PARAM.T = T;
PARAM.P = P;

% Force first update
do_update = false;

if PARAM.PT_update_count == 1
    do_update = true;
end

% Periodic update
if update_every <= 1
    do_update = true;
elseif mod(PARAM.PT_update_count - 1,update_every) == 0
    do_update = true;
end

% Force update if MODEL.pars does not exist
if ~isfield(MODEL,'pars') || isempty(MODEL.pars)
    do_update = true;
end

% Force update if stored phase pars do not exist
if ~isfield(PARAM,'pars_phase_PT') || isempty(PARAM.pars_phase_PT)
    do_update = true;
end

% ------------------------------------------------------------
% Optional force update by temperature / pressure difference
% ------------------------------------------------------------
if isfield(PARAM,'PT_update_dT') && isfield(MODEL,'T')
    if abs(T - MODEL.T) >= PARAM.PT_update_dT
        do_update = true;
    end
end

if isfield(PARAM,'PT_update_dP') && isfield(MODEL,'P')
    if abs(P - MODEL.P) >= PARAM.PT_update_dP
        do_update = true;
    end
end

% ------------------------------------------------------------
% Skip expensive thermodynamic rebuild
% ------------------------------------------------------------
PARAM.PT_updated = do_update;

if ~do_update
    return
end

% ------------------------------------------------------------
% Rebuild only phases present in the current PFM model
% ------------------------------------------------------------
pars_phase = Build_Pars_Phases(MODEL.phs_name,MODEL.Cname,MODEL.solmod,T,P,PHYS.E_sc,vref);

% ------------------------------------------------------------
% Expand phase pars to grain-resolved MODEL.pars
% ------------------------------------------------------------
if max(MODEL.phase_index) > numel(pars_phase)
    error('Update_Model_PT: MODEL.phase_index exceeds number of present phases.')
end

MODEL.pars = pars_phase(MODEL.phase_index);

% ------------------------------------------------------------
% Store P-T state used by MODEL.pars
% ------------------------------------------------------------
MODEL.T = T;
MODEL.P = P;

PARAM.pars_phase_PT = pars_phase;
PARAM.pars_T        = T;
PARAM.pars_P        = P;

end