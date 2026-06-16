function [MODEL,PARAM] = Update_Model_PT(MODEL,PARAM,PHYS,T,P)

%UPDATE_MODEL_PT Rebuild thermodynamic pars for present phases only.
%
% Full-name-only version.
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
% Rebuild only phases present in the current PFM model
% ------------------------------------------------------------
pars_phase = Build_Pars_Phases( ...
    MODEL.phs_name,MODEL.Cname,MODEL.solmod,T,P,PHYS.E_sc,vref);

% ------------------------------------------------------------
% Expand phase pars to grain-resolved MODEL.pars
% ------------------------------------------------------------
if max(MODEL.phase_index) > numel(pars_phase)
    error('Update_Model_PT: MODEL.phase_index exceeds number of present phases.')
end

MODEL.pars = pars_phase(MODEL.phase_index);

% ------------------------------------------------------------
% Store current P-T state
% ------------------------------------------------------------
MODEL.T = T;
MODEL.P = P;

PARAM.T = T;
PARAM.P = P;

PARAM.pars_phase_PT = pars_phase;

end