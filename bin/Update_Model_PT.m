function [MODEL,PARAM] = Update_Model_PT(MODEL,PARAM,T,P)

%UPDATE_MODEL_PT Update MODEL.pars for current pressure-temperature.

pars_phase = Build_Pars_Phases_PT(PARAM.PT,T,P);

MODEL.pars = pars_phase(MODEL.phase_index);

PARAM.T = T;
PARAM.P = P;

end