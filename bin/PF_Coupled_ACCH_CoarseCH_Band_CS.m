function [STATE_NEW,DIAG] = PF_Coupled_ACCH_CoarseCH_Band_CS(STATE_REF,PARAM,MODEL,GRID,PHYS,NUM,STATE_COEF)
%PF_COUPLED_ACCH_COARSECH_BAND_CS
%
% One-step predictor:
%   1. coarse-grid global CH background
%   2. fine-grid local ACCH correction in an interface band
%
% This is a convenience wrapper around:
%   PF_CH_Coarse_Background_CS
%   PF_Coupled_ACCH_LETangent_BandCorrection_CS

if nargin < 7 || isempty(STATE_COEF)
    STATE_COEF = STATE_REF;
end

[STATE_BG,DIAG_BG]   = PF_CH_Coarse_Background_CS(STATE_REF,PARAM,MODEL,GRID,PHYS,NUM,STATE_COEF);

[STATE_NEW,DIAG_LOC] = PF_Coupled_ACCH_LETangent_BandCorrection_CS(STATE_REF,STATE_BG,PARAM,MODEL,GRID,PHYS,NUM,STATE_BG);

DIAG.background = DIAG_BG;
DIAG.local      = DIAG_LOC;

DIAG.max_dphi   = DIAG_LOC.max_dphi;
DIAG.max_dmu    = DIAG_LOC.max_dmu;
DIAG.max_dE     = DIAG_LOC.max_dE;
DIAG.relres     = DIAG_LOC.relres;
DIAG.matrix_size= DIAG_LOC.matrix_size;
DIAG.nnz        = DIAG_LOC.nnz;

end
