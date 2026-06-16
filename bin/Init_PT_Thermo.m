function PT = Init_PT_Thermo(phs_name,PT)

%INIT_PT_THERMO Initialize thermodynamic phase objects for later P-T update.
%
% This function does the expensive / structural part once:
%   init_thermo
%   endmember matrix perturbation
%
% Later Build_Pars_Phases_PT only updates:
%   T, P, g0

Nphase = numel(phs_name);

PT.td_base    = cell(1,Nphase);
PT.phase_full = cell(1,Nphase);

for ip = 1:Nphase

    name_short = phs_name{ip};

    id = find(strcmp(PT.phase_short,name_short),1,'first');

    if isempty(id)
        id = find(strcmp(PT.phase_all,name_short),1,'first');
    end

    if isempty(id)
        error('Phase name %s was not found in PT.phase_short or PT.phase_all.',name_short)
    end

    name_full = PT.phase_all{id};

    td                 = init_thermo({name_full},PT.Cname,PT.solmod);
    n                  = td.n_em(:,1:end-1);
    td.n_em(:,1:end-1) = n;
    td.n               = n;
    td.E_sc            = PT.E_sc;
    td.vref            = PT.vref;
    td.phase_name      = {name_full};

    PT.td_base{ip}     = td;
    PT.phase_full{ip}  = name_full;

end

end