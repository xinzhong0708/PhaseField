function pars_phase = Build_Pars_Phases(phs_name,Cname,solmod,T,P,E_sc,vref)

%BUILD_PARS_PHASES Build pars for present phases at given T and P.
%
% Full-name-only version.
%
% Recomputes:
%   td = init_thermo(...)
%   g0 = tl_g0(T,P,td)
%   n  = td.n_em(:,1:end-1)

Nphase     = numel(phs_name);
pars_phase = cell(1,Nphase);

for ip = 1:Nphase

    % ------------------------------------------------------------
    % Phase
    % ------------------------------------------------------------
    phase_name = phs_name(ip);

    td = init_thermo(phase_name,Cname,solmod);
    g0 = cell2mat(tl_g0(T,P,td));
    n  = td.n_em(:,1:end-1);

    % ------------------------------------------------------------
    % Thermodynamic values
    % ------------------------------------------------------------
    td.n_em(:,1:end-1) = n;

    pars            = td;
    pars.n          = n;
    pars.P          = P;
    pars.T          = T;
    pars.g0         = g0;
    pars.E_sc       = E_sc;
    pars.vref       = vref;
    pars.phase_name = phase_name;

    pars_phase{ip}  = pars;

end

end