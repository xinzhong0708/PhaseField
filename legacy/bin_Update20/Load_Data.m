function [pars] = Load_Data(phs_name)
%Load thermodynamic data
for ip = 1:length(phs_name)
    load(['Data_',phs_name{ip}])
    pars{ip}  =  par;
end
end
