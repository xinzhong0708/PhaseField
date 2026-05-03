function [e] = Calc_e(pars,c)
%Calculate e from c
Ny  =  size(c{1}{1},1);
Np  =  length(c);
c   =  Unpack_c(c);
for ip = 1:Np
    c_temp  =  cell2mat(c{ip}');
    %For pure phase
    if size(pars{ip}.n,1)==1
        c_temp  =  c_temp;
    else
        c_temp  = [c_temp ; 1-sum(c_temp,1)];
    end
    n           = [pars{ip}.n ; pars{ip}.nN];
    e_temp      = (c_temp'*n)';
    e_temp      =  e_temp./sum(e_temp,1);
    e{ip}       =  num2cell(e_temp(1:end-1,:),2)';
end
e = Pack_c(e,Ny);
end



