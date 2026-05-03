function [omg,G] = Calc_Omega(c,e,mu_e,F)
%Calculate energy and grand potential
G   = zeros(size(c{1}{1},1),size(c{1}{1},2),length(c));
omg = zeros(size(c{1}{1},1),size(c{1}{1},2),length(c));
for ip = 1:length(c)
    G(  :,:,ip)  =  F.g{ip}(c{ip});
    omg(:,:,ip)  =  real(G(:,:,ip));
    for ie = 1:length(e{1})
        omg(:,:,ip) = omg(:,:,ip) - (mu_e{ie}.*e{ip}{ie});
    end
end
end