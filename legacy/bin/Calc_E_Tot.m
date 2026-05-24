function [E] = Calc_E_Tot(e,p)
%Size
Np   = size(p,3);
ny   = size(p,1); 
nx   = size(p,2); 
Ne   = length(e{1});
E    = repmat({zeros(ny,nx)}, 1, Ne); 
for ip = 1:Np
    for ie = 1:Ne
        E{ie} = E{ie} + (p(:,:,ip).*e{ip}{ie});
    end
end
end