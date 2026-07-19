function [chi] = Pack_chi(chi,ny)
for i = 1:size(chi,1)
    for j = 1:size(chi,2)
        chi{i,j} = reshape(chi{i,j},ny,[]);
    end
end
end