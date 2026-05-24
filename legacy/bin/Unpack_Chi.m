function [chi] = Unpack_Chi(chi)
for i = 1:size(chi,1)
    for j = 1:size(chi,2)
        chi{i,j} = chi{i,j}(:)';
    end
end
end