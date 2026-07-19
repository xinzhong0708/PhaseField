function [cc] = Slice_c(c,ip,id)
%c is phase endmember fraction
% ip is the number of phase wanted, id is binary matrix
for i = 1:length(ip)
    for j = 1:length(c{ip(i)})
        cc{i}{j} = c{ip(i)}{j}(id==1);
    end
end
end
