function [c] = Unpack_c(c)
for ip = 1:length(c)
    for ic = 1:length(c{ip})
        c{ip}{ic} = c{ip}{ic}(:)';
    end
end
end