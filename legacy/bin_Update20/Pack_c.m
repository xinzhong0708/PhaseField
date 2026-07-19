function [c] = Pack_c(c,ny)
for ip = 1:length(c)
    for ic = 1:length(c{ip})
        c{ip}{ic} = reshape(c{ip}{ic},ny,[]);
    end
end
end