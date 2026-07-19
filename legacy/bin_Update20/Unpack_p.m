function [pp] = Unpack_p(p)
pp = zeros(1,size(p,1)*size(p,2),size(p,3));
for ip = 1:size(p,3)
    pp(:,:,ip) = reshape(p(:,:,ip),1,[]);
end
end