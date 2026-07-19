function [E] = Unpack_E(E)
for ie = 1:length(E)
    E{ie} = E{ie}(:)';
end
end