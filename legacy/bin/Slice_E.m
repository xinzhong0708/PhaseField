function [EE] = Slice_E(E,id)
for ie = 1:length(E)
    EE{ie} = E{ie}(id);
end
end