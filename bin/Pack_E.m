function E = Pack_E(E,ny)
for ie = 1:length(E)
    E{ie} = reshape(E{ie},ny,[]);
end
end