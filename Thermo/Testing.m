%Testing
I   = 4;
J   = 1;
xx1 = {rand(4,4)/5 rand(4,4)/5 rand(4,4)/5 rand(4,4)/5};
dc  = 1e-6;
for i = 1:4
    for j = 1:4
        xx2        = xx1;
        xx2{j}     = xx2{j}+dc;
        MU1        = mu(xx1);
        MU2        = mu(xx2);
        dmudc(i,j) = (MU2{i}(I,J)-MU1{i}(I,J))/dc;
    end
end
A=S(xx1);
for i = 1:4
    for j = 1:4
        Sval(i,j) = A{i,j}(I,J);
    end
end
inv(dmudc)-Sval