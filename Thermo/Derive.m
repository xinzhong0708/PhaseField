clear
Np = 4;
Nc = 3;
T  = 0;
c  = sym('c',[Np,Nc]);
for i = 1:Np
    f(i) = T*sum(c(i,:).*log(c(i,:))) + T*(1-sum(c(i,:)))*log(1-sum(c(i,:))) + rand*c(i,1)*c(i,2) + rand*c(i,2)*c(i,3) + rand*c(i,1)*c(i,3);
end
for i = 1:Np
    for j = 1:Nc
        mu(i,j) = diff(f(i),c(i,j));
    end
end
syms C1 C2 C3
p     = rand(1,Np); p=p/sum(p);
Eq(1) = C1 - p*c(:,1);
Eq(2) = C2 - p*c(:,2);
Eq(3) = C3 - p*c(:,3);
Eq(4) = mu(1,1)-mu(2,1);
Eq(5) = mu(1,1)-mu(3,1);
Eq(6) = mu(1,1)-mu(4,1);
Eq(7) = mu(1,2)-mu(2,2);
Eq(8) = mu(1,2)-mu(3,2);
Eq(9) = mu(1,2)-mu(4,2);
Eq(10) = mu(1,3)-mu(2,3);
Eq(11) = mu(1,3)-mu(3,3);
Eq(12) = mu(1,3)-mu(4,3);

sol = solve(Eq,c(:));

simplify(sol.c1_1)


