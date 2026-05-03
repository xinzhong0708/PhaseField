clear;figure(1);clf
%Load data
Phase = {'Cpx','Olv','Melt'};
for i = 1:length(Phase)
    load(['Data_',Phase{i}])
end
%Symbolics
syms C1 C2
a = sym('a',[3,6]);
c = sym('c',[3,2]);
p = sym('p',[3,1]);

%Free energy formulation
for i = 1:size(a,1)
    F(i) = a(i,1)+a(i,2).*c(i,1)+a(i,3).*c(i,2)+a(i,4)*c(i,1).^2+a(i,5).*c(i,1).*c(i,2)+a(i,6).*c(i,2).^2;
end

%Chemical potential
for i = 1:size(a,1)
    mu(i,1) = diff(F(i),c(i,1));
    mu(i,2) = diff(F(i),c(i,2));
end

%Equation 1: Fe conservation
Eq1 = C1 - p(1)*c(1,1) - p(2)*c(2,1) - p(3)*c(3,1);

%Equation 2: Mg conservation
Eq2 = C2 - p(1)*c(1,2) - p(2)*c(2,2) - p(3)*c(3,2);

%Equation 3:
Eq3 = mu(1,1) - mu(2,1);

%Equation 4:
Eq4 = mu(1,1) - mu(3,1);

%Equation 5:
Eq5 = mu(1,2) - mu(2,2);

%Equation 6:
Eq6 = mu(1,2) - mu(3,2);

%Collect terms
Eq  = [Eq1 Eq2 Eq3 Eq4 Eq5 Eq6];
c   = c(:);

%Formulate Jacobian
for i = 1:length(Eq)
    for j = 1:length(Eq)
        Term = coeffs(Eq(i),c(j));
        if length(Term)>1
            Jac_Func(i,j) = Term(2);
        end
    end
end

%Right hand side
Rhs_Func = simplify(Jac_Func*c(:) - Eq(:));

%Set values
C1    =  0.3;
C2    =  0.3;
c1_1  =  0.2;
c1_2  =  0.3;
c2_1  =  0.33;
c2_2  =  0.33;
c3_1  =  0.3;
c3_2  =  0.3;

%Proportion
phi   = [0.8 0.1 0.1];

p1    =  0.;%phi(1)^2/sum(phi.^2);
p2    =  1;%phi(2)^2/sum(phi.^2);
p3    =  1-p1-p2;

factor=  logspace(-4,0,100);
for ii = 1:length(factor)
    %Save old
    sol_o        = reshape(eval(c),3,[]);
    %Asign phases
    c1_int{1}    = c_Fe_Cpx;
    c2_int{1}    = c_Mg_Cpx;
    pv_int{1}    = pval_Cpx;
    c1_int{2}    = c_Fe_Olv;
    c2_int{2}    = c_Mg_Olv;
    pv_int{2}    = pval_Olv;
    c1_int{3}    = c_Fe_Melt;
    c2_int{3}    = c_Mg_Melt;
    pv_int{3}    = pval_Melt;
    a1_1         = griddata(c1_int{1},c2_int{1},pv_int{1}{1},c1_1,c1_2);
    a1_2         = griddata(c1_int{1},c2_int{1},pv_int{1}{2},c1_1,c1_2);
    a1_3         = griddata(c1_int{1},c2_int{1},pv_int{1}{3},c1_1,c1_2);
    a1_4         = griddata(c1_int{1},c2_int{1},pv_int{1}{4},c1_1,c1_2);
    a1_5         = griddata(c1_int{1},c2_int{1},pv_int{1}{5},c1_1,c1_2);
    a1_6         = griddata(c1_int{1},c2_int{1},pv_int{1}{6},c1_1,c1_2);

    a2_1         = griddata(c1_int{2},c2_int{2},pv_int{2}{1},c2_1,c2_2);
    a2_2         = griddata(c1_int{2},c2_int{2},pv_int{2}{2},c2_1,c2_2);
    a2_3         = griddata(c1_int{2},c2_int{2},pv_int{2}{3},c2_1,c2_2);
    a2_4         = griddata(c1_int{2},c2_int{2},pv_int{2}{4},c2_1,c2_2);
    a2_5         = griddata(c1_int{2},c2_int{2},pv_int{2}{5},c2_1,c2_2);
    a2_6         = griddata(c1_int{2},c2_int{2},pv_int{2}{6},c2_1,c2_2);

    a3_1         = griddata(c1_int{3},c2_int{3},pv_int{3}{1},c3_1,c3_2);
    a3_2         = griddata(c1_int{3},c2_int{3},pv_int{3}{2},c3_1,c3_2);
    a3_3         = griddata(c1_int{3},c2_int{3},pv_int{3}{3},c3_1,c3_2);
    a3_4         = griddata(c1_int{3},c2_int{3},pv_int{3}{4},c3_1,c3_2);
    a3_5         = griddata(c1_int{3},c2_int{3},pv_int{3}{5},c3_1,c3_2);
    a3_6         = griddata(c1_int{3},c2_int{3},pv_int{3}{6},c3_1,c3_2);


    %Solution update
    sol_n        = reshape(eval(Jac_Func)\eval(Rhs_Func),3,[]);
    sol          = sol_o*(1-factor(ii))+sol_n*factor(ii);
    
    %Update phase concentration
    c1_1         = sol(1,1);
    c1_2         = sol(1,2);
    c2_1         = sol(2,1);
    c2_2         = sol(2,2);
    c3_1         = sol(3,1);
    c3_2         = sol(3,2);

    plot(ii,sol(:),'+');hold on;drawnow

end

sol

%Grand potential
eval(F') - sum(eval(mu).*sol,2)

