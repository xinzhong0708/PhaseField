clear,figure(1);clf;addpath ./bin ./Utilities/ ./Solutions/ ./EOS
run_name = 'Mantle';
T        = [1500  1500  ] + 273.15;
P        = [1     1+1e-5] * 1e9;
solmod   = 'solution_models_H18';
dz       = 1/4;
eps_solv = 2;
Cname    = {'Si'  ,'Al' , 'Cr',    'Ti'    ,'Fe'   ,'Mn',    'Mg',   'Ca',   'Na',    'K',    'H'     ,'O' };
Oxname   = {'SiO2','Al2O3','Cr2O3','TiO2'  ,'FeO'  ,'MnO',   'MgO',  'CaO',  'Na2O',  'K2O',  'H2O'        };
noxy     = [2      3      3        2        1       1         1        1       1        1       1          ];
ncat     = [1      2      2        1        1       1         1        1       2        2       2          ];
molmOx   = [60.084  101.961 151.9904 79.8658  71.844  70.93744  40.304  56.077 61.97894 94.196  18.01528   ];
mass     = [28.085, 26.9815, 51.9961, 47.867, 55.845, 54.9380 , 24.305, 40.078,22.9898, 39.0983,1.0079, 15.999];
% Choose possible phases to consider in the equilibrium calculation (in the Gibbs minimization)
phs_name = {'Clinopyroxene'};

% Calculate free energy
td            = init_thermo(phs_name,Cname,solmod);
[g0,v0]       = tl_g0(T(1),P(1),td);
p             = props_generate(td);     % generate endmember proportions
[g,Npc,pc_id] = tl_gibbs_energy(T(1),P(1),phs_name,td,p,g0,v0);    

% Calculate free energy
[g0n,v0n]     = tl_g0(T(2),P(2),td);
[gn]          = tl_gibbs_energy(T(2),P(2),phs_name,td,p,g0n,v0n);    

%Find only Fe, Si, Mg
id            =(Npc(2,:)+Npc(3,:)+Npc(4,:)+Npc(6,:)+Npc(8,:)+Npc(9,:)+Npc(10,:)+Npc(11,:))>0;
Npc(:,id==1)  = [];
g(id)         = [];
gn(id)        = [];

%Normalize to 1 mole of cation
g             = g ./(Npc(1,:)+Npc(5,:)+Npc(7,:))';
gn            = gn./(Npc(1,:)+Npc(5,:)+Npc(7,:))';

%Calculate volume per mole Si
v             = 1.8886e-05;

%Energy per volume (J/m3)
g             =  g/v/1e5;
gn            = gn/v/1e5;

%Density
Npc           = Npc./repmat(Npc(1,:)+Npc(5,:)+Npc(7,:),12,1);
rho           = sum(Npc.*repmat(mass',1,length(Npc)))'./v/1e3;

%Concentration
X_Fe          = Npc(5,:)./sum(Npc(1:end-1,:),1);
X_Mg          = Npc(7,:)./sum(Npc(1:end-1,:),1);
X_Si          = Npc(1,:)./sum(Npc(1:end-1,:),1);

%Find minimal and sorting the raw data
[a,b]         = unique(round(X_Fe*1e12)/1e12);
for i = 1:length(a)
    id        = (X_Fe>a(i)-1e-3) .* (X_Fe<a(i)+1e-3);
    gu(i)     = min(g(id==1));
end
X_Fe          = X_Fe(b);
X_Mg          = X_Mg(b);
X_Si          = X_Si(b);
g             = gu;
[~,b]         = sort(X_Fe);
X_Fe          = X_Fe(b);
X_Mg          = X_Mg(b);
X_Si          = X_Si(b);
g             =    g(b);

%Construct free energy and add penalty
X_Fe_2d       = linspace(-0.2,0.7,51);
X_Mg_2d       = linspace(-0.2,0.7,51); [X_Fe_2d,X_Mg_2d]=ndgrid(X_Fe_2d,X_Mg_2d);
X_Fe_2d       = [X_Fe_2d(:) ; linspace(0.0001,0.4999,0)'];
X_Mg_2d       = [X_Mg_2d(:) ; linspace(0.4999,0.0001,0)'];
w             = 1e8;
g_2d          = interp1(X_Fe,g,X_Fe_2d)/2 + interp1(X_Mg,g,X_Mg_2d)/2;
g_2d          = g_2d + w*(0.5-X_Fe_2d-X_Mg_2d).^2;

%Make flag for good points
flag_good     = isfinite(g_2d);

%Fit good points
for i = 1:size(g_2d,1)
    %If good points: then only use the good point, do not use the protection points
    r               =  sqrt( (X_Fe_2d(i) - X_Fe_2d).^2 + (X_Mg_2d(i) - X_Mg_2d).^2 );
    if flag_good(i)==1
        %Find good points within the range
        idv         =  find(flag_good);
        %Sort the good points
        [~,rid]     =  sort(r(idv));
        rid         =  idv(rid);
        rid         =  rid(1:64);
        X_Fit       = [X_Fe_2d(rid,:) , X_Mg_2d(rid,:)];
        Y_Fit       =  g_2d(rid);
        %If bad point (protect ghost node)
        p           =  fit(X_Fit,Y_Fit,'poly22');
        pval{1}(i)  =  p.p00;
        pval{2}(i)  =  p.p10;
        pval{3}(i)  =  p.p01;
        pval{4}(i)  =  p.p20;
        pval{5}(i)  =  p.p11;
        pval{6}(i)  =  p.p02;
    end
end

%Fit bad points
w0    = 1e6;
for i = 1:size(g_2d,1)
    %If bad points: add protection term directly on p value
    if flag_good(i)==0
        r           =  sqrt((X_Fe_2d(i)-X_Fe_2d).^2 + (X_Mg_2d(i)-X_Mg_2d).^2);
        [~,rid]     =  sort(r);
        id          =  find(isfinite(g_2d(rid)));
        %The ID of the good point
        %g=p00+p10*X_Fe+p01*X_Mg+p20*X_Fe^2+p11*X_Fe*X_Mg+p02*X_Mg^2
        rid         =  rid(id(1));
        x           =  X_Fe_2d(i);
        y           =  X_Mg_2d(i);
        x0          =  X_Fe_2d(rid);
        y0          =  X_Mg_2d(rid);
        w           =  w0*r(rid);
        pval{1}(i)  =  pval{1}(rid) + w*(x0^2+y0^2);
        pval{2}(i)  =  pval{2}(rid) + w*(-2*x0);
        pval{3}(i)  =  pval{3}(rid) + w*(-2*y0);
        pval{4}(i)  =  pval{4}(rid) + w;
        pval{5}(i)  =  pval{5}(rid) + 0;        
        pval{6}(i)  =  pval{6}(rid) + w;
    end
end


%Testing and compare
X_Fe_test             = linspace(-0.1,0.6,51);
X_Mg_test             = linspace(-0.1,0.6,51);
[X_Fe_test,X_Mg_test] = ndgrid(X_Fe_test,X_Mg_test);
g_2d_high             = griddata(X_Fe_2d,X_Mg_2d,g_2d,X_Fe_test,X_Mg_test,'cubic');
subplot(121);
surf(X_Fe_test,X_Mg_test,g_2d_high);colorbar;shading interp
title('Energy');xlabel('Fe');ylabel('Mg')

%Evaluate
p00    = griddata(X_Fe_2d,X_Mg_2d,pval{1}(:,:),X_Fe_test,X_Mg_test);
p10    = griddata(X_Fe_2d,X_Mg_2d,pval{2}(:,:),X_Fe_test,X_Mg_test);
p01    = griddata(X_Fe_2d,X_Mg_2d,pval{3}(:,:),X_Fe_test,X_Mg_test);
p20    = griddata(X_Fe_2d,X_Mg_2d,pval{4}(:,:),X_Fe_test,X_Mg_test);
p11    = griddata(X_Fe_2d,X_Mg_2d,pval{5}(:,:),X_Fe_test,X_Mg_test);
p02    = griddata(X_Fe_2d,X_Mg_2d,pval{6}(:,:),X_Fe_test,X_Mg_test);
g_test = p00 + p10.*X_Fe_test + p01.*X_Mg_test + p20.*X_Fe_test.^2 + p11.*X_Fe_test.*X_Mg_test + p02.*X_Mg_test.^2;
hold on
plot3(X_Fe_test(:),X_Mg_test(:),g_test(:),'k.')

%Calculate derivatives
for i = 1:length(X_Fe_2d)
    dmu1dc1       = 2*pval{4}(i);
    dmu1dc2       =   pval{5}(i);
    dmu2dc1       =   pval{5}(i);
    dmu2dc2       = 2*pval{6}(i);
    D             = inv([dmu1dc1 dmu1dc2 ; dmu2dc1 dmu2dc2]);
    dcdmu{1,1}(i) = D(1,1);
    dcdmu{1,2}(i) = D(1,2);
    dcdmu{2,1}(i) = D(2,1);
    dcdmu{2,2}(i) = D(2,2);
end

subplot(122);
pcolor(X_Fe_test,X_Mg_test,(g_test-g_2d_high)./g_2d_high);colorbar;shading interp
title('Error');xlabel('Fe');ylabel('Mg')
% clim([-1 1]/10)

%Final result
%p order: p00 p10 p01 p20 p11 p02
%Equation: g=p00+p10*X_Fe+p01*X_Mg+p20*X_Fe^2+p11*X_Fe*X_Mg+p02*X_Mg^2
pval_Cpx = pval;
c_Fe_Cpx = X_Fe_2d;
c_Mg_Cpx = X_Mg_2d;
f_Cpx    = g_2d;
dcdmu_Cpx= dcdmu;
save('Data_Cpx','c_Fe_Cpx','c_Mg_Cpx','pval_Cpx','f_Cpx','dcdmu_Cpx')


% clf
% plot3(X_Fe_2d,X_Mg_2d,pval{1}(:,:),'.')
% % plot3(X_Fe_2d,X_Mg_2d,g_2d,'.')
% plot3(X_Fe_test,X_Mg_test,g_test,'.')







