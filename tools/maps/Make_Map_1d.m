clear;figure(1);clf;addpath([cd,'\bin'])

%Scaling
E_sc            =  1e8;
L_sc            =  1;
eta             =  1000e10/E_sc;

%Thermodynamic data
pars            =  Load_Data({'Olv','Cpx','Grt','Qtz'});
% pars            =  Load_Data({'Olv','Cpx'});
% pars            =  Load_Data({'Cpx'});

%Make F
F.p_fun         =  @(a  ,phi)  phi(:,:,a).^2./sum(phi.^2,3);
F.dpdphi        =  @(a,b,phi) (a==b)*2*phi(:,:,b)./sum(phi.^2,3) - 2*phi(:,:,a).*phi(:,:,b).^2./sum(phi.^2,3).^2;  %dp(b)dphi(a)

%Make map
Lx              =  1e-5;
nx              =  400;
ny              =  1;
x               =  linspace(0,Lx,nx);

% %Composition
c{1}{1}         =  0.0063*ones(ny,nx);
c{1}{2}         =  0.2597*ones(ny,nx);
c{1}{3}         =  0.5133*ones(ny,nx);

c{2}{1}         =  0.0200*ones(ny,nx);
c{2}{2}         =  0.1200*ones(ny,nx);
c{2}{3}         =  0.0300*ones(ny,nx);
c{2}{4}         =  0.6000*ones(ny,nx);


c{3}{1}         =  0.3800*ones(ny,nx);
c{3}{2}         =  0.4800*ones(ny,nx);

c{4}{1}         =  1.0000*ones(ny,nx);



phi                  =  zeros(1,nx,length(c));
nn                   =  nx/4;
phi(1,1:nn,1)        =  1;
phi(1,nn+1:2*nn,2)   =  1;
phi(1,2*nn+1:3*nn,4) =  1;
phi(1,3*nn+1:4*nn,3) =  1;

% phi(1,1:100,1)  =  1;
% phi(1,101:200,2) =  1;

% c{1}{1}(300:nx) = c{1}{1}(300:nx) + 0.06;
% c{1}{2}(300:nx) = c{1}{2}(300:nx) + 0.06;
% c{1}{3}(300:nx) = c{1}{3}(300:nx) - 0.06;
% c{1}{4}(300:nx) = c{1}{4}(300:nx) - 0.06;

%Calculate p
p               =  Calc_p(F,phi);

%Calculate element concentration e for each phase
e               =  Calc_e(pars,c);
E               =  Calc_E_Tot(e,p);

%Scale model
dx              =  x(2)-x(1);
dx              =  dx/L_sc;
x               =   x/L_sc;

%Phase number and element number
Np              =  size(phi,3);
Ne              =  length(e{1});

%Local equilibrium
mu_e            =  repmat({zeros(1,nx)}, 1,Ne);
chi             =  repmat({zeros(1,nx)},Ne,Ne);

%Calculate element concentration e for each phase
e               =  Calc_e(pars,c);
E               =  Calc_E_Tot(e,p);

for ip = 1:Np
    id              =  find(p(:,:,ip)==1);
    pid_other       =  1:Np;
    pid_other(ip)   =  [];
    for io = pid_other
        pp          =  p;
        pp(:,id,ip) =  1-3e-2;
        pp(:,id,io) =  3e-2;
        e_test      =  Calc_e(pars,c);
        E_test      =  Calc_E_Tot(e_test,pp);
        %Slice everything
        p_slice     =  pp(:,id,:);
        c_slice     =  Slice_c(c,1:Np,p(:,:,ip)==1);
        E_slice     =  Slice_E(E_test,id);
        chi_temp    =  repmat({zeros(1,length(id))},Ne,Ne);
        mue_temp    =  repmat({zeros(1,length(id))}, 1,Ne);
        ctemp       =  LE_Run(pars,p_slice,c_slice,E_slice,mue_temp,chi_temp,eta*ones(1,length(id)),[0.5 5000],[0.5 5000]);
        %Put back c
        for ic = 1:length(ctemp{io})
            c{io}{ic}(:,id) = ctemp{io}{ic};
        end
        for ic = 1:length(ctemp{ip})
            c{ip}{ic}(:,id) = ctemp{ip}{ic};
        end
    end
end
p               =  Calc_p(F,phi);
e               =  Calc_e(pars,c);
E               =  Calc_E_Tot(e,p);

%Save
save('Map1d','E_sc','Lx','pars','c','E','e','p','phi','eta','mu_e','chi','x','dx','nx','ny','L_sc','F','Np','Ne')


%Optional plot
for ip = 1:Np
    subplot(1,Np,ip); hold on
    for ic = 1:length(c{ip})
        plot(x, c{ip}{ic})
    end
    title(['Phase ',num2str(ip)])
    xlabel('x')
    ylabel('c')
    hold off
end
drawnow

