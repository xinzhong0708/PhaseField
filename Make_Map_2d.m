clear;figure(1);clf;addpath([cd,'\bin'])

%Scaling
E_sc            =  1e8;
L_sc            =  1;
eta             =  1200e10/E_sc;

%Thermodynamic data
pars            =  Load_Data({'Olv','Cpx','Grt','Qtz'});
% pars            =  Load_Data({'Olv','Cpx'});

%Make F
F.p_fun         =  @(a  ,phi)  phi(:,:,a).^2./sum(phi.^2,3);
F.dpdphi        =  @(a,b,phi) (a==b)*2*phi(:,:,b)./sum(phi.^2,3) - 2*phi(:,:,a).*phi(:,:,b).^2./sum(phi.^2,3).^2;  %dp(b)dphi(a)

%Make map
Lx              =  1e-5;
Ly              =  1e-5;
nx              =  120;
ny              =  120;
x               =  linspace(0,Lx,nx);
y               =  linspace(0,Ly,ny);

%Composition
c{1}{1}         =  0.0063*ones(ny,nx);
c{1}{2}         =  0.2597*ones(ny,nx);
c{1}{3}         =  0.5133*ones(ny,nx);

c{2}{1}         =  0.0200*ones(ny,nx);
c{2}{2}         =  0.1200*ones(ny,nx);
c{2}{3}         =  0.0300*ones(ny,nx);
c{2}{4}         =  0.6000*ones(ny,nx);

c{3}{1}         =  0.4818*ones(ny,nx);
c{3}{2}         =  0.3943*ones(ny,nx);

c{4}{1}         =  1.0000*ones(ny,nx);


phi             =  zeros(ny,nx,length(pars));
md              =  ny/2;

phi(   1:md,    1:md  ,1)=  1;
phi(md+1:nx,    1:md  ,2)=  1;
phi(1:md-10, md+1:nx  ,3)=  1;
phi(md-9:nx, md+1:nx  ,4)=  1;

% phi(   1:md,    1:nx  ,1)=  1;
% phi(md+1:ny,    1:nx  ,2)=  1;
% phi(1:md-10, md+1:nx  ,3)=  1;
% phi(md-9:nx, md+1:nx  ,4)=  1;

%Calculate p
p               =  Calc_p(F,phi);

%Calculate element concentration e for each phase
e               =  Calc_e(pars,c);
E               =  Calc_E_Tot(e,p);

%Scale model
dx              =  x(2)-x(1);
dy              =  y(2)-y(1);
dx              =  dx/L_sc;
dy              =  dy/L_sc;
x               =   x/L_sc;
y               =   y/L_sc;

%Phase number and element number
Np              =  size(phi,3);
Ne              =  length(e{1});

%Local equilibrium
mu_e            =  repmat({zeros(ny,nx)}, 1,Ne);
chi             =  repmat({zeros(ny,nx)},Ne,Ne);
eta             =  eta*ones(ny,nx);

% Pair-pair initialization in pure phase regions for 2D
tol_pure        =  1e-12;
p_pair          =  5e-2;

for ip = 1:Np

    % Mask of pure ip region
    mask = p(:,:,ip) > 1 - tol_pure;

    if ~any(mask(:))
        continue
    end

    id    = find(mask);
    Nmask = numel(id);

    pid_other     = 1:Np;
    pid_other(ip) = [];

    for io = pid_other

        disp(['Checking: ',num2str([ip io])])

        % Build artificial pair-phase field only in the masked region
        pp          = p;
        pp_ip       = pp(:,:,ip);
        pp_io       = pp(:,:,io);
        pp_ip(mask) = 1 - p_pair;
        pp_io(mask) = p_pair;
        pp(:,:,ip)  = pp_ip;
        pp(:,:,io)  = pp_io;

        % Build local bulk composition from current c and artificial pp
        e_test      = Calc_e(pars,c);
        E_test      = Calc_E_Tot(e_test,pp);

        % Slice p into 1 x Nmask x Np form
        p_slice     = zeros(1,Nmask,Np);
        for jp = 1:Np
            tmp             = pp(:,:,jp);
            p_slice(1,:,jp) = tmp(mask).';
        end

        % Slice c and E using the same 2D mask
        c_slice  = Slice_c(c,1:Np,mask);
        E_slice  = Slice_E(E_test,mask);
        chi_temp = repmat({zeros(1,Nmask)},Ne,Ne);
        mue_temp = repmat({zeros(1,Nmask)},1,Ne);

        % Run local equilibrium only for the masked points
        ctemp    =  LE_Run(pars,p_slice,c_slice,E_slice,mue_temp,chi_temp,eta,[0.5 5000],[0.5 5000]);

        % Put back only the updated pair compositions
        for ic = 1:length(ctemp{io})
            tmp = c{io}{ic};
            tmp(mask) = ctemp{io}{ic}(:);
            c{io}{ic} = tmp;
        end
        for ic = 1:length(ctemp{ip})
            tmp = c{ip}{ic};
            tmp(mask) = ctemp{ip}{ic}(:);
            c{ip}{ic} = tmp;
        end
    end
end

% Recalculate p and E after pair-pair initialization
p = Calc_p(F,phi);
e = Calc_e(pars,c);
E = Calc_E_Tot(e,p);

%Save
save('Map2d','pars','E_sc','Lx','Ly','c','E','e','p','phi','eta','mu_e','chi','x','dx','nx','y','dy','ny','L_sc','F','Np','Ne')
