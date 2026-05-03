%Clear and restart
clear;figure(1);clf;addpath([cd,'\bin']);addpath([cd,'\Thermo']);addpath([cd,'\Thermo\Solutions'])

%Load map
load Map1d.mat

%Scales (E_sc from map)
t_sc            =  1;                         % Time scale vary depending on D
L_sc            =  1;                         % Length scale, fixed to 1m (due to g=J/m3)

%Physical parameters
l               =  Lx/100/L_sc;               % interface thickness (m)
sigma           =   1.0/E_sc*L_sc^2;          % surface energy (J/m^2)
kappa           =  1e-7/E_sc*L_sc;            % 4th order term, can be set to 0 if no solvus
M0              =  1e-16*E_sc/L_sc^5/t_sc;    % Diffusion mobility

%Dependent variables
m               =  6*sigma/l;
kap             =  3/4*sigma*l;
dceq            =  0.2;
L               =  4*m/3/kap/(4*dceq^2/M0)/1e2;
disp('Time scale of interface')
disp(l^2/L/kap)
disp(1/L/m)
pause(1)

%Time step
dt_phy          =   1e-9/t_sc;
dt_max          =    1e2/t_sc;
dt_min          =  1e-16/t_sc;
t_tot           =    1e5/t_sc;
dE_target       =  2e-2;
dp_target       =  2e-2;
dmu_target      =  1e5;
time            =  0;

%Load map
load Map1d.mat

%Load interface
F.dgdphi        =  @(phi) 2*m*phi.*(phi - 1).^2 + m*phi.^2.*(2.*phi - 2);

%Phase transition L*kappa matrix
LL              =   L*ones(1,nx);
Lm              =  LL.*m;
LK              =   L*kap*ones(1,nx);

%Initialize Source germ
omg             =  zeros(1,nx,Np);
omg_t           =  omg;
s               =  Calc_S_AllenCahn(phi,p,LL,F,omg);
src             =  Calc_S_Diffusion(p,p,e,dt_phy);
[c_t,mu_e]      =  LE_Run(pars, p, c, E, mu_e, chi, eta*ones(1,nx), [1 1000], [1 1000]);
%Mobility
M               =  M0;

%Time loop
dt_phy_0        =  dt_phy;
LE_state        =  [];
%Display E
disp('Composition')
disp([mean(E{1}) mean(E{2}) mean(E{3}) mean(E{4}) ])

% Simple stable timestep controller
dt_good_count   = 0;

dt_grow_after   = 8;      % grow only after 8 very good accepted steps
dt_grow_fac     = 1.05;   % slow growth
dt_shrink_fac   = 0.5;    % fast shrink after rejection
err_grow        = 0.25;   % only count as "good" if error < 25% of target

for it = 1:1e6

    tic


    % Accepted old state
    Eo        = E;
    co        = c;
    mu_eo     = mu_e;
    chio      = chi;
    po        = p;
    phio      = phi;
    omgo      = omg;
    dt_try    = dt_phy;

    % Save LE state before trial
    LE_state_old = LE_state;


    % TRIAL: Allen-Cahn
    % Build Allen-Cahn source from accepted old state
    s     = Calc_S_AllenCahn(phio, po, LL, F, omgo);
    phi_t = zeros(size(phio));
    for ip = 1:Np
        phi_t(:,:,ip) = Solver_1D_AllenCahn_Periodic(s{ip}, phio(:,:,ip), LK, dx, dt_try);
    end
    phi_t = Norm_Phi(phi_t);
    p_t = Calc_p(F, phi_t);


    % Variable penalty
    eta_vec = kappa_from_p_smooth_full(p_t, eta, eta*0.1);


    % First LE: use trial p_t but old bulk composition Eo
    [c_t, mu_e_t, chi_t, omg_t, LE_state_try] = LE_Run(pars, p_t, co, Eo, mu_eo, chio, eta_vec, [1 1000], [1 1000], LE_state_old);


    % CH trial: diffusion using LE-updated state
    e_t = Calc_e(pars, c_t);
    src = Calc_S_Diffusion(p_t, po, e_t, dt_try);
    [mu_t, E_t] = Solver_1D_Diffusion(Eo, mu_e_t, M, chi_t, kappa, dx, dt_try, src);


    % Final LE correction BEFORE timestep decision
    [c_corr, mu_corr, chi_corr, omg_corr, LE_state_corr] = LE_Run(pars, p_t, c_t, E_t, mu_t, chi_t, eta_vec,[1 1000], [1 1000], LE_state_try);


    % Diagnostics on the actual accepted-candidate state
    dE  = max(abs(cell2mat(E_t)      - cell2mat(Eo)),    [], 'all');
    dp  = max(abs(p_t(:) - po(:)));
    dmu = max(abs(cell2mat(mu_corr)  - cell2mat(mu_eo)), [], 'all');

    % Normalized timestep error
    err_E  = dE  / max(dE_target,  eps);
    err_p  = dp  / max(dp_target,  eps);
    err_mu = dmu / max(dmu_target, eps);

    err = max([err_E, err_p, err_mu]);


    % Simple stable timestep controller
    dt_next     = dt_try;
    reject_step = err > 1.0;

    if reject_step
        % Bad step: reject and shrink immediately
        dt_next = max(dt_try * dt_shrink_fac, dt_min);
        dt_good_count = 0;
    else
        % Accepted step
        if err < err_grow
            dt_good_count = dt_good_count + 1;
        else
            dt_good_count = 0;
        end
        % Grow only after several very good accepted steps
        if dt_good_count >= dt_grow_after
            dt_next = min(dt_try * dt_grow_fac, dt_max);
            dt_good_count = 0;
        else
            dt_next = dt_try;
        end

    end


    % ============================================================
    % Accept or reject
    % ============================================================

    if reject_step
        % Reject: keep accepted state unchanged
        dt_phy   = dt_next;
        LE_state = LE_state_old;
    else
        % Accept final LE-corrected state
        phi      = phi_t;
        p        = p_t;
        c        = c_corr;
        mu_e     = mu_corr;
        chi      = chi_corr;
        omg      = omg_corr;
        E        = E_t;
        e        = Calc_e(pars, c);
        LE_state = LE_state_corr;
        time     = time + dt_try;
        dt_phy   = dt_next;
    end



    toc


    % Plotting / history
    TIME(it)   = time;
    PHASE1(it) = mean(p(:,:,1),'all');
    PHASE2(it) = mean(p(:,:,2),'all');

    if mod(it,5)==0
        disp(dt_phy)
        disp([mean(p(:,:,1),'all'), mean(E_t{1}(:))])
        subplot(341); plot(x,E{1}(1,:),x,E{2}(1,:),x,E{3}(1,:),x,E{end}(1,:)); title('C')
        subplot(342); plot(x,mu_e{1},x,mu_e{2},x,mu_e{3},x,mu_e{end}); title('mu')
        subplot(343); plot(x,phi(1,:,1),'.-',x,phi(1,:,2),'.-',x,phi(1,:,3),'.-',x,phi(1,:,end),'.-'); title(num2str([mean(p(:,:,1),'all'),mean(p(:,:,2),'all')]))
        subplot(345); plot(x,omg(1,:,1)-omg(1,:,2),'.-',x,omg(1,:,1)-omg(1,:,end),'.-',x,omg(1,:,2)-omg(1,:,3),'.-',x,omg(1,:,3)-omg(1,:,4),'.-'); title('\omega')
        subplot(346); plot(x,c{1}{1}(1,:),x,c{1}{2}(1,:),x,c{1}{3}(1,:)); title('c ph1')
        subplot(347); plot(x,c{2}{1}(1,:),x,c{2}{2}(1,:),x,c{2}{3}(1,:),x,c{2}{4}(1,:)); title('c ph2')
        subplot(348); plot(it,dt_phy,'b.');hold on; title('dt')
        subplot(3,4,10); plot(x,eta_vec,'.-')
        subplot(3,4,11); plot(TIME,PHASE1,'.-')
        subplot(3,4,12); plot(TIME,PHASE2,'.-')
        drawnow
    end
end




function [mu,C] = Solver_1D_Diffusion(co, muo, M, chi, kappa, dx, dt, s, chiRelFloor, chiAbsFloor)
% Mu-only fully implicit solver with periodic BC and Chi closure:
%
%   C^n = C^o + Chi*(mu^n - mu^o)
%
% Diffusion + gradient term + source:
%   (C^n - C^o)/dt = M * lap(mu^n) - M*kappa*lap^2(C^n) - s
%
% We eliminate C^n using the closure and solve only for mu^n, then recover C^n.
%
% s{a}(i) is an externally-provided source term (units of dC_a/dt)

if nargin < 8 || isempty(s), s = []; end
if nargin < 9 || isempty(chiRelFloor), chiRelFloor = 0e-12; end
if nargin < 10 || isempty(chiAbsFloor), chiAbsFloor = 0e-14; end

n  = length(co{1});
Nc = numel(co);
Nd = Nc;

idf = @(ii) (ii>0)*(ii<=Nd*n)*ii + (ii<=0)*(ii+Nd*n) + (ii>Nd*n)*(ii-Nd*n);

L   = sparse(Nd*n, Nd*n);
R   = zeros(Nd*n, 1);

dx2 = dx*dx;
dx4 = dx^4;

% ---- normalize/interpret s ----
% Accept:
%   - [] (no source)
%   - cell(1,Nc) with vectors length n
%   - numeric Nc x n or n x Nc
if ~isempty(s) && ~iscell(s)
    if ismatrix(s)
        if size(s,1) == Nc && size(s,2) == n
            sMat = s;
        elseif size(s,1) == n && size(s,2) == Nc
            sMat = s.'; % -> Nc x n
        else
            error('s must be cell(1,Nc) or numeric (Nc x n) or (n x Nc).');
        end
        sCell = cell(1,Nc);
        for a = 1:Nc
            sCell{a} = sMat(a,:);
        end
        s = sCell;
    else
        error('s must be cell(1,Nc) or numeric matrix.');
    end
end
if isempty(s)
    s = cell(1,Nc);
    for a = 1:Nc, s{a} = zeros(1,n); end
else
    for a = 1:Nc
        s{a} = reshape(s{a}, 1, []);
        if numel(s{a}) ~= n
            error('s{%d} must have length n=%d.', a, n);
        end
    end
end

% ---- regularize Chi at every grid point ----
ChiReg = cell(Nc,Nc);
for a = 1:Nc
    for b = 1:Nc
        ChiReg{a,b} = zeros(1,n);
    end
end

for i = 1:n
    ChiMat = zeros(Nc,Nc);
    for a = 1:Nc
        for b = 1:Nc
            ChiMat(a,b) = chi{a,b}(i);
        end
    end
    ChiMat = 0.5*(ChiMat + ChiMat.');  % symmetrize

    [V,D]  = eig(ChiMat);
    lam    = diag(D);
    smax   = max(abs(lam));
    floorv = max(chiAbsFloor, chiRelFloor * max(smax, eps));
    lam2   = sign(lam) .* max(abs(lam), floorv);
    ChiLoc = V * diag(lam2) * V.';

    for a = 1:Nc
        for b = 1:Nc
            ChiReg{a,b}(i) = ChiLoc(a,b);
        end
    end
end

% ---- precompute B = C^o - ChiReg * mu^o ----
B = cell(1,Nc);
for a = 1:Nc
    B{a} = co{a};
    for b = 1:Nc
        B{a} = B{a} - ChiReg{a,b} .* muo{b};
    end
end

% Constant operator coefficients
dC  =  2*M/dx2;
dL  = -1*M/dx2;
dR  = -1*M/dx2;

aC  =  1/dt + M*kappa*( 6)/dx4;
aL  =            M*kappa*(-4)/dx4;
aR  =            M*kappa*(-4)/dx4;
aL2 =            M*kappa*( 1)/dx4;
aR2 =            M*kappa*( 1)/dx4;

for i = 1:n
    base = Nd*(i-1);

    im1 = i - 1; if im1 < 1, im1 = n; end
    ip1 = i + 1; if ip1 > n, ip1 = 1; end
    im2 = i - 2; if im2 < 1, im2 = im2 + n; end
    ip2 = i + 2; if ip2 > n, ip2 = ip2 - n; end

    % --------------------------------------------------------
    % Mu-only equations after eliminating C
    % --------------------------------------------------------
    for a = 1:Nc
        id_m = base + a;

        % RHS = co/dt - s - A_C * B
        R(id_m) = co{a}(i)/dt - s{a}(i) ...
                - aC  * B{a}(i) ...
                - aL  * B{a}(im1) ...
                - aR  * B{a}(ip1) ...
                - aL2 * B{a}(im2) ...
                - aR2 * B{a}(ip2);

        % Diffusion part on mu_a: -M*lap(mu_a)
        L(id_m, idf(id_m      )) = L(id_m, idf(id_m      )) + dC;
        L(id_m, idf(id_m - Nd )) = L(id_m, idf(id_m - Nd )) + dL;
        L(id_m, idf(id_m + Nd )) = L(id_m, idf(id_m + Nd )) + dR;

        % A_C * Chi * mu contribution
        for b = 1:Nc
            id_b = base + b;

            L(id_m, idf(id_b        )) = L(id_m, idf(id_b        )) + aC  * ChiReg{a,b}(i);
            L(id_m, idf(id_b - Nd   )) = L(id_m, idf(id_b - Nd   )) + aL  * ChiReg{a,b}(im1);
            L(id_m, idf(id_b + Nd   )) = L(id_m, idf(id_b + Nd   )) + aR  * ChiReg{a,b}(ip1);
            L(id_m, idf(id_b - 2*Nd )) = L(id_m, idf(id_b - 2*Nd )) + aL2 * ChiReg{a,b}(im2);
            L(id_m, idf(id_b + 2*Nd )) = L(id_m, idf(id_b + 2*Nd )) + aR2 * ChiReg{a,b}(ip2);
        end
    end
end

% Solve mu-only system
sol = L \ R;

% Unpack mu
mu = cell(1,Nc);
for a = 1:Nc
    mu{a} = sol(a:Nd:end).';
end

% Reconstruct C from closure:
%   C^n = B + ChiReg * mu^n
C = cell(1,Nc);
for a = 1:Nc
    C{a} = B{a};
    for b = 1:Nc
        C{a} = C{a} + ChiReg{a,b} .* mu{b};
    end
end

% Mass normalization:
% For periodic BC, diffusion and kappa terms conserve total mass.
% Therefore:
%   mean(C^n) = mean(C^o) - dt * mean(s)
for a = 1:Nc
    mean_target = mean(co{a}) - dt * mean(s{a});
    mean_now    = mean(C{a});
    C{a}        = C{a} + (mean_target - mean_now);
end

end


function [phi] = Solver_1D_AllenCahn_Periodic(S, phio, LK, dx, dt)
% 1D Allen–Cahn with periodic boundary conditions:
%   (phi^n - phi^o)/dt = d/dx( LK * dphi^n/dx ) + S
%
% LK can be:
%   - scalar
%   - vector of length n

n   = length(phio);
L   = sparse(n,n);
R   = zeros(n,1);
dx2 = dx*dx;

% make row vectors
phio = reshape(phio, 1, []);
S    = reshape(S,    1, []);

% normalize LK
if isscalar(LK)
    LK = LK * ones(1,n);
else
    LK = reshape(LK, 1, []);
    if numel(LK) ~= n
        error('LK must be a scalar or a vector of length n.');
    end
end

for i = 1:n
    % periodic neighbors
    im1 = i - 1;  if im1 < 1, im1 = n; end
    ip1 = i + 1;  if ip1 > n, ip1 = 1; end

    % face mobilities / coefficients
    LK_L = 0.5 * (LK(i) + LK(im1));
    LK_R = 0.5 * (LK(i) + LK(ip1));

    % implicit discretization:
    % phi^n/dt - [ LK_R*(phi_{i+1}-phi_i) - LK_L*(phi_i-phi_{i-1}) ] / dx^2
    %   = phi^o/dt + S
    L(i,i)   = 1/dt + (LK_L + LK_R)/dx2;
    L(i,im1) = -LK_L/dx2;
    L(i,ip1) = -LK_R/dx2;

    R(i)     = S(i) + phio(i)/dt;
end

sol = L \ R;
phi = sol.';
end

