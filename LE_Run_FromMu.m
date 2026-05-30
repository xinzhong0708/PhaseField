function [STATE,DIAG] = LE_Run_FromMu(STATE,PARAM,MODEL)
%LE_RUN_FROMMU Nonlinear thermodynamic projection at prescribed mu_e and p.
%
% This is the mu-primary counterpart of LE_Run. It retains the incoming
% STATE.mu_e and STATE.p, solves active phase compositions from
%
%       d/dc [ g_phase(c) - mu_e' * e_phase(c) ] = 0,
%
% then reconstructs the finite-eta bulk composition and susceptibility:
%
%       E   = sum_phase p_le * e_phase(mu_e) + mu_e / eta
%       chi = sum_phase p_le * de_phase/dmu_e + I / eta.
%
% Minimal changes relative to LE_Run (PFM_Update3):
%   - Grain collapse, p_th smoothing, hysteretic active sets, Pmax,
%     W-scaling option, unpack/pack logic and omega calculation are retained.
%   - The fixed-E call LE_Calculator(...) is replaced by the local helper
%     LE_Calculator_FromMu(...), which does not change mu_e.
%   - STATE.E is reconstructed from the prescribed mu_e after projection.
%
% IMPORTANT:
%   - This function does NOT enforce the spatial mean of E. If required,
%     correct the mean outside this function through a spatially uniform
%     correction to mu_e and call this function again.
%   - For strict comparison with the current LE_Run, this function uses p_le
%     (not geometric p) in E and chi reconstruction, because current LE_Run
%     uses p_le in its inverse constitutive solve.
%
% Usage:
%       [STATE,DIAG] = LE_Run_FromMu(STATE,PARAM,MODEL);
%
% DIAG.max_stationarity measures the residual of
%       mu_c - J' * mu_e
% for active solution phases after the projection.

%Unpack structured variables
pars        = MODEL.pars;
phase_index = MODEL.phase_index(:).';

p           = STATE.p;
c           = STATE.c;
E_in        = STATE.E; %#ok<NASGU> %Retained only for diagnostics outside this function
mu_e        = STATE.mu_e;
chi         = STATE.chi;
eta_vec     = PARAM.eta(:);

if isfield(STATE,'LE_state') && ~isempty(STATE.LE_state)
    LE_state = STATE.LE_state;
else
    LE_state = struct();
end

%Reshape 2D fields into 1D LE format
ny          = size(p,1);
nx          = size(p,2); %#ok<NASGU>

c           = Unpack_c(c);
p           = Unpack_p(p);
mu_e        = Unpack_E(mu_e);
chi         = Unpack_Chi(chi);

%Collapse repeated grains to thermodynamic phases
[pars,c,p,grain_to_phase] = Collapse_LE_Phases(pars,c,p,phase_index);

%Local-equilibrium controls: unchanged from LE_Run
Pmax        = 4;

p_tail      = 1e-3;
p_full      = 5e-2;
p_on        = 5e-3;
p_off       = 2e-3;       %Must be smaller than p_on for hysteresis

%Damping factors for 1-, 2-, 3- and 4-phase coexistence
alpha_LE    = [0.5  0.5   0.4   0.3 ];
iter_LE     = [100  1000  1000  1000];

Np          = numel(c);
N           = numel(c{1}{1});

%Thermodynamic weights before active-set filtering
p_th        = Calc_Thermo_p(p,p_tail,p_full);
pThMat      = reshape(p_th,N,Np);

%Hysteretic active set: unchanged from LE_Run
if isempty(LE_state) || ~isfield(LE_state,'active') || ~isequal(size(LE_state.active),[N,Np])
    active_old = false(N,Np);
else
    active_old = LE_state.active;
end

active = ( active_old & pThMat > p_off) | (~active_old & pThMat > p_on );

%Ensure at least one and at most Pmax active phases at every node
repair = find(~any(active,2) | sum(active,2) > Pmax).';
for i = repair
    if ~any(active(i,:))
        [~,idmax]       = max(pThMat(i,:));
        active(i,idmax) = true;
    end
    if sum(active(i,:)) > Pmax
        score           = pThMat(i,:) + 0.5*p_on*active_old(i,:);
        [~,ord]         = sort(score,'descend');
        active(i,:)     = false;
        active(i,ord(1:Pmax)) = true;
    end
end

%Renormalize LE weights after inactive tails are removed
p_le            = p_th .* reshape(active,1,N,Np);
p_le            = p_le ./ max(sum(p_le,3),eps);

LE_state.active = active;
LE_state.p_th   = p_th;
LE_state.p_le   = p_le;

%Optional interface W scaling, evaluated only when used
if PARAM.use_WScale
    pars_inter = Apply_WScale_FromP(pars,p,0.85,0.99,1);
end

%Diagnostics in unpacked node ordering
stat_node = zeros(1,N);
iter_node = zeros(1,N);

%Process only phase combinations present in the grid
[active_sets,~,set_id] = unique(active,'rows');

for iset = 1:size(active_sets,1)

    ph_act = find(active_sets(iset,:));
    k      = numel(ph_act);
    mask   = (set_id == iset).';

    c_cur  = Slice_c(c,ph_act,mask);
    mu_cur = Slice_E(mu_e,mask);       %Prescribed coupled-solver mu_e
    E_cur  = Slice_E(E_in,mask);       %Prescribed coupled-solver mu_e
    eta    = eta_vec(mask);
    p_cur  = p_le(:,mask,ph_act);

    if PARAM.use_WScale && k > 1
        pars_cur = Slice_Pars_WScale(pars_inter,ph_act,mask);
    else
        pars_cur = pars(ph_act);
    end

    [c_tmp,chi_tmp,stat_tmp,niter_tmp] = LE_Calculator_FromMu( ...
        pars_cur,p_cur,c_cur,mu_cur,eta,E_cur,[alpha_LE(k),iter_LE(k)]);

    %Reuse current unpack/assignment route; mu_cur is written back unchanged.
    [c,mu_e,chi] = Assign_LE_Back(c,mu_e,chi,c_tmp,mu_cur,chi_tmp,ph_act,mask);

    stat_node(mask) = stat_tmp;
    iter_node(mask) = niter_tmp;

end

%Forward finite-eta closure from prescribed mu_e and active p_le.
%E is reconstructed here; it is not supplied to the local composition solve.
E       = Calc_E_FromMu(pars,c,p_le,mu_e,eta_vec);

%Pack collapsed thermodynamic phases
c_col   = Pack_c(c,ny);
E       = Pack_E(E,ny);
mu_e    = Pack_E(mu_e,ny);
chi     = Pack_chi(chi,ny);
e_col   = Calc_e(pars,c_col);

%Calculate grand potentials from original excess-energy thermodynamics.
%This retains the behaviour of LE_Run when PARAM.use_WScale is active.
Ne      = numel(mu_e);
omg_col = zeros(ny,size(STATE.p,2),Np);

for ip = 1:Np
    omg_col(:,:,ip) = reshape(PhaseG(pars{ip},c{ip}),ny,[]);
    for ie = 1:Ne
        omg_col(:,:,ip) = omg_col(:,:,ip) - e_col{ip}{ie} .* mu_e{ie};
    end
end

%Copy collapsed results back to grain-resolved variables
Ngrain  = numel(grain_to_phase);
c_out   = cell(1,Ngrain);
omg     = zeros(ny,size(STATE.p,2),Ngrain);

for ig = 1:Ngrain
    iph         = grain_to_phase(ig);
    c_out{ig}   = c_col{iph};
    omg(:,:,ig) = omg_col(:,:,iph);
end

%Store collapsed active-set information
LE_state.phase_index    = phase_index;
LE_state.grain_to_phase = grain_to_phase;

%Write back to STATE. STATE.p and STATE.phi are not changed.
STATE.c        = c_out;
STATE.e        = Calc_e(MODEL.pars,c_out);
STATE.E        = E;
STATE.mu_e     = mu_e;       %Preserved input chemical potential field
STATE.chi      = chi;
STATE.omg      = omg;
STATE.LE_state = LE_state;

%Optional diagnostics
DIAG.max_stationarity = max(stat_node);
DIAG.stat_node        = reshape(stat_node,ny,[]);
DIAG.iter_max         = max(iter_node);
DIAG.iter_node        = reshape(iter_node,ny,[]);
DIAG.active           = active;
DIAG.p_le             = p_le;

end


function [c,chi,stat_node,niter] = LE_Calculator_FromMu(pars,p,c,mu_e,eta,E_in,level)
%LE_CALCULATOR_FROMMU Solve phase compositions at prescribed mu_e.
%
% The fixed-mu stationarity equation for each active phase is
%       mu_c(c) - J(c)'*mu_e = 0.
% The Newton direction is the same c-block used in LE_Calculator, but
% mu_e is no longer obtained by inversion of the tangent-reconstructed E.

%Prepare
Np           = length(c);
Ne           = length(mu_e);
N            = numel(mu_e{1});
alpha        = level(1);
Miter        = level(2);
c_tol        = 1e-5;

%Line-search & damping parameters matched to LE_Calculator
MaxLS        = 10;
amin         = 1e-7;
energy_tol   = 1e-9;
lam_c        = 0;

%Check whether any phase has internal composition degrees of freedom
has_dof      = false(1,Np);
for ip = 1:Np
    R0          = PhaseThermo(pars{ip},c{ip});
    has_dof(ip) = ~(isempty(R0.mu_c) || isempty(R0.H_c) || isempty(R0.Jac));
end

niter = 0;

%If all phases are pure/no-DOF phases, only the finite-eta tangent is needed
if ~any(has_dof)
    chi       = Calc_Chi_FromMu(pars,p,c,eta,Ne,0);
    stat_node = zeros(1,N);
    return
end

%Newton iteration at fixed prescribed mu_e
for it = 1:Miter

    niter = it;
    c_old = c;

    %Grand-potential objective at the supplied mu_e
    F_old  = LE_Mu_Objective(pars,p,c_old,mu_e);

    %Newton composition direction at fixed mu_e
    % dc_all = LE_Mu_Step(pars,c_old,mu_e,lam_c);
    lam_proj = 1e-10;
    wmu_proj = 1;
    dc_all   = LE_Mu_EProjected_Step(pars,p,c_old,mu_e,E_in,eta,lam_proj,wmu_proj);

    %If proposed step is essentially zero, stop
    dcmax = MaxAbsStep(dc_all);
    if dcmax < c_tol
        break
    end

    %Backtracking line search independently for each grid point
    good_node = false(1,N);
    alpha_try = alpha*ones(1,N);
    alpha_acc = zeros(1,N);

    for ils = 1:MaxLS

        c_try = AddStep(c_old,dc_all,alpha_try);
        F_try = LE_Mu_Objective(pars,p,c_try,mu_e);

        good = isfinite(F_try) & ...
               (F_try <= F_old + energy_tol.*max(1,abs(F_old)));

        good_new            = good & ~good_node;
        alpha_acc(good_new) = alpha_try(good_new);
        good_node(good_new) = true;

        bad            = ~good_node;
        alpha_try(bad) = 0.2*alpha_try(bad);

        if all(good_node | alpha_try < amin)
            break
        end

    end

    %Retain the LE_Calculator behaviour: nodes without accepted descent do
    %not move in this iteration; successful nodes remain vectorized.
    if any(good_node)
        c = AddStep(c_old,dc_all,alpha_acc);
    else
        c = c_old;
        break
    end

    %Check convergence of composition update
    cchg = 0;
    for ip = 1:Np
        if ~has_dof(ip)
            continue
        end
        for ic = 1:length(c{ip})
            cchg = max(cchg,max(abs(c{ip}{ic}-c_old{ip}{ic})));
        end
    end

    if cchg < c_tol
        break
    end

end

%Constitutive tangent on the final fixed-mu composition state
% chi       = Calc_Chi_FromMu(pars,p,c,eta,Ne,0);

chi       = Calc_Chi_FromMu_EProjected(pars,p,c,eta,Ne,lam_proj,wmu_proj);

stat_node = Mu_Stationarity_Residual(pars,c,mu_e);

end


function chi = Calc_Chi_FromMu_EProjected(pars,p,c,eta,Ne,lam_proj,wmu_proj)
% Effective bounded susceptibility for E-projected fixed-mu LE.
% Avoids direct H^{-1}. Uses:
% dc/dmu = A^{-1} * wmu * H' * J'
% A = J'J + wmu H'H + lam I
%
% de/dmu = J dc/dmu

Np = length(c);
N  = numel(c{1}{1});

if isscalar(eta)
    eta_vec = eta*ones(1,N);
else
    eta_vec = eta(:).';
end

Cmix = zeros(Ne,Ne,N);

for ip = 1:Np

    R    = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);

    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        C_phase = zeros(Ne,Ne,N);
    else
        J  = R.Jac;
        H  = 0.5*(R.H_c + permute(R.H_c,[2 1 3]));
        Nc = size(H,1);

        JT  = permute(J,[2 1 3]);
        HT  = permute(H,[2 1 3]);

        JTJ = pagemtimes(JT,J);
        HTH = pagemtimes(HT,H);

        A = JTJ + wmu_proj*HTH + repmat(eye(Nc),1,1,N).*lam_proj;

        % dc/dmu = A^{-1} * wmu * H' * J'
        B = wmu_proj * pagemtimes(HT,JT);
        dc_dmu = pagemldivide(A,B);

        % de/dmu = J * dc/dmu
        C_phase = pagemtimes(J,dc_dmu);
    end

    Cmix = Cmix + C_phase .* reshape(p_ip,1,1,N);
end

IpagesE  = repmat(eye(Ne),1,1,N);
chi_page = Cmix + IpagesE .* reshape(1./eta_vec,1,1,N);

chi = cell(Ne,Ne);
for i = 1:Ne
    for j = 1:Ne
        chi{i,j} = reshape(chi_page(i,j,:),1,[]);
    end
end

end
function dc_all = LE_Mu_Step(pars,c,mu_e,lam_c)
%LE_MU_STEP Vectorized Newton direction at prescribed chemical potential.

Np     = length(c);
Ne     = length(mu_e);
N      = numel(mu_e{1});
mu_mat = cell2mat(mu_e(:));

%Initialize step
dc_all = c;
for ip = 1:Np
    for ic = 1:length(c{ip})
        dc_all{ip}{ic} = zeros(size(c{ip}{ic}));
    end
end

for ip = 1:Np

    R = PhaseThermo(pars{ip},c{ip});

    %Pure phase: no internal composition variable to adjust
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    mu_c   = cell2mat(R.mu_c(:));
    Hc     = R.H_c;
    J      = R.Jac;
    Nc     = size(mu_c,1);
    Hreg   = RegularizeHessian(Hc,lam_c);
    JT     = permute(J,[2 1 3]);
    JTmu3  = pagemtimes(JT,reshape(mu_mat,Ne,1,N));
    JTmu   = reshape(JTmu3,Nc,N);
    rhs_dc = JTmu - mu_c;
    dc3    = pagemldivide(Hreg,reshape(rhs_dc,Nc,1,N));
    dc     = reshape(dc3,Nc,N);

    for ic = 1:length(c{ip})
        dc_all{ip}{ic} = dc(ic,:);
    end

end

end


function chi = Calc_Chi_FromMu(pars,p,c,eta,Ne,lam_c)
%CALC_CHI_FROMMU Compute chi = sum p*J*H^{-1}*J' + I/eta.
% This duplicates the final tangent definition of LE_Calculator so the
% coupled AC-CH solver receives the same kind of stabilized susceptibility.

Np = length(c);
N  = numel(c{1}{1});

if isscalar(eta)
    eta_vec = eta*ones(1,N);
else
    eta_vec = eta(:).';
end

Cmix = zeros(Ne,Ne,N);

for ip = 1:Np

    R    = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);

    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        C_phase = zeros(Ne,Ne,N);
    else
        Hreg     = RegularizeHessian(R.H_c,lam_c);
        JT       = permute(R.Jac,[2 1 3]);
        Hinv_JT  = pagemldivide(Hreg,JT);
        C_phase  = pagemtimes(R.Jac,Hinv_JT);
    end

    Cmix = Cmix + C_phase .* reshape(p_ip,1,1,N);

end

IpagesE  = repmat(eye(Ne),1,1,N);
chi_page = Cmix + IpagesE .* reshape(1./eta_vec,1,1,N);
chi      = cell(Ne,Ne);

for i = 1:Ne
    for j = 1:Ne
        chi{i,j} = reshape(chi_page(i,j,:),1,[]);
    end
end

end


function E = Calc_E_FromMu(pars,c,p,mu_e,eta)
%CALC_E_FROMMU Forward finite-eta constitutive reconstruction:
%       E = sum_ip p_ip*e_ip(c_ip(mu_e)) + mu_e/eta.

Np = length(c);
Ne = length(mu_e);
N  = numel(mu_e{1});

if isscalar(eta)
    eta_vec = eta*ones(1,N);
else
    eta_vec = eta(:).';
end

Emix = zeros(Ne,N);

for ip = 1:Np
    R     = PhaseThermo(pars{ip},c{ip});
    e_ip  = cell2mat(R.e(:));
    p_ip  = reshape(p(:,:,ip),1,N);
    Emix  = Emix + e_ip .* p_ip;
end

E = cell(1,Ne);
for ie = 1:Ne
    E{ie} = Emix(ie,:) + reshape(mu_e{ie},1,N)./eta_vec;
end

end


function stat_node = Mu_Stationarity_Residual(pars,c,mu_e)
%MU_STATIONARITY_RESIDUAL max |mu_c - J'*mu_e| over active phases.

Np        = length(c);
Ne        = length(mu_e);
N         = numel(mu_e{1});
mu_mat    = cell2mat(mu_e(:));
stat_node = zeros(1,N);

for ip = 1:Np

    R = PhaseThermo(pars{ip},c{ip});

    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    mu_c  = cell2mat(R.mu_c(:));
    JT     = permute(R.Jac,[2 1 3]);
    JTmu3  = pagemtimes(JT,reshape(mu_mat,Ne,1,N));
    JTmu   = reshape(JTmu3,size(mu_c,1),N);
    stat_node = max(stat_node,max(abs(mu_c-JTmu),[],1));

end

end


function F = LE_Mu_Objective(pars,p,c,mu_e)
%LE_MU_OBJECTIVE Weighted local grand-potential objective at prescribed mu.
%
% The p weights preserve the same active-mixture weighting convention used
% in LE_Run. For every active phase with p > 0, the stationary composition
% is nevertheless determined by its own g - mu'*e condition.

Np     = length(c);
Ne     = length(mu_e);
N      = numel(mu_e{1});
mu_mat = cell2mat(mu_e(:));
F      = zeros(1,N);

try
    for ip = 1:Np
        R     = PhaseThermo(pars{ip},c{ip});
        g     = R.g(:).';
        e     = cell2mat(R.e(:));
        p_ip  = reshape(p(:,:,ip),1,N);
        omg   = g - sum(e .* mu_mat,1);
        F     = F + p_ip .* omg;
    end
catch
    F = inf(1,N);
end

end


function Hreg = RegularizeHessian(Hc,lam_c)
%REGULARIZEHESSIAN Same Hessian regularization used in LE_Calculator.
[Nc,~,N] = size(Hc);
Hreg     = zeros(Nc,Nc,N);
for i = 1:N
    H           = 0.5*(Hc(:,:,i)+Hc(:,:,i).');
    scale       = max(1,norm(H,'fro')/max(1,Nc));
    evmin       = min(eig(H));
    shift       = max(0,1e-20*scale-evmin);
    Hreg(:,:,i) = H + (shift+lam_c*scale)*eye(Nc);
end
end


function c_new = AddStep(c_old,dc_all,alpha_node)
%ADDSTEP Composition update with nodewise line-search coefficient.

c_new = c_old;
for ip = 1:length(c_old)
    for ic = 1:length(c_old{ip})
        c_new{ip}{ic} = c_old{ip}{ic} + alpha_node.*dc_all{ip}{ic};
    end
end

end


function dcmax = MaxAbsStep(dc_all)
%MAXABSSTEP Maximum absolute composition Newton step.

dcmax = 0;
for ip = 1:length(dc_all)
    for ic = 1:length(dc_all{ip})
        dcmax = max(dcmax,max(abs(dc_all{ip}{ic})));
    end
end

end


function pars_cur = Slice_Pars_WScale(pars_inter,ph_act,mask)
%SLICE_PARS_WSCALE Extract nodewise W scaling for one active subset.

pars_cur = pars_inter(ph_act);
for ia = 1:numel(pars_cur)
    if isfield(pars_cur{ia},'w_scale') && ~isempty(pars_cur{ia}.w_scale)
        pars_cur{ia}.w_scale = pars_cur{ia}.w_scale(mask);
    end
end

end


function p_th = Calc_Thermo_p(p,p_tail,p_full)
%CALC_THERMO_P Smooth small geometric phase tails for LE calculations.

if p_full <= p_tail
    error('p_full must be larger than p_tail')
end

x     = (p-p_tail)./(p_full-p_tail);
x     = min(max(x,0),1);

w     = p .* x.^2 .* (3-2*x);
wsum  = sum(w,3);
p_th  = w ./ max(wsum,eps);

%If all smoothed weights vanish, retain the dominant geometric phase
bad = wsum <= eps;

if any(bad(:))

    [~,idmax] = max(p,[],3);

    for ip = 1:size(p,3)
        tmp                    = p_th(:,:,ip);
        tmp(bad & idmax == ip) = 1;
        p_th(:,:,ip)           = tmp;
    end

end

end


function [pars_c,c_c,p_c,grain_to_phase] = Collapse_LE_Phases(pars,c,p,phase_index)
%COLLAPSE_LE_PHASES Sum grains of the same thermodynamic phase for LE.

phase_id       = unique(phase_index,'stable');
Ngrain         = numel(c);
Nphase         = numel(phase_id);
N              = size(p,2);
grain_to_phase = zeros(1,Ngrain);

for iph = 1:Nphase
    grain_to_phase(phase_index == phase_id(iph)) = iph;
end

%Representative thermodynamic data
pars_c = cell(1,Nphase);

for iph = 1:Nphase
    ig          = find(grain_to_phase == iph,1,'first');
    pars_c{iph} = pars{ig};
end

%Collapse p by summing grains of the same thermodynamic phase
p_c = zeros(1,N,Nphase);

for ig = 1:Ngrain
    iph          = grain_to_phase(ig);
    p_c(:,:,iph) = p_c(:,:,iph) + p(:,:,ig);
end

%Collapse c by p-weighted averaging
c_c = cell(1,Nphase);

for iph = 1:Nphase

    grains   = find(grain_to_phase == iph);
    ig0      = grains(1);
    Nc       = numel(c{ig0});
    c_c{iph} = cell(1,Nc);
    den      = reshape(sum(p(:,:,grains),3),size(c{ig0}{1}));
    good     = den > eps;

    for ic = 1:Nc

        num = zeros(size(c{ig0}{ic}));

        for ig = grains
            num = num + reshape(p(:,:,ig),size(num)).*c{ig}{ic};
        end

        tmp             = c{ig0}{ic};
        tmp(good)       = num(good)./den(good);
        c_c{iph}{ic}    = tmp;

    end

end

end
