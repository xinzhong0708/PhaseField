function [STATE] = LE_Run_Mode_New(STATE,PARAM,MODEL)
%This function perform local equilibrium
%Two modes are available
%1) 'LE': Classical KKS type fixed E mode with fixed p and E
%2) 'GP': Grand potential based local equilibrium with fixed mu_e

%Choose the mode
mode       =  PARAM.LE_mode;

%p mode on and off
p_on       =  1e-3;
p_off      =  2e-4;

%Maximal number of active phases for LE/GP
Pmax       =  3;

%Local equilibrium parameters
alpha_LE   = [0.7 0.3 0.3 0.2];
iter_LE    = [100 100 100 100];
iter_GP    = [100 100 100 100];

%Interface excess damping
wscale_p0  =  0.60;
wscale_p1  =  0.99;
wscale_gam =  1;

%Unpack fields and prepare the calculation
ny         =  size(STATE.p,1);
nx         =  size(STATE.p,2);
N          =  ny*nx;
p          =  UnpackP(STATE.p);
c          =  UnpackC(STATE.c);
E          =  UnpackCell(STATE.E);
mu_e       =  UnpackCell(STATE.mu_e);
chi        =  UnpackChi(STATE.chi);
eta        =  reshape(PARAM.eta,1,[]);
pars       =  MODEL.pars;

%Phase index of the grains
phase_index= MODEL.phase_index(:).';

% Collapse repeated grains into thermodynamic phases
[pars,c,p,grain_to_phase] = CollapsePhases(pars,c,p,phase_index);

%Phase number after collapse
Np         =  size(p,3);
Ne         =  numel(E);
p2         =  reshape(p,N,Np);

% Damping of the interface excess energy
pars_orig      = pars;
if isfield(PARAM,'use_WScale') && PARAM.use_WScale
    pars_inter = Apply_WScale_FromP(pars,p,wscale_p0,wscale_p1,wscale_gam);
else
    pars_inter = pars;
end

%Find the active phase with hysteresis
if isfield(STATE,'LE_state') && isfield(STATE.LE_state,'active') && isequal(size(STATE.LE_state.active),[N,Np])
    active_old = STATE.LE_state.active;
else
    active_old = false(N,Np);
end
active     = (active_old & (p2 > p_off)) | (~active_old & (p2 > p_on));

%Make sure every grid has at least one active phase
for i = 1:N
    if ~any(active(i,:))
        [~,im] = max(p2(i,:));
        active(i,im) = true;
    end
end

%Limit maximum number of active phases
active  =  Limit_Pmax(active,p2,N,Np,Pmax);

%Group grids by active phase set, this will be doing for LE or GP
[active_sets,~,set_id] = unique(active,'rows');

%MAIN PART OF THE LE
for iset = 1:size(active_sets,1)

    %Active phase
    ph_act       = find(active_sets(iset,:));
    k            = numel(ph_act);
    mask         = (set_id == iset).';

    %Slice active local problem
    if PARAM.use_WScale && k > 1
        pars_loc = SliceParsWScale(pars_inter,ph_act,mask);
    else
        pars_loc = pars(ph_act);
    end
    c_loc        = SliceC(c,ph_act,mask);
    p_loc        = SliceP(p,mask,ph_act);
    eta_loc      = eta(mask);

    %Choose damping based on number of active phases
    kk           = min(k,numel(alpha_LE));
    aa           = alpha_LE(kk);

    % %MODE 1: KKS TYPE LOCAL EQUILIBRIUM WITH E AND P FIXED
    % if strcmpi(mode,'LE')
    %     mm                     =  iter_LE(kk);
    %     [~,~,E_loc]            =  SliceActiveMassBalance(pars_orig,c,p,E,ph_act,mask);
    %     [c_loc,mu_loc,chi_loc] =  LE_Calculator(pars_loc,p_loc,c_loc,E_loc,eta_loc,[aa,mm]);
    % end 
    % 
    % %MODE 2: GRAND POTENTIAL (GP) WITH MU_E FIXED
    % if strcmpi(mode,'GP')
    %     [~,~,E_loc]            =  SliceActiveMassBalance(pars_orig,c,p,E,ph_act,mask);
    %     %For one phase, use LE because E is the conserved variable.
    %     %For multiphase regions, use mass-consistent GP as smooth closure.
    %     if k == 1
    %         mm                     =  iter_LE(kk);
    %         [c_loc,mu_loc,chi_loc] =  LE_Calculator(pars_loc,p_loc,c_loc,E_loc,eta_loc,[aa,mm]);
    %     else
    %         mm                     =  iter_GP(kk);
    %         mu_loc                 =  SliceCell(mu_e,mask);
    %         [c_loc,chi_loc]        =  GP_Calculator(pars_loc,p_loc,c_loc,mu_loc,E_loc,eta_loc,aa,mm);
    %     end
    % end



    [~,~,E_loc] = SliceActiveMassBalance(pars_orig,c,p,E,ph_act,mask);

    %Pure phase: use LE because E is conserved and one-phase closure is well-defined
    if k == 1
        mm                     =  iter_LE(kk);
        [c_loc,mu_loc,chi_loc] =  LE_Calculator(pars_loc,p_loc,c_loc,E_loc,eta_loc,[aa,mm]);

        %Coexisting grid: use mass-consistent GP to avoid branch-jumpy interface LE
    else
        mm                     =  iter_GP(kk);
        mu_loc                 =  SliceCell(mu_e,mask);
        [c_loc,chi_loc]        =  GP_Calculator(pars_loc,p_loc,c_loc,mu_loc,E_loc,eta_loc,aa,mm);
    end


    %Assign back
    c    = AssignC(c,ph_act,mask,c_loc);
    mu_e = AssignCell(mu_e,mask,mu_loc);
    chi  = AssignChi(chi,mask,chi_loc);

end

% -------------------------------------------------------------------------
% Pack output
% -------------------------------------------------------------------------
c_col       = PackC(c,ny);
E           = PackCell(E,ny);
mu_e        = PackCell(mu_e,ny);
chi         = PackChi(chi,ny);

% Always calculate e from original thermodynamics
e_col       = Calc_e(pars_orig,c_col);

% Raw omega from full original thermodynamics
omg_raw_col = CalcOmegaLocal(pars_orig,c,e_col,mu_e,ny,nx,Ne,Np);
omg_col     = omg_raw_col;

% Expand thermodynamic phases back to grains
Ngrain      = numel(grain_to_phase);
c_out       = cell(1,Ngrain);
omg         = zeros(ny,nx,Ngrain);
omg_raw     = zeros(ny,nx,Ngrain);

for ig = 1:Ngrain
    iph             = grain_to_phase(ig);
    c_out{ig}       = c_col{iph};
    omg(:,:,ig)     = omg_col(:,:,iph);
    omg_raw(:,:,ig) = omg_raw_col(:,:,iph);
end

% Save active-set memory
LE_state                 =  struct();
LE_state.active          =  active;
LE_state.phase_index     =  phase_index;
LE_state.grain_to_phase  =  grain_to_phase;

% Update STATE
STATE.c        = c_out;
STATE.e        = Calc_e(MODEL.pars,c_out);
STATE.E        = E;
STATE.mu_e     = mu_e;
STATE.chi      = chi;
STATE.omg      = omg;
STATE.omg_raw  = omg_raw;
STATE.LE_state = LE_state;
end


% ALL HELPER FUNCTIONS BELOW
% =========================================================================
% Omega
% =========================================================================
function omg_col = CalcOmegaLocal(pars,c,e_col,mu_e,ny,nx,Ne,Np)
omg_col = zeros(ny,nx,Np);
for ip = 1:Np
    omg_col(:,:,ip) = reshape(PhaseG(pars{ip},c{ip}),ny,[]);
    for ie = 1:Ne
        omg_col(:,:,ip) = omg_col(:,:,ip) - e_col{ip}{ie}.*mu_e{ie};
    end
end
end


% =========================================================================
% WScale slicing
% =========================================================================
function pars_cur = SliceParsWScale(pars_inter,ph_act,mask)
%SLICEPARSWSCALE Extract nodewise w_scale for one active-set block.
pars_cur = pars_inter(ph_act);
for ia = 1:numel(pars_cur)

    if isfield(pars_cur{ia},'w_scale') && ~isempty(pars_cur{ia}.w_scale)
        ws = reshape(pars_cur{ia}.w_scale,1,[]);
        if isscalar(ws)
            pars_cur{ia}.w_scale = ws;

        elseif numel(ws) == numel(mask)
            pars_cur{ia}.w_scale = ws(mask);

        elseif numel(ws) == nnz(mask)
            pars_cur{ia}.w_scale = ws;
        end
    end
end
end

% =========================================================================
% Pack / unpack
% =========================================================================
function p = UnpackP(p)
[ny,nx,np] = size(p);
p          = reshape(p,1,ny*nx,np);
end

function c = UnpackC(c)
for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        c{ip}{ic} = reshape(c{ip}{ic},1,[]);
    end
end
end

function f = UnpackCell(f)
for i = 1:numel(f)
    f{i} = reshape(f{i},1,[]);
end
end


function ch = UnpackChi(ch)
for i = 1:size(ch,1)
    for j = 1:size(ch,2)
        ch{i,j} = reshape(ch{i,j},1,[]);
    end
end
end

function c = PackC(c,ny)

for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        c{ip}{ic} = reshape(c{ip}{ic},ny,[]);
    end
end

end


function f = PackCell(f,ny)

for i = 1:numel(f)
    f{i} = reshape(f{i},ny,[]);
end

end


function ch = PackChi(ch,ny)

for i = 1:size(ch,1)
    for j = 1:size(ch,2)
        ch{i,j} = reshape(ch{i,j},ny,[]);
    end
end

end


% =========================================================================
% Collapse grains
% =========================================================================
function [pars_c,c_c,p_c,grain_to_phase] = CollapsePhases(pars,c,p,phase_index)

phase_id       = unique(phase_index,'stable');
Ngrain         = numel(c);
Nphase         = numel(phase_id);
N              = size(p,2);
grain_to_phase = zeros(1,Ngrain);

for iph = 1:Nphase
    grain_to_phase(phase_index == phase_id(iph)) = iph;
end

pars_c = cell(1,Nphase);
p_c    = zeros(1,N,Nphase);
c_c    = cell(1,Nphase);

for iph = 1:Nphase
    grains      = find(grain_to_phase == iph);
    ig0         = grains(1);
    pars_c{iph} = pars{ig0};
    for ig = grains
        p_c(:,:,iph) = p_c(:,:,iph) + p(:,:,ig);
    end
    Nc          = numel(c{ig0});
    den         = reshape(p_c(:,:,iph),1,N);
    good        = den > eps;
    c_c{iph}    = cell(1,Nc);

    for ic = 1:Nc
        num     = zeros(1,N);
        for ig  = grains
            num = num + reshape(p(:,:,ig),1,N).*reshape(c{ig}{ic},1,N);
        end
        tmp = reshape(c{ig0}{ic},1,N);
        tmp(good) = num(good)./den(good);
        c_c{iph}{ic} = tmp;
    end
end

end


%Slice active phases and subtract inactive phase contribution from E
function [c_blk,p_blk,E_blk] = SliceActiveMassBalance(pars,c,p,E,ph_act,mask)
Np    = numel(c);
Ne    = numel(E);
Nloc  = nnz(mask);
c_blk = SliceC(c,ph_act,mask);
p_blk = SliceP(p,mask,ph_act);
E_blk = SliceCell(E,mask);

inactive = setdiff(1:Np,ph_act);
Efix     = zeros(Ne,Nloc);

for ii = 1:numel(inactive)

    ip    = inactive(ii);
    c_tmp = SliceC(c,ip,mask);
    c_ip  = c_tmp{1};

    R     = PhaseThermo(pars{ip},c_ip);
    p_ip  = reshape(p(:,mask,ip),1,Nloc);
    e_ip  = StackFields(R.e,Nloc);

    Efix  = Efix + e_ip.*p_ip;
end

for ie = 1:Ne
    E_blk{ie} = E_blk{ie} - Efix(ie,:);
end

end


function c2 = SliceC(c,ph,mask)
c2 = cell(1,numel(ph));
for a = 1:numel(ph)
    ip    = ph(a);
    c2{a} = cell(size(c{ip}));
    for ic = 1:numel(c{ip})
        x          = reshape(c{ip}{ic},1,[]);
        c2{a}{ic} = x(mask);
    end
end
end

function p2 = SliceP(p,mask,ph)
p2 = p(:,mask,ph);
p2 = reshape(p2,1,nnz(mask),numel(ph));
end

function f2 = SliceCell(f,mask)
f2 = cell(size(f));
for i = 1:numel(f)
    x     = reshape(f{i},1,[]);
    f2{i}= x(mask);
end
end


function c = AssignC(c,ph,mask,c_sub)
for a = 1:numel(ph)
    ip = ph(a);
    for ic = 1:numel(c{ip})
        x              = reshape(c{ip}{ic},1,[]);
        x(mask)        = reshape(c_sub{a}{ic},1,[]);
        c{ip}{ic}      = x;
    end
end
end


function f = AssignCell(f,mask,f_sub)
for i = 1:numel(f)
    x       = reshape(f{i},1,[]);
    x(mask) = reshape(f_sub{i},1,[]);
    f{i}    = x;
end
end

function ch = AssignChi(ch,mask,ch_sub)
Ne = size(ch,1);
for i = 1:Ne
    for j = 1:Ne
        x          = reshape(ch{i,j},1,[]);
        x(mask)    = reshape(ch_sub{i,j},1,[]);
        ch{i,j}    = x;
    end
end
end

function A = StackFields(fields,N)
A = zeros(numel(fields),N);
for i = 1:numel(fields)
    A(i,:) = reshape(fields{i},1,N);
end
end

function [active] = Limit_Pmax(active,p2,N,Np,Pmax)
nact = sum(active,2);
bad  = nact > Pmax;
if any(bad)
    [~,ord]       = sort(p2,2,'descend');
    keep          = false(N,Np);
    row           = repmat((1:N)',1,Pmax);
    col           = ord(:,1:Pmax);
    ind           = sub2ind([N,Np],row,col);
    keep(ind)     = true;
    active(bad,:) = keep(bad,:);
end
end


% =========================================================================
% GP calculator
% =========================================================================
function [c,chi] = GP_Calculator(pars,p,c,mu_e,E_in,eta,alpha,Miter)
c_tol = 1e-6;
wmu   = 10;
lam   = 1e-8;
for it = 1:Miter
    dc = GP_Step_MassConsistent(pars,p,c,mu_e,E_in,eta,lam,wmu);
    dc_step = ScaleStep(dc,alpha);
    if MaxStep(dc_step) < c_tol
        break
    end
    c  = AddStep(c,dc,alpha);
end
if it == Miter
    disp('GP fails'); 
end
Ne  = numel(mu_e);
chi = Chi_FromMu_Projected(pars,p,c,eta,Ne,lam,wmu);
end

function dc = ScaleStep(dc,alpha)
for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        dc{ip}{ic} = alpha.*dc{ip}{ic};
    end
end
end

% =========================================================================
% GP step: fixed mu_e
% =========================================================================
function dc = GP_Step_MassConsistent(pars,p,c,mu_e,E_in,eta,lam,wmu)
Np     = numel(c);
Ne     = numel(mu_e);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);
E_mat  = StackFields(E_in,N);
eta_v  = EtaVector(eta,N);
dc     = ZeroStepLike(c);
Emix   = zeros(Ne,N);
Rall   = cell(1,Np);
%Current mixture composition
for ip = 1:Np
    Rall{ip} = PhaseThermo(pars{ip},c{ip});
    p_ip     = reshape(p(:,:,ip),1,N);
    Emix     = Emix + StackFields(Rall{ip}.e,N).*p_ip;
end
%Mass residual: E = Emix + mu/eta
R_E_global = E_mat - Emix - mu_mat./eta_v;

%Residual distribution weight
psum2 = zeros(1,N);
for jp = 1:Np
    p_j   = reshape(p(:,:,jp),1,N);
    psum2 = psum2 + p_j.^2;
end
psum2 = max(psum2,eps);

for ip = 1:Np
    R = Rall{ip};
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end
    J    = R.Jac;
    H    = SymPages(R.H_c);
    JT   = permute(J,[2 1 3]);
    HT   = permute(H,[2 1 3]);
    Nc   = size(H,1);
    mu_c = StackFields(R.mu_c,N);
    %Chemical-potential residual: %   mu_c = J' * mu_e
    R_mu = reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),Nc,N) - mu_c;
    %This phase contributes to the mass residual in proportion to p_ip
    p_ip = reshape(p(:,:,ip),1,N);
    R_E  = R_E_global.*p_ip./psum2;
    %Solve mixed local projection: J dc  ~= R_E  and   H dc  ~= R_mu
    A = pagemtimes(JT,J) + wmu*pagemtimes(HT,H) + repmat(eye(Nc),1,1,N).*lam;
    b = pagemtimes(JT,reshape(R_E,Ne,1,N)) + wmu*pagemtimes(HT,reshape(R_mu,Nc,1,N));
    d = reshape(pagemldivide(A,b),Nc,N);
    for ic = 1:Nc
        dc{ip}{ic} = d(ic,:);
    end
end
end


% =========================================================================
% GP susceptibility: chi = dE/dmu_e = sum_ip p_ip * de_ip/dmu_e + I/eta
% =========================================================================
function chi = Chi_FromMu_Projected(pars,p,c,eta,Ne,lam,wmu)
Np    = numel(c);
N     = numel(c{1}{1});
eta_v = EtaVector(eta,N);
Cmix  = zeros(Ne,Ne,N);
for ip = 1:Np
    R    = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);
    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        Cphase = zeros(Ne,Ne,N);
    else
        J      = R.Jac;
        H      = SymPages(R.H_c);
        JT     = permute(J,[2 1 3]);
        HT     = permute(H,[2 1 3]);
        Nc     = size(H,1);
        A      = pagemtimes(JT,J) + wmu*pagemtimes(HT,H) + repmat(eye(Nc),1,1,N).*lam;
        B      = wmu*pagemtimes(HT,JT);
        Cphase = pagemtimes(J,pagemldivide(A,B));
    end
    Cmix = Cmix + Cphase.*reshape(p_ip,1,1,N);
end
I   = repmat(eye(Ne),1,1,N);
chi = PagesToChi(Cmix + I.*reshape(1./eta_v,1,1,N),Ne);
end


% =========================================================================
% GP helpers
% =========================================================================
function H = SymPages(H)
H = 0.5*(H + permute(H,[2 1 3]));
end


function dc = ZeroStepLike(c)
dc = c;
for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        dc{ip}{ic} = zeros(size(c{ip}{ic}));
    end
end
end


function m = MaxStep(dc)
m = 0;
for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        m = max(m,max(abs(dc{ip}{ic})));
    end
end
end


function c_new = AddStep(c,dc,alpha)
c_new = c;
for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        c_new{ip}{ic} = c{ip}{ic} + alpha.*dc{ip}{ic};
    end
end
end

function eta_v = EtaVector(eta,N)
if isscalar(eta)
    eta_v = eta*ones(1,N);
else
    eta_v = reshape(eta,1,N);
end
end

function chi = PagesToChi(A,Ne)
chi = cell(Ne,Ne);
for i = 1:Ne
    for j = 1:Ne
        chi{i,j} = reshape(A(i,j,:),1,[]);
    end
end
end