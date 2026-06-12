function [STATE,DIAG] = LE_Run_Mode(STATE,PARAM,MODEL)
%LE_RUN_MODE Simple LE / GP / Auto local thermodynamic update.
%
% Modes:
%   PARAM.LE_mode = 'LE'   : E-primary LE_Calculator
%   PARAM.LE_mode = 'GP'   : mu-primary projected GP closure
%   PARAM.LE_mode = 'Auto' : GP where H_c is positive definite, LE otherwise
%
% Main stabilizations:
%   1. Active set is selected using p_le, but local mass balance uses
%      physical p and subtracts inactive phase contributions from E.
%   2. GP composition update is p-damped, capped, and bounded.
%   3. GP does not write conserved E unless PARAM.LE_GP_write_E = true.
%   4. Optional WScale can be used for local solve and optionally for omega.

% -------------------------------------------------------------------------
% Direct options
% -------------------------------------------------------------------------
mode = 'LE';
if isfield(PARAM,'LE_mode') && ~isempty(PARAM.LE_mode)
    mode = PARAM.LE_mode;
end
mode = lower(char(mode));

Pmax   = 4;
p_tail = 1e-4;
p_full = 1e-2;
p_on   = 1e-3;
p_off  = 5e-4;

alpha_LE = [0.7 0.4 0.3 0.2];
iter_LE  = [100 1000 1000 1000];
iter_GP  = [50 50 50 50];

h_tol     = 1e-9;
lam_gp    = 1e-9;
wmu_gp    = 100;
c_tol     = 1e-7;
write_GP_E = false;

use_WScale       = false;
use_WScale_omega = false;
wscale_p0        = 0.85;
wscale_p1        = 0.99;
wscale_gam       = 1;

gp_dc_cap   = 0.02;
gp_c_floor  = 1e-10;
gp_p_freeze = 1e-4;
gp_p_solve  = 1e-2;

% In GP mode, keep pure one-phase blocks on LE.
% This avoids fixed-mu GP projection creating dE spikes in pure regions.
gp_pure_LE = true;

if isfield(PARAM,'use_WScale'),       use_WScale       = PARAM.use_WScale;       end
if isfield(PARAM,'use_WScale_omega'), use_WScale_omega = PARAM.use_WScale_omega; end
if isfield(PARAM,'WScale_p0'),        wscale_p0        = PARAM.WScale_p0;        end
if isfield(PARAM,'WScale_p1'),        wscale_p1        = PARAM.WScale_p1;        end
if isfield(PARAM,'WScale_gam'),       wscale_gam       = PARAM.WScale_gam;       end

if isfield(PARAM,'LE_Pmax'),        Pmax        = PARAM.LE_Pmax;        end
if isfield(PARAM,'LE_p_tail'),      p_tail      = PARAM.LE_p_tail;      end
if isfield(PARAM,'LE_p_full'),      p_full      = PARAM.LE_p_full;      end
if isfield(PARAM,'LE_p_on'),        p_on        = PARAM.LE_p_on;        end
if isfield(PARAM,'LE_p_off'),       p_off       = PARAM.LE_p_off;       end
if isfield(PARAM,'LE_alpha'),       alpha_LE    = PARAM.LE_alpha;       end
if isfield(PARAM,'LE_iter_LE'),     iter_LE     = PARAM.LE_iter_LE;     end
if isfield(PARAM,'LE_iter_GP'),     iter_GP     = PARAM.LE_iter_GP;     end
if isfield(PARAM,'LE_h_tol'),       h_tol       = PARAM.LE_h_tol;       end
if isfield(PARAM,'LE_GP_lam_proj'), lam_gp      = PARAM.LE_GP_lam_proj; end
if isfield(PARAM,'LE_GP_wmu_proj'), wmu_gp      = PARAM.LE_GP_wmu_proj; end
if isfield(PARAM,'LE_GP_c_tol'),    c_tol       = PARAM.LE_GP_c_tol;    end
if isfield(PARAM,'LE_GP_write_E'),  write_GP_E  = PARAM.LE_GP_write_E;  end

if isfield(PARAM,'LE_GP_dc_cap'),   gp_dc_cap   = PARAM.LE_GP_dc_cap;   end
if isfield(PARAM,'LE_GP_c_floor'),  gp_c_floor  = PARAM.LE_GP_c_floor;  end
if isfield(PARAM,'LE_GP_p_freeze'), gp_p_freeze = PARAM.LE_GP_p_freeze; end
if isfield(PARAM,'LE_GP_p_solve'),  gp_p_solve  = PARAM.LE_GP_p_solve;  end
if isfield(PARAM,'LE_GP_pure_LE'),  gp_pure_LE  = PARAM.LE_GP_pure_LE;  end

% -------------------------------------------------------------------------
% Unpack fields to 1-by-N form
% -------------------------------------------------------------------------
ny = size(STATE.p,1);
nx = size(STATE.p,2);

p    = UnpackP(STATE.p);
c    = UnpackC(STATE.c);
E    = UnpackCell(STATE.E);
mu_e = UnpackCell(STATE.mu_e);
chi  = UnpackChi(STATE.chi);
eta  = reshape(PARAM.eta,1,[]);

pars = MODEL.pars;

if isfield(MODEL,'phase_index') && ~isempty(MODEL.phase_index)
    phase_index = MODEL.phase_index(:).';
else
    phase_index = 1:numel(pars);
end

if isfield(STATE,'LE_state') && ~isempty(STATE.LE_state)
    LE_state = STATE.LE_state;
else
    LE_state = struct();
end

% Collapse repeated grains into thermodynamic phases
[pars,c,p,grain_to_phase] = CollapsePhases(pars,c,p,phase_index);

pars_orig = pars;

if use_WScale
    pars_inter = Apply_WScale_FromP(pars,p,wscale_p0,wscale_p1,wscale_gam);
else
    pars_inter = pars;
end

Np = numel(c);
Ne = numel(E);
N  = numel(E{1});

% -------------------------------------------------------------------------
% Active set
% -------------------------------------------------------------------------
[p_le,active,LE_state] = BuildActiveSet(p,LE_state,Pmax,p_tail,p_full,p_on,p_off);
[active_sets,~,set_id] = unique(active,'rows');

mode_map = zeros(1,N);   % 0 = LE, 1 = GP
iter_map = zeros(1,N);

% -------------------------------------------------------------------------
% Loop over active phase combinations
% -------------------------------------------------------------------------
for iset = 1:size(active_sets,1)

    ph_act = find(active_sets(iset,:));
    k      = numel(ph_act);
    mask   = (set_id == iset).';

    if use_WScale && k > 1
        pars_cur = SliceParsWScale(pars_inter,ph_act,mask);
    else
        pars_cur = pars(ph_act);
    end

    % Active unknowns and active mass balance.
    % Important:
    %   p_le selected the active set, but the local solve uses physical p.
    %   Inactive phase contributions are subtracted from E.
    [c_blk,p_blk,E_solve_blk,Efix_blk] = SliceActiveMassBalance( ...
        pars_orig,c,p,E,ph_act,mask);

    % Full conserved E is kept separately for assignment back.
    E_blk   = SliceCell(E,mask);
    mu_blk  = SliceCell(mu_e,mask);
    chi_blk = SliceChi(chi,mask);
    eta_blk = eta(mask);

    nblk = nnz(mask);

    if strcmp(mode,'gp')
        gp_local = true(1,nblk);
    elseif strcmp(mode,'auto')
        gp_local = AutoSafeMask(pars_cur,c_blk,h_tol);
    else
        gp_local = false(1,nblk);
    end

    % Minimal mixed-mode rule:
    % In GP mode, pure one-phase active-set blocks are solved by LE.
    % Multiphase blocks still use GP.
    if gp_pure_LE && strcmp(mode,'gp') && k == 1
        gp_local(:) = false;
    end

    le_local = ~gp_local;

    % ---------------------------------------------------------------------
    % GP nodes: smooth fixed-mu projected closure
    % ---------------------------------------------------------------------
    if any(gp_local)

        kk = min(k,numel(alpha_LE));
        aa = alpha_LE(kk);
        mm = iter_GP(min(k,numel(iter_GP)));

        c_gp   = SliceC(c_blk,1:k,gp_local);
        p_gp   = SliceP(p_blk,gp_local,1:k);
        mu_gp  = SliceCell(mu_blk,gp_local);
        E_gp0  = SliceCell(E_solve_blk,gp_local);
        eta_gp = eta_blk(gp_local);

        pars_gp = SliceParsLocal(pars_cur,gp_local);

        [c_gp,E_gp_active,chi_gp,niter] = GP_Project( ...
            pars_gp,p_gp,c_gp,mu_gp,E_gp0,eta_gp,aa,mm, ...
            lam_gp,wmu_gp,c_tol,gp_dc_cap,gp_c_floor,gp_p_freeze,gp_p_solve);

        c_blk   = AssignC(c_blk,1:k,gp_local,c_gp);
        chi_blk = AssignChi(chi_blk,gp_local,chi_gp);

        if write_GP_E
            Efix_gp = Efix_blk(:,gp_local);
            E_gp_full = AddEfixToE(E_gp_active,Efix_gp);
            E_blk = AssignCell(E_blk,gp_local,E_gp_full);
        end

        loc = find(mask);
        mode_map(loc(gp_local)) = 1;
        iter_map(loc(gp_local)) = niter;
    end

    % ---------------------------------------------------------------------
    % LE nodes: conserved-E local equilibrium
    % ---------------------------------------------------------------------
    if any(le_local)

        kk = min(k,numel(alpha_LE));
        aa = alpha_LE(kk);
        mm = iter_LE(min(k,numel(iter_LE)));

        c_le   = SliceC(c_blk,1:k,le_local);
        p_le2  = SliceP(p_blk,le_local,1:k);
        E_le   = SliceCell(E_solve_blk,le_local);
        eta_le = eta_blk(le_local);

        pars_le = SliceParsLocal(pars_cur,le_local);

        [c_le,mu_le,chi_le,Dle] = LE_Calculator( ...
            pars_le,p_le2,c_le,E_le,eta_le,[aa,mm]);

        c_blk   = AssignC(c_blk,1:k,le_local,c_le);
        mu_blk  = AssignCell(mu_blk,le_local,mu_le);
        chi_blk = AssignChi(chi_blk,le_local,chi_le);

        if isstruct(Dle) && isfield(Dle,'iter')
            loc = find(mask);
            iter_map(loc(le_local)) = Dle.iter;
        end
    end

    % Assign block back
    c    = AssignC(c,ph_act,mask,c_blk);
    E    = AssignCell(E,mask,E_blk);
    mu_e = AssignCell(mu_e,mask,mu_blk);
    chi  = AssignChi(chi,mask,chi_blk);
end

% -------------------------------------------------------------------------
% Pack output
% -------------------------------------------------------------------------
c_col = PackC(c,ny);
E     = PackCell(E,ny);
mu_e  = PackCell(mu_e,ny);
chi   = PackChi(chi,ny);

% Always calculate e from original thermodynamics
e_col = Calc_e(pars_orig,c_col);

% Raw omega from full original thermodynamics
omg_raw_col = CalcOmegaLocal(pars_orig,c,e_col,mu_e,ny,nx,Ne,Np);

% Omega used by AC
if use_WScale && use_WScale_omega
    omg_col = CalcOmegaLocal(pars_inter,c,e_col,mu_e,ny,nx,Ne,Np);
else
    omg_col = omg_raw_col;
end

% Expand thermodynamic phases back to grains
Ngrain = numel(grain_to_phase);
c_out  = cell(1,Ngrain);
omg    = zeros(ny,nx,Ngrain);
omg_raw = zeros(ny,nx,Ngrain);

for ig = 1:Ngrain
    iph            = grain_to_phase(ig);
    c_out{ig}      = c_col{iph};
    omg(:,:,ig)    = omg_col(:,:,iph);
    omg_raw(:,:,ig)= omg_raw_col(:,:,iph);
end

LE_state.phase_index    = phase_index;
LE_state.grain_to_phase = grain_to_phase;
LE_state.active         = active;
LE_state.p_le           = p_le;

STATE.c        = c_out;
STATE.e        = Calc_e(MODEL.pars,c_out);
STATE.E        = E;
STATE.mu_e     = mu_e;
STATE.chi      = chi;
STATE.omg      = omg;
STATE.omg_raw  = omg_raw;
STATE.LE_state = LE_state;

DIAG.mode              = mode;
DIAG.mode_map          = reshape(mode_map,ny,[]);
DIAG.gp_fraction       = mean(mode_map == 1);
DIAG.le_fraction       = mean(mode_map == 0);
DIAG.iter_map          = reshape(iter_map,ny,[]);
DIAG.active            = active;
DIAG.p_le              = p_le;
DIAG.use_WScale        = use_WScale;
DIAG.use_WScale_omega  = use_WScale_omega;
DIAG.gp_pure_LE        = gp_pure_LE;
DIAG.WScale_p0         = wscale_p0;
DIAG.WScale_p1         = wscale_p1;
DIAG.WScale_gam        = wscale_gam;

end


% =========================================================================
% GP projected update
% =========================================================================

function [c,E,chi,niter] = GP_Project( ...
    pars,p,c,mu_e,E_in,eta,alpha,Miter,lam,wmu,c_tol, ...
    dc_cap,c_floor,p_freeze,p_solve)

niter = 0;

for it = 1:Miter

    niter = it;

    dc = GP_Step(pars,p,c,mu_e,E_in,eta,lam,wmu);

    % Stabilize disappearing spinodal phases
    dc = DampStepByPhaseFraction(dc,p,p_freeze,p_solve);
    dc = CapCellStep(dc,dc_cap);

    if MaxStep(dc) < c_tol
        break
    end

    c = AddStepBounded(c,dc,alpha,c_floor);
end

E   = E_FromMu(pars,c,p,mu_e,eta);
chi = Chi_FromMu(pars,p,c,eta,numel(mu_e),lam,wmu);

end


function dc = GP_Step(pars,p,c,mu_e,E_in,eta,lam,wmu)

Np     = numel(c);
Ne     = numel(mu_e);
N      = numel(mu_e{1});
mu_mat = StackFields(mu_e,N);
E_mat  = StackFields(E_in,N);
eta_v  = EtaVector(eta,N);

dc = ZeroStepLike(c);

Emix = zeros(Ne,N);
Rall = cell(1,Np);

for ip = 1:Np
    Rall{ip} = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + StackFields(Rall{ip}.e,N).*p_ip;
end

R_E_global = E_mat - Emix - mu_mat./eta_v;

for ip = 1:Np

    R = Rall{ip};

    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    J  = R.Jac;
    H  = SymPages(R.H_c);
    JT = permute(J,[2 1 3]);
    HT = permute(H,[2 1 3]);
    Nc = size(H,1);

    mu_c = StackFields(R.mu_c,N);
    R_mu = reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),Nc,N) - mu_c;

    p_ip = reshape(p(:,:,ip),1,N);
    R_E  = R_E_global.*p_ip;

    A = pagemtimes(JT,J) + wmu*pagemtimes(HT,H) + ...
        repmat(eye(Nc),1,1,N).*lam;

    b = pagemtimes(JT,reshape(R_E,Ne,1,N)) + ...
        wmu*pagemtimes(HT,reshape(R_mu,Nc,1,N));

    d = reshape(pagemldivide(A,b),Nc,N);

    for ic = 1:Nc
        dc{ip}{ic} = d(ic,:);
    end
end

end


function E = E_FromMu(pars,c,p,mu_e,eta)

Np    = numel(c);
Ne    = numel(mu_e);
N     = numel(mu_e{1});
eta_v = EtaVector(eta,N);

Emix = zeros(Ne,N);

for ip = 1:Np
    R    = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);
    Emix = Emix + StackFields(R.e,N).*p_ip;
end

E = cell(1,Ne);
for ie = 1:Ne
    E{ie} = Emix(ie,:) + reshape(mu_e{ie},1,N)./eta_v;
end

end


function chi = Chi_FromMu(pars,p,c,eta,Ne,lam,wmu)

Np    = numel(c);
N     = numel(c{1}{1});
eta_v = EtaVector(eta,N);

Cmix = zeros(Ne,Ne,N);

for ip = 1:Np

    R    = PhaseThermo(pars{ip},c{ip});
    p_ip = reshape(p(:,:,ip),1,N);

    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)

        Cphase = zeros(Ne,Ne,N);

    else

        J  = R.Jac;
        H  = SymPages(R.H_c);
        JT = permute(J,[2 1 3]);
        HT = permute(H,[2 1 3]);
        Nc = size(H,1);

        A = pagemtimes(JT,J) + wmu*pagemtimes(HT,H) + ...
            repmat(eye(Nc),1,1,N).*lam;

        B = wmu*pagemtimes(HT,JT);

        Cphase = pagemtimes(J,pagemldivide(A,B));
    end

    Cmix = Cmix + Cphase.*reshape(p_ip,1,1,N);
end

I   = repmat(eye(Ne),1,1,N);
chi = PagesToChi(Cmix + I.*reshape(1./eta_v,1,1,N),Ne);

end


% =========================================================================
% Auto mode
% =========================================================================

function safe = AutoSafeMask(pars,c,h_tol)

N = numel(c{1}{1});
safe = true(1,N);

for ip = 1:numel(c)

    R = PhaseThermo(pars{ip},c{ip});

    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    H = SymPages(R.H_c);

    for i = 1:N
        ev = eig(H(:,:,i));
        sc = max(1,max(abs(ev)));
        safe(i) = safe(i) && min(ev) > h_tol*sc;
    end
end

end


% =========================================================================
% Active-set construction
% =========================================================================

function [p_le,active,LE_state] = BuildActiveSet(p,LE_state,Pmax,p_tail,p_full,p_on,p_off)

N  = size(p,2);
Np = size(p,3);

p_th = CalcThermoP(p,p_tail,p_full);
p2   = reshape(p_th,N,Np);

if isfield(LE_state,'active') && isequal(size(LE_state.active),[N,Np])
    active_old = LE_state.active;
else
    active_old = false(N,Np);
end

active = (active_old & p2 > p_off) | (~active_old & p2 > p_on);

for i = 1:N

    if ~any(active(i,:))
        [~,im] = max(p2(i,:));
        active(i,im) = true;
    end

    if sum(active(i,:)) > Pmax
        score = p2(i,:) + 0.5*p_on*active_old(i,:);
        [~,ord] = sort(score,'descend');

        active(i,:) = false;
        active(i,ord(1:Pmax)) = true;
    end
end

p_le = p_th.*reshape(active,1,N,Np);
p_le = p_le./max(sum(p_le,3),eps);

LE_state.active = active;
LE_state.p_th   = p_th;
LE_state.p_le   = p_le;

end


function p_th = CalcThermoP(p,p_tail,p_full)

x = (p - p_tail)./(p_full - p_tail);
x = min(max(x,0),1);

w    = p.*x.^2.*(3 - 2*x);
wsum = sum(w,3);

p_th = w./max(wsum,eps);

bad = wsum <= eps;

if any(bad(:))
    [~,idmax] = max(p,[],3);

    for ip = 1:size(p,3)
        tmp = p_th(:,:,ip);
        tmp(bad & idmax == ip) = 1;
        p_th(:,:,ip) = tmp;
    end
end

end


% =========================================================================
% Active mass-balance slicing
% =========================================================================

function [c_blk,p_blk,E_blk,Efix] = SliceActiveMassBalance(pars,c,p,E,ph_act,mask)
%SLICEACTIVEMASSBALANCE
% The active set is selected using p_le, but the local mass balance uses
% physical p. Inactive phase contributions are subtracted from conserved E.
%
% E_total = sum_active p_i e_i + sum_inactive p_i e_i
% E_active = E_total - sum_inactive p_i e_i

Np   = numel(c);
Ne   = numel(E);
Nloc = nnz(mask);

c_blk = SliceC(c,ph_act,mask);
p_blk = SliceP(p,mask,ph_act);

E_blk = SliceCell(E,mask);

inactive = setdiff(1:Np,ph_act);
Efix = zeros(Ne,Nloc);

for ii = 1:numel(inactive)

    ip = inactive(ii);

    c_tmp = SliceC(c,ip,mask);
    c_ip  = c_tmp{1};

    R = PhaseThermo(pars{ip},c_ip);

    p_ip = reshape(p(:,mask,ip),1,Nloc);
    e_ip = StackFields(R.e,Nloc);

    Efix = Efix + e_ip.*p_ip;
end

for ie = 1:Ne
    E_blk{ie} = E_blk{ie} - Efix(ie,:);
end

end


function E = AddEfixToE(E_active,Efix)

E = E_active;

for ie = 1:numel(E)
    E{ie} = reshape(E{ie},1,[]) + Efix(ie,:);
end

end


% =========================================================================
% Slice / assign helpers
% =========================================================================

function c2 = SliceC(c,ph,mask)

c2 = cell(1,numel(ph));

for a = 1:numel(ph)
    ip = ph(a);
    c2{a} = cell(size(c{ip}));

    for ic = 1:numel(c{ip})
        x = reshape(c{ip}{ic},1,[]);
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
    x = reshape(f{i},1,[]);
    f2{i} = x(mask);
end

end


function ch2 = SliceChi(ch,mask)

Ne  = size(ch,1);
ch2 = cell(Ne,Ne);

for i = 1:Ne
    for j = 1:Ne
        x = reshape(ch{i,j},1,[]);
        ch2{i,j} = x(mask);
    end
end

end


function c = AssignC(c,ph,mask,c_sub)

for a = 1:numel(ph)
    ip = ph(a);

    for ic = 1:numel(c{ip})
        x = reshape(c{ip}{ic},1,[]);
        x(mask) = reshape(c_sub{a}{ic},1,[]);
        c{ip}{ic} = x;
    end
end

end


function f = AssignCell(f,mask,f_sub)

for i = 1:numel(f)
    x = reshape(f{i},1,[]);
    x(mask) = reshape(f_sub{i},1,[]);
    f{i} = x;
end

end


function ch = AssignChi(ch,mask,ch_sub)

Ne = size(ch,1);

for i = 1:Ne
    for j = 1:Ne
        x = reshape(ch{i,j},1,[]);
        x(mask) = reshape(ch_sub{i,j},1,[]);
        ch{i,j} = x;
    end
end

end


% =========================================================================
% Pack / unpack
% =========================================================================

function p = UnpackP(p)

[ny,nx,np] = size(p);
p = reshape(p,1,ny*nx,np);

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
% Collapse repeated grains with same thermodynamic phase
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

    grains = find(grain_to_phase == iph);
    ig0    = grains(1);

    pars_c{iph} = pars{ig0};

    for ig = grains
        p_c(:,:,iph) = p_c(:,:,iph) + p(:,:,ig);
    end

    Nc   = numel(c{ig0});
    den  = reshape(p_c(:,:,iph),1,N);
    good = den > eps;

    c_c{iph} = cell(1,Nc);

    for ic = 1:Nc
        num = zeros(1,N);

        for ig = grains
            num = num + reshape(p(:,:,ig),1,N).*reshape(c{ig}{ic},1,N);
        end

        tmp = reshape(c{ig0}{ic},1,N);
        tmp(good) = num(good)./den(good);

        c_c{iph}{ic} = tmp;
    end
end

end


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
% GP stabilization helpers
% =========================================================================

function dc = DampStepByPhaseFraction(dc,p,p_freeze,p_solve)

for ip = 1:numel(dc)

    p_ip = reshape(p(:,:,ip),1,[]);

    w = (p_ip - p_freeze)./(p_solve - p_freeze);
    w = min(max(w,0),1);
    w = w.^2.*(3 - 2*w);

    for ic = 1:numel(dc{ip})
        dc{ip}{ic} = dc{ip}{ic}.*w;
    end
end

end


function dc = CapCellStep(dc,dc_cap)

if dc_cap <= 0 || isinf(dc_cap)
    return
end

amp = zeros(size(dc{1}{1}));

for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        amp = max(amp,abs(dc{ip}{ic}));
    end
end

fac = ones(size(amp));
bad = amp > dc_cap;
fac(bad) = dc_cap./amp(bad);

for ip = 1:numel(dc)
    for ic = 1:numel(dc{ip})
        dc{ip}{ic} = dc{ip}{ic}.*fac;
    end
end

end


function c_new = AddStepBounded(c,dc,alpha,c_floor)

c_new = c;

for ip = 1:numel(c)

    Nc = numel(c{ip});
    N  = numel(c{ip}{1});
    X  = zeros(Nc,N);

    for ic = 1:Nc
        X(ic,:) = c{ip}{ic} + alpha.*dc{ip}{ic};
    end

    X(~isfinite(X)) = c_floor;
    X = max(X,c_floor);

    s = sum(X,1);
    bad = s > 1 - c_floor;

    if any(bad)
        X(:,bad) = X(:,bad)./s(bad).*(1 - c_floor);
    end

    for ic = 1:Nc
        c_new{ip}{ic} = reshape(X(ic,:),size(c{ip}{ic}));
    end
end

end


% =========================================================================
% Small helpers
% =========================================================================

function H = SymPages(H)

H = 0.5*(H + permute(H,[2 1 3]));

end


function A = StackFields(fields,N)

A = zeros(numel(fields),N);

for i = 1:numel(fields)
    A(i,:) = reshape(fields{i},1,N);
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

        else
            error('SliceParsWScale: w_scale size mismatch. numel(w_scale)=%d, numel(mask)=%d, nnz(mask)=%d.', ...
                numel(ws),numel(mask),nnz(mask));
        end
    end
end

end


function pars_out = SliceParsLocal(pars_in,local_mask)
%SLICEPARSLOCAL Further slice nodewise w_scale inside one active-set block.

pars_out = pars_in;

for ip = 1:numel(pars_out)

    if isfield(pars_out{ip},'w_scale') && ~isempty(pars_out{ip}.w_scale)

        ws = reshape(pars_out{ip}.w_scale,1,[]);

        if isscalar(ws)
            pars_out{ip}.w_scale = ws;

        elseif numel(ws) == numel(local_mask)
            pars_out{ip}.w_scale = ws(local_mask);

        elseif numel(ws) == nnz(local_mask)
            pars_out{ip}.w_scale = ws;

        else
            error('SliceParsLocal: w_scale size mismatch. numel(w_scale)=%d, numel(mask)=%d, nnz(mask)=%d.', ...
                numel(ws),numel(local_mask),nnz(local_mask));
        end
    end
end

end