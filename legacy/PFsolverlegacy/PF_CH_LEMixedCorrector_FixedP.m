function [STATE_OUT,DIAG] = PF_CH_LEMixedCorrector_FixedP(STATE_OLD,STATE_BASE,PARAM,MODEL,GRID,PHYS,NUM)
%PF_CH_LEMIXEDCORRECTOR_FIXEDP  Fixed-p mixed CH-LE corrector.
%
% First implementation of the selectively-uncondensed spinodal idea.
%
% Use after your AC+CH predictor and a thermodynamic response:
%
%   [STATE_RAW,DIAG] = PF_Coupled_ACCH_LETangent(...);
%   STATE_LE0        = LE_Run_FromMu(STATE_RAW,PARAM,MODEL);
%   [STATE_TRIAL,DIAG_MIX] = PF_CH_LEMixedCorrector_FixedP( ...
%        STATE_OLD,STATE_LE0,PARAM,MODEL,GRID,PHYS,NUM);
%
% Stable phase/node:
%   dE = chi_safe*dmu
%
% Spinodal or near-singular phase/node:
%   do not condense into chi.  Add dc as unknown and solve
%       Hc*dc - J'*dmu = -r
%   and enter CH through
%       dE += p*J*dc
%
% This function keeps phi and p fixed.

[ny,nx,Np] = size(STATE_BASE.p);
Ne         = numel(STATE_BASE.E);
N          = nx*ny;
dt         = NUM.dt_phy;

dx = GRID.dx; dy = GRID.dy;
dx2 = dx^2; dy2 = dy^2;
dx4 = dx2^2; dy4 = dy2^2;

if isfield(PARAM,'kappa_eff') && ~isempty(PARAM.kappa_eff)
    kappa = PARAM.kappa_eff;
    if isscalar(kappa), kappa = kappa*ones(ny,nx); end
else
    kappa = PHYS.kappa*ones(ny,nx);
end

if isscalar(PARAM.eta)
    eta = PARAM.eta*ones(ny,nx);
else
    eta = PARAM.eta;
end
eta_vec = eta(:).';

spin_on   = getfield_default(NUM,'MIX_spinodal_on',0);
rcond_min = getfield_default(NUM,'MIX_rcond_min',1e-10);
dc_limit  = getfield_default(NUM,'MIX_dc_limit',inf);
mu_limit  = getfield_default(NUM,'MIX_mu_limit',inf);

% -------------------------------------------------------------------------
% Thermodynamics and unsafe node detection
% -------------------------------------------------------------------------
Rph = cell(1,Np);
Nc_phase = zeros(1,Np);
unsafe = false(Np,N);
lam_min = inf(Np,N);
rc_H    = ones(Np,N);

for ip = 1:Np
    cvec = pack_c_phase(STATE_BASE.c{ip});
    Rph{ip} = PhaseThermo(MODEL.pars{ip},cvec);

    if isempty(Rph{ip}.H_c) || isempty(Rph{ip}.Jac)
        continue
    end

    Nc_phase(ip) = size(Rph{ip}.H_c,1);
    p_ip = reshape(STATE_BASE.p(:,:,ip),1,N);

    for q = 1:N
        H = 0.5*(Rph{ip}.H_c(:,:,q)+Rph{ip}.H_c(:,:,q).');

        if any(~isfinite(H(:)))
            unsafe(ip,q) = true;
            lam_min(ip,q) = -inf;
            rc_H(ip,q) = 0;
        else
            ev = eig(H);
            lam_min(ip,q) = min(ev);
            rc_H(ip,q) = rcond(H);
            unsafe(ip,q) = p_ip(q) > 1e-12 && ...
                (lam_min(ip,q) <= spin_on || rc_H(ip,q) < rcond_min);
        end
    end
end

% -------------------------------------------------------------------------
% Unknown IDs
% -------------------------------------------------------------------------
idMu = cell(1,Ne);
for ie = 1:Ne
    idMu{ie} = ((ie-1)*N + (1:N)).';
end
Nmu = Ne*N;

[idC,Nc_unknown] = allocate_c_ids(unsafe,Nc_phase,N,Nmu);
Ntot = Nmu + Nc_unknown;

% -------------------------------------------------------------------------
% Safe condensed chi and unsafe pJ
% -------------------------------------------------------------------------
chi_safe = repmat(eye(Ne),1,1,N).*reshape(1./eta_vec,1,1,N);

for ip = 1:Np
    if isempty(Rph{ip}.chi), continue; end
    S = 0.5*(Rph{ip}.chi + permute(Rph{ip}.chi,[2 1 3]));
    p_ip = reshape(STATE_BASE.p(:,:,ip),1,1,N);
    safe = reshape(~unsafe(ip,:),1,1,N);
    chi_safe = chi_safe + S.*p_ip.*safe;
end

% -------------------------------------------------------------------------
% Grid indexing, matching current reflected convention
% -------------------------------------------------------------------------
[I,J] = ndgrid(1:ny,1:nx);
ii = I(:); jj = J(:);
refI = @(a,s) reflect_index_mixed(a+s,ny);
refJ = @(a,s) reflect_index_mixed(a+s,nx);

idx_c = sub2ind([ny,nx],ii,jj);
idx_L = sub2ind([ny,nx],ii,refJ(jj,-1));
idx_R = sub2ind([ny,nx],ii,refJ(jj,+1));
idx_U = sub2ind([ny,nx],refI(ii,-1),jj);
idx_D = sub2ind([ny,nx],refI(ii,+1),jj);
idx_L2 = sub2ind([ny,nx],ii,refJ(jj,-2));
idx_R2 = sub2ind([ny,nx],ii,refJ(jj,+2));
idx_U2 = sub2ind([ny,nx],refI(ii,-2),jj);
idx_D2 = sub2ind([ny,nx],refI(ii,+2),jj);
idx_UR = sub2ind([ny,nx],refI(ii,-1),refJ(jj,+1));
idx_DR = sub2ind([ny,nx],refI(ii,+1),refJ(jj,+1));
idx_UL = sub2ind([ny,nx],refI(ii,-1),refJ(jj,-1));
idx_DL = sub2ind([ny,nx],refI(ii,+1),refJ(jj,-1));

% -------------------------------------------------------------------------
% Sparse assembly
% -------------------------------------------------------------------------
maxnnz = Ne*N*(5 + 13*Ne + 13*max(1,sum(Nc_phase))) + ...
         max(1,Nc_unknown)*(Ne+max(1,max(Nc_phase))+1) + 1000;
rows = zeros(maxnnz,1); cols = rows; vals = rows;
b = zeros(Ntot,1);
kk = 1;

for l = 1:Ne

    row = idMu{l}(idx_c);

    M = PARAM.M{l};
    Mc = M(idx_c);
    ML = M(idx_L); MR = M(idx_R);
    MU = M(idx_U); MD = M(idx_D);

    dL = -(ML+Mc)/2/dx2;
    dR = -(MR+Mc)/2/dx2;
    dU = -(MU+Mc)/2/dy2;
    dD = -(MD+Mc)/2/dy2;
    dC = -(dL+dR+dU+dD);

    MK = Mc.*kappa(idx_c);

    aL  = MK.*(-4/dx4 - 4/(dx2*dy2));
    aR  = aL;
    aU  = MK.*(-4/dy4 - 4/(dx2*dy2));
    aD  = aU;
    aL2 = MK.*(1/dx4); aR2 = aL2;
    aU2 = MK.*(1/dy4); aD2 = aU2;
    aUR = MK.*(2/(dx2*dy2)); aDR = aUR; aUL = aUR; aDL = aUR;
    aC  = 1/dt + MK.*(6/dx4 + 6/dy4 + 8/(dx2*dy2));

    E = STATE_BASE.E{l};
    mu = STATE_BASE.mu_e{l};

    AE_E = aC.*E(idx_c) + aL.*E(idx_L) + aR.*E(idx_R) + ...
           aU.*E(idx_U) + aD.*E(idx_D) + ...
           aL2.*E(idx_L2) + aR2.*E(idx_R2) + ...
           aU2.*E(idx_U2) + aD2.*E(idx_D2) + ...
           aUR.*E(idx_UR) + aDR.*E(idx_DR) + ...
           aUL.*E(idx_UL) + aDL.*E(idx_DL);

    D_mu = dC.*mu(idx_c) + dL.*mu(idx_L) + dR.*mu(idx_R) + ...
           dU.*mu(idx_U) + dD.*mu(idx_D);

    b(row) = STATE_OLD.E{l}(idx_c)/dt - AE_E - D_mu;

    % D*dmu_l
    [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{l}(idx_c),dC);
    [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{l}(idx_L),dL);
    [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{l}(idx_R),dR);
    [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{l}(idx_U),dU);
    [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{l}(idx_D),dD);

    % AE*(chi_safe*dmu)
    for m = 1:Ne
        Ch = reshape(chi_safe(l,m,:),ny,nx);

        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_c), aC.*Ch(idx_c));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_L), aL.*Ch(idx_L));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_R), aR.*Ch(idx_R));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_U), aU.*Ch(idx_U));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_D), aD.*Ch(idx_D));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_L2), aL2.*Ch(idx_L2));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_R2), aR2.*Ch(idx_R2));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_U2), aU2.*Ch(idx_U2));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_D2), aD2.*Ch(idx_D2));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_UR), aUR.*Ch(idx_UR));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_DR), aDR.*Ch(idx_DR));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_UL), aUL.*Ch(idx_UL));
        [rows,cols,vals,kk] = add_block_mixed(rows,cols,vals,kk,row,idMu{m}(idx_DL), aDL.*Ch(idx_DL));
    end

    % AE*(p*J*dc) for unsafe phases
    for ip = 1:Np
        if ~any(unsafe(ip,:)), continue; end
        Jip = Rph{ip}.Jac;
        Nc = Nc_phase(ip);
        p_ip = STATE_BASE.p(:,:,ip);

        for ic = 1:Nc
            idmap = idC{ip,ic};
            if isempty(idmap) || ~any(idmap(:)), continue; end

            P = reshape(p_ip(:).'.*reshape(Jip(l,ic,:),1,N).*unsafe(ip,:),ny,nx);

            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_c, aC.*P(idx_c));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_L, aL.*P(idx_L));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_R, aR.*P(idx_R));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_U, aU.*P(idx_U));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_D, aD.*P(idx_D));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_L2, aL2.*P(idx_L2));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_R2, aR2.*P(idx_R2));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_U2, aU2.*P(idx_U2));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_D2, aD2.*P(idx_D2));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_UR, aUR.*P(idx_UR));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_DR, aDR.*P(idx_DR));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_UL, aUL.*P(idx_UL));
            [rows,cols,vals,kk] = add_active_mixed(rows,cols,vals,kk,row,idmap,idx_DL, aDL.*P(idx_DL));
        end
    end
end

% Local thermodynamic rows Hc*dc - J'*dmu = -r
mu_base_mat = pack_E_cell(STATE_BASE.mu_e,N);

for ip = 1:Np
    if ~any(unsafe(ip,:)), continue; end

    Ri = Rph{ip};
    Nc = Nc_phase(ip);
    mu_c = cell2mat(Ri.mu_c(:));
    Jip = Ri.Jac;
    Hc  = Ri.H_c;

    JT = permute(Jip,[2 1 3]);
    JTmu3 = pagemtimes(JT,reshape(mu_base_mat,Ne,1,N));
    JTmu = reshape(JTmu3,Nc,N);
    res = mu_c - JTmu;

    nodes = find(unsafe(ip,:));

    for q = nodes
        for a = 1:Nc
            row = idC{ip,a}(q);
            b(row) = -res(a,q);

            for bb = 1:Nc
                col = idC{ip,bb}(q);
                [rows,cols,vals,kk] = add_scalar_mixed(rows,cols,vals,kk,row,col,Hc(a,bb,q));
            end

            for m = 1:Ne
                col = idMu{m}(q);
                [rows,cols,vals,kk] = add_scalar_mixed(rows,cols,vals,kk,row,col,-Jip(m,a,q));
            end
        end
    end
end

rows = rows(1:kk-1); cols = cols(1:kk-1); vals = vals(1:kk-1);
A = sparse(rows,cols,vals,Ntot,Ntot);

qord = colamd(A);
y = A(:,qord)\b;
sol = zeros(Ntot,1);
sol(qord) = y;

relres = norm(A*sol-b)/max(norm(b),eps);

% Unpack dmu
dmu = cell(1,Ne);
for ie = 1:Ne
    dmu{ie} = reshape(sol(idMu{ie}),ny,nx);
end

% Unpack dc unsafe
dc = make_zero_c(STATE_BASE.c);
for ip = 1:Np
    Nc = Nc_phase(ip);
    for ic = 1:Nc
        idmap = idC{ip,ic};
        if isempty(idmap) || ~any(idmap(:)), continue; end
        tmp = zeros(ny,nx);
        mask = idmap > 0;
        tmp(mask) = sol(idmap(mask));
        dc{ip}{ic} = tmp;
    end
end

% Build output
STATE_OUT = STATE_BASE;

for ie = 1:Ne
    STATE_OUT.mu_e{ie} = STATE_BASE.mu_e{ie} + dmu{ie};
end

% First apply unsafe dc
for ip = 1:Np
    for ic = 1:numel(STATE_OUT.c{ip})
        STATE_OUT.c{ip}{ic} = STATE_BASE.c{ip}{ic} + dc{ip}{ic};
    end
end

% Then apply stable condensed dc
STATE_OUT.c = apply_stable_dc(STATE_OUT.c,STATE_BASE.c,STATE_BASE.mu_e,Rph,unsafe,dmu,Nc_phase,Ne,N,ny,nx);

STATE_OUT.e = Calc_e(MODEL.pars,STATE_OUT.c);
STATE_OUT.E = calc_E_from_c_mu(STATE_OUT.e,STATE_OUT.p,STATE_OUT.mu_e,eta);
STATE_OUT.omg = calc_omega_local(MODEL.pars,STATE_OUT.c,STATE_OUT.e,STATE_OUT.mu_e);

% Return safe chi for next condensed predictor. Unsafe contribution excluded.
[Rnew,unsafe_new,lam_new,rc_new,~] = local_reclassify(MODEL.pars,STATE_OUT.c,STATE_OUT.p,spin_on,rcond_min);
chi_safe_new = build_safe_chi_from_R(Rnew,STATE_OUT.p,unsafe_new,eta_vec,Ne,N);
STATE_OUT.chi = unpack_chi(chi_safe_new,ny,nx);

% Diagnostics
max_dmu = 0;
for ie = 1:Ne
    max_dmu = max(max_dmu,max(abs(dmu{ie}(:))));
end
max_dc = 0;
for ip = 1:Np
    for ic = 1:numel(dc{ip})
        max_dc = max(max_dc,max(abs(dc{ip}{ic}(:))));
    end
end

DIAG.relres = relres;
DIAG.matrix_size = Ntot;
DIAG.nnz = nnz(A);
DIAG.Nmu = Nmu;
DIAG.Nc_unsafe = Nc_unknown;
DIAG.max_dmu = max_dmu;
DIAG.max_dc = max_dc;
DIAG.unsafe = reshape(any(unsafe,1),ny,nx);
DIAG.unsafe_frac = nnz(any(unsafe,1))/N;
DIAG.lam_min = reshape(min(lam_min,[],1),ny,nx);
DIAG.rcond_min = reshape(min(rc_H,[],1),ny,nx);
DIAG.lam_min_new = reshape(min(lam_new,[],1),ny,nx);
DIAG.rcond_new = reshape(min(rc_new,[],1),ny,nx);
DIAG.rejected = false;
DIAG.reject_reason = '';

if ~isfinite(relres) || ~isfinite(max_dmu) || ~isfinite(max_dc)
    DIAG.rejected = true;
    DIAG.reject_reason = 'nonfinite mixed solve';
elseif max_dmu > mu_limit
    DIAG.rejected = true;
    DIAG.reject_reason = 'dmu limit exceeded';
elseif max_dc > dc_limit
    DIAG.rejected = true;
    DIAG.reject_reason = 'dc limit exceeded';
end

end

% ============================== helpers ==============================
function cvec = pack_c_phase(cphase)
cvec = cell(size(cphase));
for i = 1:numel(cphase)
    cvec{i} = reshape(cphase{i},1,[]);
end
end

function E = pack_E_cell(C,N)
E = zeros(numel(C),N);
for i = 1:numel(C)
    E(i,:) = reshape(C{i},1,N);
end
end

function [idC,Nc_unknown] = allocate_c_ids(unsafe,Nc_phase,N,Noffset)
Np = size(unsafe,1);
maxNc = max(max(Nc_phase),1);
idC = cell(Np,maxNc);
counter = Noffset;
for ip = 1:Np
    nodes = find(unsafe(ip,:));
    for ic = 1:Nc_phase(ip)
        idmap = zeros(N,1);
        if ~isempty(nodes)
            ids = counter + (1:numel(nodes)).';
            idmap(nodes) = ids;
            counter = counter + numel(nodes);
        end
        idC{ip,ic} = idmap;
    end
end
Nc_unknown = counter - Noffset;
end

function [R,unsafe,lam,rc,Nc_phase] = local_reclassify(pars,c,p,spin_on,rcond_min)
Np = numel(c); N = numel(c{1}{1});
R = cell(1,Np); unsafe = false(Np,N); lam = inf(Np,N); rc = ones(Np,N); Nc_phase = zeros(1,Np);
for ip = 1:Np
    R{ip} = PhaseThermo(pars{ip},pack_c_phase(c{ip}));
    if isempty(R{ip}.H_c) || isempty(R{ip}.Jac), continue; end
    Nc_phase(ip) = size(R{ip}.H_c,1);
    p_ip = reshape(p(:,:,ip),1,N);
    for q = 1:N
        H = 0.5*(R{ip}.H_c(:,:,q)+R{ip}.H_c(:,:,q).');
        if any(~isfinite(H(:)))
            unsafe(ip,q)=true; lam(ip,q)=-inf; rc(ip,q)=0;
        else
            ev=eig(H); lam(ip,q)=min(ev); rc(ip,q)=rcond(H);
            unsafe(ip,q)=p_ip(q)>1e-12 && (lam(ip,q)<=spin_on || rc(ip,q)<rcond_min);
        end
    end
end
end

function chi_safe = build_safe_chi_from_R(R,p,unsafe,eta_vec,Ne,N)
chi_safe = repmat(eye(Ne),1,1,N).*reshape(1./eta_vec,1,1,N);
for ip = 1:numel(R)
    if isempty(R{ip}.chi), continue; end
    S = 0.5*(R{ip}.chi+permute(R{ip}.chi,[2 1 3]));
    chi_safe = chi_safe + S.*reshape(p(:,:,ip),1,1,N).*reshape(~unsafe(ip,:),1,1,N);
end
end

function C = make_zero_c(c)
C = c;
for ip = 1:numel(c)
    for ic = 1:numel(c{ip})
        C{ip}{ic}=zeros(size(c{ip}{ic}));
    end
end
end

function cnew = apply_stable_dc(cnew,cbase,mu_base,R,unsafe,dmu,Nc_phase,Ne,N,ny,nx)
dmu_mat = pack_E_cell(dmu,N);
mu_mat = pack_E_cell(mu_base,N);
for ip = 1:numel(cbase)
    Nc = Nc_phase(ip);
    if Nc==0, continue; end
    Ri=R{ip}; J=Ri.Jac; Hc=Ri.H_c; mu_c=cell2mat(Ri.mu_c(:));
    JT=permute(J,[2 1 3]);
    JTdm=reshape(pagemtimes(JT,reshape(dmu_mat,Ne,1,N)),Nc,N);
    JTmu=reshape(pagemtimes(JT,reshape(mu_mat,Ne,1,N)),Nc,N);
    res=mu_c-JTmu;
    dc=zeros(Nc,N);
    for q=1:N
        if unsafe(ip,q), continue; end
        H=0.5*(Hc(:,:,q)+Hc(:,:,q).');
        if rcond(H)<1e-12 || any(~isfinite(H(:))), continue; end
        dc(:,q)=H\(JTdm(:,q)-res(:,q));
    end
    for ic=1:Nc
        cnew{ip}{ic}=cnew{ip}{ic}+reshape(dc(ic,:),ny,nx);
    end
end
end

function E = calc_E_from_c_mu(e,p,mu,eta)
Ne=numel(mu); [~,~,Np]=size(p); E=cell(1,Ne);
for ie=1:Ne
    tmp=zeros(size(mu{ie}));
    for ip=1:Np
        tmp=tmp+p(:,:,ip).*e{ip}{ie};
    end
    E{ie}=tmp+mu{ie}./eta;
end
end

function omg = calc_omega_local(pars,c,e,mu)
Np=numel(c); Ne=numel(mu); [ny,nx]=size(mu{1}); omg=zeros(ny,nx,Np);
for ip=1:Np
    omg(:,:,ip)=PhaseG(pars{ip},c{ip});
    for ie=1:Ne
        omg(:,:,ip)=omg(:,:,ip)-e{ip}{ie}.*mu{ie};
    end
end
end

function chi=unpack_chi(pg,ny,nx)
[Ne,~,~]=size(pg); chi=cell(Ne,Ne);
for i=1:Ne
    for j=1:Ne
        chi{i,j}=reshape(pg(i,j,:),ny,nx);
    end
end
end

function [rows,cols,vals,k]=add_block_mixed(rows,cols,vals,k,row,col,val)
n=numel(row); need=k+n-1;
if need>numel(rows)
    grow=max(numel(rows),n+1000);
    rows(end+grow)=0; cols(end+grow)=0; vals(end+grow)=0;
end
rows(k:need)=row(:); cols(k:need)=col(:); vals(k:need)=val(:); k=need+1;
end

function [rows,cols,vals,k]=add_scalar_mixed(rows,cols,vals,k,row,col,val)
if row==0 || col==0, return; end
if k>numel(rows)
    grow=max(numel(rows),1000); rows(end+grow)=0; cols(end+grow)=0; vals(end+grow)=0;
end
rows(k)=row; cols(k)=col; vals(k)=val; k=k+1;
end

function [rows,cols,vals,k]=add_active_mixed(rows,cols,vals,k,row,idmap,idx,val)
col=idmap(idx); mask=col>0;
if any(mask)
    [rows,cols,vals,k]=add_block_mixed(rows,cols,vals,k,row(mask),col(mask),val(mask));
end
end

function idx=reflect_index_mixed(idx,n)
if n==1, idx(:)=1; return; end
period=2*n-2; r=mod(idx-1,period); idx=1+min(r,period-r);
end

function v=getfield_default(S,name,default)
if isfield(S,name) && ~isempty(S.(name)), v=S.(name); else, v=default; end
end
