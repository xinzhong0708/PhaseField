function [Result] = PhaseThermo(pars,c)
% Thermodynamic calculator with optional penalty endmembers.
%
% Full endmember order is ALWAYS:
%   [ real_1 ... real_nReal , pen_1 ... pen_nPen ]
%
% Dependent endmember is ALWAYS the LAST full endmember:
%   - if nPen = 0, last real is dependent
%   - if nPen > 0, last penalty is dependent
%
% Input c can be either:
%   (1) full composition, size = nAll x N
%   (2) independent composition, size = (nAll-1) x N
%
% Output mu_c is always with respect to the last full endmember.

% -------------------------------------------------------------------------
% unpack input
% -------------------------------------------------------------------------
if iscell(c)
    c = cell2mat(c');
end

% -------------------------------------------------------------------------
% load parameters
% -------------------------------------------------------------------------
P      = pars.P/1e8;
T      = pars.T;
RT     = 8.3144*T;

mtpl   = pars.mtpl;
zt     = pars.zt;
alp    = pars.alp;
w      = pars.w;

g0     = pars.g0(:);          % nReal x 1
n_real = pars.n;              % nReal x NeCat

if isfield(pars,'gN') && ~isempty(pars.gN)
    g_pen = pars.gN(:);       % nPen x 1
else
    g_pen = zeros(0,1);
end

if isfield(pars,'nN') && ~isempty(pars.nN)
    n_pen = pars.nN;          % nPen x NeCat
else
    n_pen = zeros(0,size(n_real,2));
end

if isfield(pars,'penalty') && ~isempty(pars.penalty)
    penalty = pars.penalty;
else
    penalty = 1e9;
end

nReal = numel(g0);
nPen  = numel(g_pen);
nAll  = nReal + nPen;

if size(n_pen,1) ~= nPen
    error('PhaseThermo: pars.gN and pars.nN size mismatch.');
end
if size(n_pen,2) ~= size(n_real,2)
    error('PhaseThermo: pars.nN must have same number of columns as pars.n.');
end

% dependent is always the LAST full endmember
dep_idx = nAll;

% -------------------------------------------------------------------------
% construct full c
% -------------------------------------------------------------------------
Npt = size(c,2);

if size(c,1) == nAll
    % full composition supplied directly
    c_full = c;
elseif size(c,1) == nAll-1
    % independent composition supplied, append dependent
    c_full = [c; 1-sum(c,1)];
else
    error('PhaseThermo: input c must have %d or %d rows, got %d.', ...
          nAll-1, nAll, size(c,1));
end

% split according to REAL/PENALTY order
c_real = c_full(1:nReal,:);           % nReal x N
c_pen  = c_full(nReal+1:end,:);       % nPen  x N

c_realt = c_real.';
c_pent  = c_pen.';

% -------------------------------------------------------------------------
% total endmember composition matrix
% -------------------------------------------------------------------------
n_all    = [n_real; n_pen];           % nAll x NeCat
aCat_all = sum(n_all,2);              % nAll x 1

% -------------------------------------------------------------------------
% free energy: REAL contribution
% -------------------------------------------------------------------------
g_mech_real = c_realt * g0;           % N x 1

z     = c_realt * zt + 1e-30;

eps0  = 1e-6;
sz    = sqrt(z.^2 + eps0.^2);
lz    = log(z + sz) - log(2);
dphi_z  = lz + z ./ sz;
d2phi_z = (z.^2 + 2*eps0.^2) ./ (sz.^3);

Azt   = zt .* log(zt + double(zt==0));

g_id  = RT * sum(mtpl .* ( z .* lz - c_realt * Azt ), 2);

alp_eff = [1 T P] * alp;                         % 1 x nReal
W       = w(:,:,1) + w(:,:,2)*T + w(:,:,3)*P;   % nReal x nReal

a       = alp_eff(:);                            % nReal x 1
M       = ((a*a.') ./ (a + a.')) .* W;           % nReal x nReal
q       = c_realt * a;                           % N x 1
invq    = 1 ./ q;
invq2   = invq.^2;
invq3   = invq.^3;

Mc      = c_realt * M;                           % N x nReal
n_nid   = sum(Mc .* c_realt, 2);                 % N x 1
B       = M + M.';                               % nReal x nReal
v       = c_realt * B.';                         % N x nReal
g_nid   = n_nid .* invq;                         % N x 1

% -------------------------------------------------------------------------
% penalty contribution
%   g_pen*c_pen + penalty*sum(c_pen^2)
% -------------------------------------------------------------------------
if nPen > 0
    g_pen_extra = c_pent * g_pen + penalty * sum(c_pent.^2, 2);   % N x 1
else
    g_pen_extra = zeros(Npt,1);
end

% -------------------------------------------------------------------------
% total raw free energy before cation normalization
% -------------------------------------------------------------------------
G_raw = g_mech_real + g_id + g_nid + g_pen_extra;   % N x 1

% -------------------------------------------------------------------------
% normalize to one cation
% -------------------------------------------------------------------------
cation  = c_full.' * aCat_all;      % N x 1
invcat  = 1 ./ cation;
invcat2 = invcat.^2;
invcat3 = invcat.^3;

scale    = 1 / (pars.E_sc * pars.vref);
G_scaled = G_raw * scale;           % N x 1
g        = G_scaled .* invcat;      % N x 1

% -------------------------------------------------------------------------
% first derivative wrt FULL endmember coordinates
% -------------------------------------------------------------------------

% REAL block dG/dc_real
mu_nid  = v .* invq - (n_nid .* invq2) .* a.';      % N x nReal
mu_mech = g0.';                                      % 1 x nReal
mu_id   = RT * ( (mtpl .* dphi_z) * zt.' - mtpl * Azt.' );

dG_real = mu_mech + mu_id + mu_nid;                  % N x nReal

% PENALTY block dG/dc_pen
if nPen > 0
    dG_pen = repmat(g_pen.', Npt, 1) + 2 * penalty * c_pent;   % N x nPen
else
    dG_pen = zeros(Npt,0);
end

% combine and scale
dG_full = [dG_real, dG_pen] * scale;                % N x nAll

% normalized chemical potentials wrt FULL coordinates
mu_full_all = dG_full .* invcat - (G_scaled .* invcat2) * aCat_all.';   % N x nAll

% -------------------------------------------------------------------------
% convert to mu_c wrt independent variables
% dependent is last full endmember
% -------------------------------------------------------------------------
mu_c = (mu_full_all(:,1:end-1) - mu_full_all(:,end)).';   % (nAll-1) x N

% -------------------------------------------------------------------------
% elemental concentration and mu_e
% -------------------------------------------------------------------------
e_all   = (c_full.' * n_all) .* invcat;      % N x NeCat

NeCat    = size(n_all,2);
nElemInd = NeCat - 1;
nVar     = nAll - 1;

% -------------------------------------------------------------------------
% Jacobian de/dc_full
% -------------------------------------------------------------------------
nT        = n_all.';
Jfull_all = repmat(nT,1,1,Npt) ...
          - reshape(e_all.',NeCat,1,Npt) .* reshape(aCat_all.',1,nAll,1);
Jfull_all = Jfull_all .* reshape(invcat.',1,1,Npt);

% map independent variables to full coordinates:
% c_full = [c_ind ; 1-sum(c_ind)]
R = [eye(nVar); -ones(1,nVar)];

R3    = repmat(R,1,1,Npt);
RT3   = repmat(R.',1,1,Npt);

% de/dc_ind
Je     = pagemtimes(Jfull_all, R3);      % NeCat x nVar x N
Je_red = Je(1:end-1,:,:);                % (NeCat-1) x nVar x N

% solve mu_c = Je_red^T * mu_e
A      = permute(Je_red,[2 1 3]);        % nVar x (NeCat-1) x N
b      = reshape(mu_c, nVar, 1, Npt);
mu3    = pagemldivide(A, b);
mu_e   = reshape(mu3, nElemInd, Npt);

% -------------------------------------------------------------------------
% Hessian wrt FULL coordinates before dependent elimination
% -------------------------------------------------------------------------

% REAL-REAL Hessian
H_id = zeros(nReal,nReal,Npt);
for is = 1:size(zt,2)
    zs    = zt(:,is) * zt(:,is).';
    coeff = RT * ( mtpl(:,is) .* d2phi_z(:,is) );
    H_id  = H_id + zs .* reshape(coeff,1,1,Npt);
end

vT      = v.';
aaT     = a * a.';

term2   = reshape(a,nReal,1,1) .* reshape(vT,1,nReal,Npt) ...
        + reshape(vT,nReal,1,Npt) .* reshape(a.',1,nReal,1);

H_nid   = B .* reshape(invq,1,1,Npt) ...
        - term2 .* reshape(invq2,1,1,Npt) ...
        + aaT .* reshape(2*n_nid .* invq3,1,1,Npt);

H_real_full = (H_id + H_nid) * scale;    % nReal x nReal x N

% PEN-PEN Hessian
H_pen_full = zeros(nPen,nPen,Npt);
if nPen > 0
    Ipen       = eye(nPen);
    H_pen_full = repmat(2 * penalty * scale * Ipen, 1, 1, Npt);
end

% FULL Hessian wrt full coordinates
H_G_full = zeros(nAll,nAll,Npt);
H_G_full(1:nReal,1:nReal,:) = H_real_full;
if nPen > 0
    H_G_full(nReal+1:end,nReal+1:end,:) = H_pen_full;
end

% normalized Hessian wrt FULL coordinates
dG_full_T = dG_full.';   % nAll x N

crossTerm = reshape(dG_full_T,nAll,1,Npt) .* reshape(aCat_all.',1,nAll,1) ...
          + reshape(aCat_all,nAll,1,1)   .* reshape(dG_full_T,1,nAll,Npt);

H_full_all = H_G_full .* reshape(invcat,1,1,Npt) ...
           - crossTerm .* reshape(invcat2,1,1,Npt) ...
           + (aCat_all*aCat_all.') .* reshape(2*G_scaled .* invcat3,1,1,Npt);

% -------------------------------------------------------------------------
% Hessian wrt independent endmember variables
% -------------------------------------------------------------------------
H_c = pagemtimes(RT3, pagemtimes(H_full_all, R3));   % nVar x nVar x N

% -------------------------------------------------------------------------
% Hessian wrt elemental variables
% -------------------------------------------------------------------------
IpagesE = repmat(eye(nElemInd),1,1,Npt);
Jinv    = pagemldivide(Je_red, IpagesE);
H_e     = pagemtimes(permute(Jinv,[2 1 3]), pagemtimes(H_c, Jinv));

% -------------------------------------------------------------------------
% chi
% -------------------------------------------------------------------------
IpagesChi = repmat(eye(nElemInd),1,1,Npt);
chi       = pagemldivide(H_e, IpagesChi);

% -------------------------------------------------------------------------
% elemental concentration output
% -------------------------------------------------------------------------
e_out = (c_full.' * n_all);
e_out = (e_out ./ sum(e_out,2)).';
e_out = e_out(1:end-1,:);

% -------------------------------------------------------------------------
% collect results
% -------------------------------------------------------------------------
Result.g      = g;
Result.mu_c   = num2cell(mu_c,2)';
Result.mu_e   = num2cell(mu_e,2)';
Result.H_e    = H_e;
Result.H_c    = H_c;
Result.chi    = chi;
Result.e      = num2cell(e_out,2)';
Result.Jac    = Je_red;

% optional diagnostics
Result.c_full = c_full;
Result.c_real = c_real;
Result.c_pen  = c_pen;
end