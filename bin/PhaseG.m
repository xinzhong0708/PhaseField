function [g] = PhaseG(pars,c)
% Thermodynamic calculator to provide g only
%
% Full endmember order is ALWAYS:
%   [ real_1 ... real_nReal , pen_1 ... pen_nPen ]
%
% Dependent endmember is ALWAYS the LAST full endmember.
%
% Input c can be either:
%   (1) full composition, size = nAll x N
%   (2) independent composition, size = (nAll-1) x N

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
    error('PhaseG: pars.gN and pars.nN size mismatch.');
end
if size(n_pen,2) ~= size(n_real,2)
    error('PhaseG: pars.nN must have same number of columns as pars.n.');
end

% -------------------------------------------------------------------------
% construct full c
% dependent is always the LAST full endmember
% -------------------------------------------------------------------------
if size(c,1) == nAll
    c_full = c;
elseif size(c,1) == nAll-1
    c_full = [c; 1-sum(c,1)];
else
    error('PhaseG: input c must have %d or %d rows, got %d.', ...
          nAll-1, nAll, size(c,1));
end

% split into real / penalty
c_real = c_full(1:nReal,:);           % nReal x N
c_pen  = c_full(nReal+1:end,:);       % nPen  x N

c_realt = c_real.';
c_pent  = c_pen.';

% -------------------------------------------------------------------------
% real-solution-model free energy
% -------------------------------------------------------------------------
% mechanical mixing over REAL endmembers
g_mech_real = c_realt * g0;           % N x 1

% regularized ideal mixing term
z     = c_realt * zt + 1e-30;

eps0  = 1e-6;
sz    = sqrt(z.^2 + eps0.^2);
lz    = log(z + sz) - log(2);

Azt   = zt .* log(zt + double(zt==0));

g_id  = RT * sum(mtpl .* ( z .* lz - c_realt * Azt ), 2);

% nonideal term (REAL only)
alp_eff = [1 T P] * alp;                         % 1 x nReal
W       = w(:,:,1) + w(:,:,2)*T + w(:,:,3)*P;   % nReal x nReal

a       = alp_eff(:);                            % nReal x 1
M       = ((a*a.') ./ (a + a.')) .* W;           % nReal x nReal
% Optional local scaling of Margules/nonideal term
if isfield(pars,'w_scale') && ~isempty(pars.w_scale)
    w_scale = pars.w_scale;
    if isscalar(w_scale)
        w_scale = w_scale * ones(size(c_full,2),1);
    else
        w_scale = w_scale(:);
    end
    w_scale = min(max(w_scale,0),1);
else
    w_scale = ones(size(c_full,2),1);
end

q       = c_realt * a;                           % N x 1
invq    = 1 ./ q;

Mc      = c_realt * M;                           % N x nReal
n_nid   = sum(Mc .* c_realt, 2);                 % N x 1
g_nid   = w_scale .* (n_nid .* invq);            % N x 1
% -------------------------------------------------------------------------
% penalty-endmember extra energy
%   g_pen*c_pen + penalty*sum(c_pen^2)
% -------------------------------------------------------------------------
if nPen > 0
    g_pen_extra = c_pent * g_pen + penalty * sum(c_pent.^2, 2);   % N x 1
else
    g_pen_extra = zeros(size(c_full,2),1);
end

% -------------------------------------------------------------------------
% total raw free energy before normalization
% -------------------------------------------------------------------------
G_raw = g_mech_real + g_id + g_nid + g_pen_extra;   % N x 1

% -------------------------------------------------------------------------
% normalize to one cation using ALL endmembers
% -------------------------------------------------------------------------
n_all    = [n_real; n_pen];         % nAll x NeCat
aCat_all = sum(n_all,2);            % nAll x 1

cation   = c_full.' * aCat_all;     % N x 1
invcat   = 1 ./ cation;

scale    = 1 / (pars.E_sc * pars.vref);
G_scaled = G_raw * scale;           % N x 1
g        = G_scaled .* invcat;      % N x 1

end