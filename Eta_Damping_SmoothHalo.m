function eta_eff = Eta_Damping_SmoothHalo(pAct, etaBulk, etaInt, q2, p02, q3, p03, etaThree, nsmooth, etaHalo, halo_cut)
%ETA_DAMPING_SMOOTHHALO Smooth interface eta damping.
%
% pAct: ny x nx x Nphase, preferably collapsed thermodynamic phase p
%
% Idea:
%   1) raw trigger:
%        small p2 immediately lowers eta to etaInt
%   2) smoothed halo:
%        eta is also reduced slightly into nearby grains
%
% This avoids a sharp 20x jump in 1/eta when etaInt = 0.05*etaBulk.
%
% Recommended:
%   etaBulk  = eta
%   etaInt   = 0.05*eta
%   etaHalo  = 0.20*eta
%   etaThree = 0.03*eta
%   q2,q3    = 4
%   p02,p03  = 3e-3 or 1e-3
%   nsmooth  = 2
%   halo_cut = 1e-4

if nargin < 4 || isempty(q2),        q2 = 4; end
if nargin < 5 || isempty(p02),      p02 = 3e-3; end
if nargin < 6 || isempty(q3),        q3 = 4; end
if nargin < 7 || isempty(p03),      p03 = 3e-3; end
if nargin < 8 || isempty(etaThree), etaThree = etaInt; end
if nargin < 9 || isempty(nsmooth),  nsmooth = 2; end
if nargin < 10, etaHalo = []; end
if nargin < 11 || isempty(halo_cut), halo_cut = 1e-4; end

ps = sort(max(pAct,0),3,'descend');

p1 = ps(:,:,1);

if size(ps,3) >= 2
    p2 = ps(:,:,2);
else
    p2 = zeros(size(p1));
end

if size(ps,3) >= 3
    p3 = ps(:,:,3);
else
    p3 = zeros(size(p1));
end

% ------------------------------------------------------------
% Expand scalar eta inputs
% ------------------------------------------------------------
if isscalar(etaBulk),  etaBulk  = etaBulk  * ones(size(p1)); end
if isscalar(etaInt),   etaInt   = etaInt   * ones(size(p1)); end
if isscalar(etaThree), etaThree = etaThree * ones(size(p1)); end

if isempty(etaHalo)
    etaHalo = sqrt(etaBulk.*etaInt);
elseif isscalar(etaHalo)
    etaHalo = etaHalo * ones(size(p1));
end

% Make sure damping never gives larger eta than bulk
etaInt   = min(etaInt,etaBulk);
etaHalo  = min(etaHalo,etaBulk);
etaThree = min(etaThree,etaBulk);

% ------------------------------------------------------------
% Local raw trigger
% ------------------------------------------------------------
w2_raw = (p2.^q2) ./ (p2.^q2 + p02^q2);
w3_raw = (p3.^q3) ./ (p3.^q3 + p03^q3);

% ------------------------------------------------------------
% Smooth halo trigger
% ------------------------------------------------------------
w2_halo = Smooth01(w2_raw,nsmooth);
w3_halo = Smooth01(w3_raw,nsmooth);

w2_halo(w2_halo < halo_cut) = 0;
w3_halo(w3_halo < halo_cut) = 0;

% ------------------------------------------------------------
% Interpolate in k = 1/eta
% ------------------------------------------------------------
kBulk  = 1 ./ etaBulk;
kHalo  = 1 ./ etaHalo;
kInt   = 1 ./ etaInt;
kThree = 1 ./ etaThree;

etaThreeHalo = sqrt(etaBulk.*etaThree);
etaThreeHalo = min(etaThreeHalo,etaBulk);
kThreeHalo   = 1 ./ etaThreeHalo;

% Start from bulk
kEff = kBulk;

% Gentle two-phase halo into the grains
kEff = kEff + max(kHalo - kEff,0).*w2_halo;

% Strong local two-phase damping where second phase is really present
kEff = kEff + max(kInt - kEff,0).*w2_raw;

% Gentle three-phase halo
kEff = kEff + max(kThreeHalo - kEff,0).*w3_halo;

% Strong local three-phase damping
kEff = kEff + max(kThree - kEff,0).*w3_raw;

eta_eff = 1 ./ kEff;

end


function A = Smooth01(A,niter)

if niter <= 0
    A = min(max(A,0),1);
    return
end

ker = [1 2 1; 2 4 2; 1 2 1];
ker = ker/sum(ker(:));

for it = 1:niter
    A = conv2(A,ker,'same');
end

A = min(max(A,0),1);

end