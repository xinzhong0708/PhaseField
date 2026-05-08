function eta_eff = Eta_Damping(pAct, etaBulk, etaInt, q2, p02, q3, p03, etaThree)
% pAct: ny x nx x Np
% etaBulk : bulk eta
% etaInt  : eta when a 2nd phase is present
% q2,p02  : sharpness and onset for 2-phase trigger
% q3,p03  : sharpness and onset for 3-phase trigger
% etaThree: eta when a true 3rd phase is present
%
% Recommended:
%   etaBulk  = eta
%   etaInt   = 0.3*eta or 0.5*eta
%   q2       = 4 or 6
%   p02      = 1e-2
%   q3       = 4 or 6
%   p03      = 1e-2
%   etaThree = 0.2*eta or 0.3*eta

if nargin < 4 || isempty(q2),             q2 = 4; end
if nargin < 5 || isempty(p02),           p02 = 1e-2; end
if nargin < 6 || isempty(q3),             q3 = 4; end
if nargin < 7 || isempty(p03),           p03 = 1e-2; end
if nargin < 8 || isempty(etaThree), etaThree = etaInt; end

ps = sort(max(pAct,0), 3, 'descend');

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

% 2-phase trigger:
% as soon as the 2nd phase appears at ~p02, eta starts dropping strongly
w2 = (p2.^q2) ./ (p2.^q2 + p02^q2);

% 3-phase trigger:
% as soon as the 3rd phase appears at ~p03, move further toward etaThree
w3 = (p3.^q3) ./ (p3.^q3 + p03^q3);

% allow scalar eta inputs
if isscalar(etaBulk),  etaBulk  = etaBulk  * ones(size(p1)); end
if isscalar(etaInt),   etaInt   = etaInt   * ones(size(p1)); end
if isscalar(etaThree), etaThree = etaThree * ones(size(p1)); end

% interpolate in kappa = 1/eta, because that is what enters the LE operator
kBulk  = 1 ./ etaBulk;
kInt   = 1 ./ etaInt;
kThree = 1 ./ etaThree;

% first: if a 2nd phase is present, soften toward etaInt
kEff = kBulk + (kInt - kBulk) .* w2;

% then: if a 3rd phase is present, soften further toward etaThree
kEff = kEff + (kThree - kEff) .* w3;

eta_eff = 1 ./ kEff;
end