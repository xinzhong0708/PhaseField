function varargout = LE_TailRegularization(action, varargin)
%LE_TAILREGULARIZATION Smooth tail regularization for KKS/LE composition.
%
% This function adds a p-dependent anchoring term for phase compositions:
%
%   F_tail_i = 1/2 * p_i * R_i(p_i) * |c_i - c_i_anchor|^2
%
% where R_i(p_i) is large for small phase fraction and fades smoothly to
% zero for larger phase fraction.
%
% Purpose:
%   Prevent tiny phase-field tails from carrying a fully free independent
%   composition, which can create c/mu/omega spikes near p -> 0.
%
% Usage inside LE_Calculator:
%
%   [Hreg, mu_eff] = LE_TailRegularization( ...
%       'apply', Hreg, mu_c, c{ip}, c_anchor{ip}, p_ip, reg, ip);
%
%   F = F + LE_TailRegularization('objective', pars, p, c, reg);
%
% Main options in reg:
%   reg.enable       : true/false
%   reg.c_anchor     : anchor composition, same structure as c
%   reg.p0           : below this p, full anchoring, default 1e-2
%   reg.p1           : above this p, no anchoring, default 5e-2
%   reg.lam_tail     : strength, scalar or vector over local phases
%   reg.scale_from_H : true/false, default true
%
% Recommended first test:
%   reg.enable       = true;
%   reg.p0           = 1e-2;
%   reg.p1           = 5e-2;
%   reg.lam_tail     = 1e-1;
%   reg.scale_from_H = true;

switch lower(action)

    case 'apply'
        [Hreg, mu_eff, info] = ApplyToQuadratic(varargin{:});
        varargout = {Hreg, mu_eff, info};

    case 'objective'
        Ftail = TailObjective(varargin{:});
        varargout = {Ftail};

    case 'weight'
        w = TailWeight(varargin{:});
        varargout = {w};

    otherwise
        error('Unknown action: %s', action);

end

end

% ============================================================
% Apply regularization to Hessian and chemical gradient
% ============================================================
function [Hreg, mu_eff, info] = ApplyToQuadratic( ...
    Hreg, mu_c, c_phase, c_anchor_phase, p_ip, reg, ip_local)

% Default: no change
mu_eff = mu_c;
info   = struct('enabled',false,'Rmax',0,'wmax',0);

if nargin < 7 || isempty(reg) || ~isfield(reg,'enable') || ~reg.enable
    return
end

if isempty(c_anchor_phase)
    return
end

[Nc,~,N] = size(Hreg);

p0       = GetPhaseValue(reg, 'p0',       ip_local, 1e-2);
p1       = GetPhaseValue(reg, 'p1',       ip_local, 5e-2);
lam_tail = GetPhaseValue(reg, 'lam_tail', ip_local, 0);

if lam_tail <= 0
    return
end

scale_from_H = GetField(reg, 'scale_from_H', true);

p_ip = reshape(p_ip,1,N);

% Tail weight:
%   w = 1 for p <= p0
%   w = 0 for p >= p1
%   smooth transition between
w = TailWeight(p_ip, p0, p1);

% Convert c and anchor into Nc x N arrays
C  = PhaseCellToMat(c_phase,        Nc, N);
Ca = PhaseCellToMat(c_anchor_phase, Nc, N);

% Hessian scale
if scale_from_H
    hscale = zeros(1,N);
    for i = 1:N
        H = 0.5*(Hreg(:,:,i) + Hreg(:,:,i).');
        hscale(i) = max(1, norm(H,'fro')/max(1,Nc));
    end
else
    hscale = ones(1,N);
end

% Tail ridge for each grid point
R = lam_tail .* hscale .* w;

% Add anchoring gradient:
%   mu_eff = mu_c + R * (c - c_anchor)
mu_eff = mu_c + (C - Ca) .* R;

% Add anchoring Hessian:
%   H_eff = H + R I
I = eye(Nc);
for i = 1:N
    if R(i) > 0
        Hreg(:,:,i) = Hreg(:,:,i) + R(i)*I;
    end
end

info.enabled = true;
info.Rmax    = max(R);
info.wmax    = max(w);

end

% ============================================================
% Objective contribution
% ============================================================
function Ftail = TailObjective(pars, p, c, reg)

Np = numel(c);
N  = numel(c{1}{1});

Ftail = zeros(1,N);

if nargin < 4 || isempty(reg) || ~isfield(reg,'enable') || ~reg.enable
    return
end

if ~isfield(reg,'c_anchor') || isempty(reg.c_anchor)
    return
end

scale_from_H = GetField(reg, 'scale_from_H', true);

for ip = 1:Np

    if ip > numel(reg.c_anchor) || isempty(reg.c_anchor{ip})
        continue
    end

    Nc = numel(c{ip});
    if Nc == 0
        continue
    end

    p0       = GetPhaseValue(reg, 'p0',       ip, 1e-2);
    p1       = GetPhaseValue(reg, 'p1',       ip, 5e-2);
    lam_tail = GetPhaseValue(reg, 'lam_tail', ip, 0);

    if lam_tail <= 0
        continue
    end

    p_ip = reshape(p(:,:,ip),1,N);
    w    = TailWeight(p_ip, p0, p1);

    if all(w == 0)
        continue
    end

    C  = PhaseCellToMat(c{ip},            Nc, N);
    Ca = PhaseCellToMat(reg.c_anchor{ip}, Nc, N);
    dC2 = sum((C - Ca).^2,1);

    if scale_from_H
        Rthermo = PhaseThermo(pars{ip}, c{ip});

        if isempty(Rthermo.H_c)
            continue
        end

        Hc = Rthermo.H_c;
        hscale = zeros(1,N);

        for i = 1:N
            H = 0.5*(Hc(:,:,i) + Hc(:,:,i).');
            hscale(i) = max(1, norm(H,'fro')/max(1,Nc));
        end
    else
        hscale = ones(1,N);
    end

    R = lam_tail .* hscale .* w;

    % Important:
    % The energy term has p_i * R, but the stationarity equation
    % divides out p_i for active phases, giving gradient R*(c-ca).
    Ftail = Ftail + 0.5 .* p_ip .* R .* dC2;

end

end

% ============================================================
% Smooth tail weight
% ============================================================
function w = TailWeight(p, p0, p1)
%TAILWEIGHT
%
% w = 1 for p <= p0
% w = 0 for p >= p1
% smoothstep transition between.

if p1 <= p0
    error('LE_TailRegularization: p1 must be larger than p0.');
end

x = (p - p0) ./ max(p1 - p0, eps);
x = min(max(x,0),1);

s = x.^2 .* (3 - 2*x);  % smoothstep
w = 1 - s;

end

% ============================================================
% Helpers
% ============================================================
function C = PhaseCellToMat(c_phase, Nc, N)

C = zeros(Nc,N);

for ic = 1:Nc
    tmp = c_phase{ic};
    C(ic,:) = reshape(tmp,1,N);
end

end

function val = GetPhaseValue(S, field, ip, default)

if ~isfield(S,field) || isempty(S.(field))
    val = default;
    return
end

x = S.(field);

if isscalar(x)
    val = x;
else
    ip = min(ip, numel(x));
    val = x(ip);
end

end

function val = GetField(S, field, default)

if isfield(S,field) && ~isempty(S.(field))
    val = S.(field);
else
    val = default;
end

end