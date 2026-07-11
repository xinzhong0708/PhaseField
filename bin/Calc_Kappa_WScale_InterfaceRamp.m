function PARAM = Calc_Kappa_WScale_InterfaceRamp(STATE,PARAM,MODEL,PHYS)
%CALC_KAPPA_WSCALE_INTERFACERAMP Build WScale and kappa_eff maps.

[ny,nx,Ng] = size(STATE.p);

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
tail_p     = GetParam(PARAM,{'WScale_tail_p_presence','WScale_p_presence'},1e-8);
w_dom_p    = GetParam(PARAM,{'WScale_domain_p'},0.5);
k_dom_p    = GetParam(PARAM,{'kappa_domain_p'},w_dom_p);
zero_w     = GetParam(PARAM,{'ramp_zero_width','kappa_zero_width'},2);
ramp_w     = GetParam(PARAM,{'ramp_width','kappa_ramp_width'},6);
kmin       = GetParam(PARAM,{'kappa_min_frac'},0.02);
fill_iter  = GetParam(PARAM,{'WScale_domain_fill_iter'},2);

tail_p    = max(tail_p,0);
w_dom_p   = min(max(w_dom_p,0),1);
k_dom_p   = min(max(k_dom_p,0),1);
zero_w    = max(0,round(zero_w));
ramp_w    = max(1,round(ramp_w));
kmin      = min(max(kmin,0),1);
fill_iter = max(0,round(fill_iter));

% -------------------------------------------------------------------------
% Collapse grains to thermodynamic phases
% -------------------------------------------------------------------------
[phase_id,grain_to_phase,p_phase,phase_name] = CollapsePhases(STATE,MODEL);

Nphase = numel(phase_id);

% -------------------------------------------------------------------------
% Phase WScale factor
% -------------------------------------------------------------------------
if isfield(PARAM,'WScale_phase_factor') && ~isempty(PARAM.WScale_phase_factor)
    wmin_phase = PhaseValue(PARAM.WScale_phase_factor,phase_id, ...
        grain_to_phase,Nphase,Ng,1,'PARAM.WScale_phase_factor');
else
    wmin_phase = ones(1,Nphase);
end

wmin_phase = min(max(wmin_phase,0),1);
active     = wmin_phase < 1 - 1e-12;

% -------------------------------------------------------------------------
% Phase kappa value
%
% Activation is controlled only by WScale_phase_factor.
% If PARAM.kappa_phase_value is absent, active phases use PHYS.kappa.
% -------------------------------------------------------------------------
if isfield(PARAM,'kappa_phase_value') && ~isempty(PARAM.kappa_phase_value)
    kappa_value = PhaseValue(PARAM.kappa_phase_value,phase_id, ...
        grain_to_phase,Nphase,Ng,PHYS.kappa,'PARAM.kappa_phase_value');
else
    kappa_value = PHYS.kappa*ones(1,Nphase);
end

kappa_phase = zeros(1,Nphase);
kappa_phase(active) = kappa_value(active);

% -------------------------------------------------------------------------
% Allocate maps
% -------------------------------------------------------------------------
kappa_eff     = zeros(ny,nx);
kappa_weight  = zeros(ny,nx,Nphase);
recovery_map  = zeros(ny,nx,Nphase);
w_scale_phase = ones(ny,nx,Nphase);

w_tail_map    = false(ny,nx,Nphase);
w_domain_mask = false(ny,nx,Nphase);
k_domain_mask = false(ny,nx,Nphase);

% -------------------------------------------------------------------------
% Build maps
% -------------------------------------------------------------------------
for iph = 1:Nphase

    pcur = p_phase(:,:,iph);
    wmin = wmin_phase(iph);
    kp   = kappa_phase(iph);

    w_dom = FillSmallHolesNoWrap(pcur > w_dom_p,fill_iter);
    k_dom = FillSmallHolesNoWrap(pcur > k_dom_p,fill_iter);

    recovery = BuildInsideRamp(w_dom,zero_w,ramp_w,ny,nx);

    ws = ones(ny,nx);
    wk = zeros(ny,nx);

    if active(iph)

        tail = (pcur > tail_p) & ~w_dom;

        ws(tail)  = wmin;
        ws(w_dom) = wmin + (1-wmin).*recovery(w_dom);

        if kp ~= 0
            wk(k_dom) = kmin + (1-kmin).*recovery(k_dom);
        end

        w_tail_map(:,:,iph) = tail;
    end

    kappa_eff = max(kappa_eff,kp.*wk);

    w_scale_phase(:,:,iph) = ws;
    kappa_weight(:,:,iph)  = wk;
    recovery_map(:,:,iph)  = recovery;
    w_domain_mask(:,:,iph) = w_dom;
    k_domain_mask(:,:,iph) = k_dom;
end

% -------------------------------------------------------------------------
% Store outputs
% -------------------------------------------------------------------------
PARAM.kappa_eff              = kappa_eff;
PARAM.kappa_weight_phase     = kappa_weight;
PARAM.kappa_weight_phase_raw = recovery_map;
PARAM.kappa_phase_collapsed  = kappa_phase;
PARAM.kappa_p_phase          = p_phase;

PARAM.w_scale_phase          = w_scale_phase;
PARAM.w_scale_active_phase   = active;
PARAM.w_scale_min_phase      = wmin_phase;
PARAM.w_scale_min_used       = wmin_phase;
PARAM.WScale_phase_factor_used = wmin_phase;

PARAM.ramp_weight_phase      = recovery_map;
PARAM.recovery_weight_phase  = recovery_map;
PARAM.purity_weight_phase    = ones(ny,nx,Nphase);

PARAM.w_tail_map             = w_tail_map;
PARAM.w_domain_mask          = w_domain_mask;
PARAM.kappa_domain_mask      = k_domain_mask;

PARAM.phase_id_collapsed     = phase_id;
PARAM.phase_name_collapsed   = phase_name;
PARAM.grain_to_phase         = grain_to_phase;

PARAM.WScale_tail_p_presence_used = tail_p;
PARAM.WScale_domain_p_used        = w_dom_p;
PARAM.kappa_domain_p_used         = k_dom_p;
PARAM.ramp_zero_width_used        = zero_w;
PARAM.ramp_width_used             = ramp_w;
PARAM.kappa_min_frac_used         = kmin;

% Compatibility fields. Eta is not changed here.
PARAM.eta_interface_scale = ones(ny,nx);
PARAM.eta_scale_phase     = ones(ny,nx,Nphase);

PARAM.use_WScale = any(active);

if any(kappa_phase ~= 0)
    mask = kappa_phase ~= 0;
    PARAM.kappa_raw_weight = max(recovery_map(:,:,mask),[],3);
    PARAM.kappa_chi_weight = max(kappa_weight(:,:,mask),[],3);
else
    PARAM.kappa_raw_weight = zeros(ny,nx);
    PARAM.kappa_chi_weight = zeros(ny,nx);
end

end

% =========================================================================
% Collapse grains to thermodynamic phases
% =========================================================================
function [phase_id,grain_to_phase,p_phase,phase_name] = CollapsePhases(STATE,MODEL)

[ny,nx,Ng] = size(STATE.p);

if isfield(MODEL,'phase_index') && ~isempty(MODEL.phase_index)
    phase_index = MODEL.phase_index(:).';
else
    phase_index = 1:Ng;
end

if numel(phase_index) ~= Ng
    error('MODEL.phase_index length must match number of grains.')
end

phase_id = unique(phase_index,'stable');
Nphase   = numel(phase_id);

grain_to_phase = zeros(1,Ng);
p_phase        = zeros(ny,nx,Nphase);
phase_name     = cell(1,Nphase);

for iph = 1:Nphase

    grains = find(phase_index == phase_id(iph));

    grain_to_phase(grains) = iph;
    p_phase(:,:,iph) = sum(STATE.p(:,:,grains),3);
    phase_name{iph}  = GetPhaseName(MODEL,phase_id(iph),iph);

end

end

% =========================================================================
% Convert scalar / grain vector / phase vector to collapsed phase vector
% =========================================================================
function v_phase = PhaseValue(v_in,phase_id,grain_to_phase,Nphase,Ng,default_value,name)

if isempty(v_in)
    v_phase = default_value*ones(1,Nphase);
    return
end

v = v_in;

if iscell(v)
    v = cellfun(@(x) x(1),v);
end

v = v(:).';

if isscalar(v)

    v_phase = v*ones(1,Nphase);

elseif numel(v) >= max(phase_id) && all(phase_id == round(phase_id))

    % Preferred case for PARAM.WScale_phase_factor:
    % one value per metadata phase, indexed by original phase_id.
    v_phase = v(phase_id);

elseif numel(v) == Nphase

    v_phase = v;

elseif numel(v) == Ng

    v_phase = default_value*ones(1,Nphase);

    for iph = 1:Nphase
        grains = find(grain_to_phase == iph);

        if ~isempty(grains)
            v_phase(iph) = min(v(grains));
        end
    end

else

    error('%s must be scalar, length Ng, length Nphase, or metadata phase length.',name)

end

end

% =========================================================================
% Build inside-domain ramp
% =========================================================================
function w = BuildInsideRamp(domain,zero_width,ramp_width,ny,nx)

if ~any(domain(:))
    w = zeros(ny,nx);
    return
end

if all(domain(:))
    w = ones(ny,nx);
    return
end

max_dist = zero_width + ramp_width + 1;

dist  = inf(ny,nx);
known = ~domain;

dist(known) = 0;

for ir = 1:max_dist
    known_new = DilateNoWrap(known);
    add = known_new & ~known;

    dist(add) = ir;
    known = known_new;
end

dist(~isfinite(dist)) = max_dist;

x = (dist - zero_width)./ramp_width;
x = min(max(x,0),1);

w = x.^2.*(3 - 2*x);
w(~domain) = 0;

end

% =========================================================================
% Optional parameter reader
% =========================================================================
function val = GetParam(PARAM,names,default_value)

val = default_value;

for i = 1:numel(names)
    fn = names{i};

    if isfield(PARAM,fn) && ~isempty(PARAM.(fn))
        val = PARAM.(fn);
        return
    end
end

end

% =========================================================================
% Phase name
% =========================================================================
function pname = GetPhaseName(MODEL,pid,iph)

pname = ['phase_',num2str(pid)];

if isfield(MODEL,'phs_name')
    if numel(MODEL.phs_name) >= pid
        pname = MODEL.phs_name{pid};
        return
    elseif numel(MODEL.phs_name) >= iph
        pname = MODEL.phs_name{iph};
        return
    end
end

if isfield(MODEL,'phase_name')
    if numel(MODEL.phase_name) >= pid
        pname = MODEL.phase_name{pid};
        return
    elseif numel(MODEL.phase_name) >= iph
        pname = MODEL.phase_name{iph};
        return
    end
end

end

% =========================================================================
% No-wrap morphology
% =========================================================================
function A = FillSmallHolesNoWrap(A,niter)

A = logical(A);

for it = 1:niter
    A = ErodeNoWrap(DilateNoWrap(A));
end

end

function B = DilateNoWrap(A)

[ny,nx] = size(A);
B = A;

if ny > 1
    B(2:ny,:)   = B(2:ny,:)   | A(1:ny-1,:);
    B(1:ny-1,:) = B(1:ny-1,:) | A(2:ny,:);
end

if nx > 1
    B(:,2:nx)   = B(:,2:nx)   | A(:,1:nx-1);
    B(:,1:nx-1) = B(:,1:nx-1) | A(:,2:nx);
end

end

function B = ErodeNoWrap(A)

[ny,nx] = size(A);
B = A;

if ny > 1
    B(2:ny,:)   = B(2:ny,:)   & A(1:ny-1,:);
    B(1:ny-1,:) = B(1:ny-1,:) & A(2:ny,:);
end

if nx > 1
    B(:,2:nx)   = B(:,2:nx)   & A(:,1:nx-1);
    B(:,1:nx-1) = B(:,1:nx-1) & A(:,2:nx);
end

end