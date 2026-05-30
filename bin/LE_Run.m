function STATE = LE_Run(STATE,PARAM,MODEL)
%LE_RUN Local-equilibrium update using structured variables.
%
% Key choices:
%   - p_th smooths small geometric phase tails for LE.
%   - active selects thermodynamic phase subsets with hysteresis.
%   - p_le renormalizes retained phases to sum to one; E remains unchanged
%     because it is the conserved bulk-composition field evolved by CH.
%   - PARAM.use_WScale controls interface W scaling in multiphase LE:
%         0 = retain original excess-energy thermodynamics;
%         1 = damp W in interface LE calculations only.
%   - omega is always calculated from original pars with excess energy.

%Unpack structured variables
pars        = MODEL.pars;
phase_index = MODEL.phase_index(:).';

p           = STATE.p;
c           = STATE.c;
E           = STATE.E;
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
nx          = size(p,2);

c           = Unpack_c(c);
E           = Unpack_E(E);
p           = Unpack_p(p);
mu_e        = Unpack_E(mu_e);
chi         = Unpack_Chi(chi);

%Collapse repeated grains to thermodynamic phases
[pars,c,p,grain_to_phase] = Collapse_LE_Phases(pars,c,p,phase_index);

%Local-equilibrium controls
Pmax        = 4;

p_tail      = 1e-3;
p_full      = 5e-2;
p_on        = 5e-3;
p_off       = 2e-3;       %Must be smaller than p_on for hysteresis

%Damping factors for 1-, 2-, 3- and 4-phase coexistence
alpha_LE    = [0.8  0.5   0.4   0.3 ];
iter_LE     = [100  1000  1000  1000];

Np          = numel(c);
N           = numel(c{1}{1});

%Thermodynamic weights before active-set filtering
p_th        = Calc_Thermo_p(p,p_tail,p_full);
pThMat      = reshape(p_th,N,Np);

%Hysteretic active set
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

%Process only phase combinations present in the grid
[active_sets,~,set_id] = unique(active,'rows');

for iset = 1:size(active_sets,1)

    ph_act = find(active_sets(iset,:));
    k      = numel(ph_act);
    mask   = (set_id == iset).';

    c_cur  = Slice_c(c,ph_act,mask);
    E_cur  = Slice_E(E,mask);
    eta    = eta_vec(mask);
    p_cur  = p_le(:,mask,ph_act);

    if PARAM.use_WScale && k > 1
        pars_cur = Slice_Pars_WScale(pars_inter,ph_act,mask);
    else
        pars_cur = pars(ph_act);
    end

    % [c_tmp,mu_tmp,chi_tmp] = LE_Calculator(pars_cur,p_cur,c_cur,E_cur,eta,[alpha_LE(k),iter_LE(k)]);


    PARAM.use_mixed_local_LE = 1;
    mu_cur                   = Slice_E(mu_e,mask);

    if isfield(PARAM,'use_mixed_local_LE') && PARAM.use_mixed_local_LE == 1
        [c_tmp,mu_tmp,chi_tmp] = LE_Calculator_MixedLocal(pars_cur,p_cur,c_cur,E_cur,mu_cur,eta,[alpha_LE(k),iter_LE(k),0,0]);
    else
        [c_tmp,mu_tmp,chi_tmp] = LE_Calculator(pars_cur,p_cur,c_cur,E_cur,eta,[alpha_LE(k),iter_LE(k)]);
    end




    [c,mu_e,chi] = Assign_LE_Back(c,mu_e,chi,c_tmp,mu_tmp,chi_tmp,ph_act,mask);

end

%Pack collapsed thermodynamic phases
c_col   = Pack_c(c,ny);
mu_e    = Pack_E(mu_e,ny);
chi     = Pack_chi(chi,ny);
e_col   = Calc_e(pars,c_col);

%Calculate grand potentials from original excess-energy thermodynamics
Ne      = numel(E);
omg_col = zeros(ny,nx,Np);

for ip = 1:Np
    omg_col(:,:,ip) = reshape(PhaseG(pars{ip},c{ip}),ny,[]);
    for ie = 1:Ne
        omg_col(:,:,ip) = omg_col(:,:,ip) - e_col{ip}{ie} .* mu_e{ie};
    end
end

%Copy collapsed results back to grain-resolved variables
Ngrain  = numel(grain_to_phase);
c_out   = cell(1,Ngrain);
omg     = zeros(ny,nx,Ngrain);

for ig = 1:Ngrain
    iph         = grain_to_phase(ig);
    c_out{ig}   = c_col{iph};
    omg(:,:,ig) = omg_col(:,:,iph);
end

%Store collapsed active-set information
LE_state.phase_index    = phase_index;
LE_state.grain_to_phase = grain_to_phase;

%Write back to STATE
STATE.c        = c_out;
STATE.e        = Calc_e(MODEL.pars,c_out);
STATE.mu_e     = mu_e;
STATE.chi      = chi;
STATE.omg      = omg;
STATE.LE_state = LE_state;

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

x     = (p - p_tail) ./ (p_full - p_tail);
x     = min(max(x,0),1);

w     = p .* x.^2 .* (3 - 2*x);
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
    iph        = grain_to_phase(ig);
    p_c(:,:,iph) = p_c(:,:,iph) + p(:,:,ig);
end

%Collapse c by p-weighted averaging
c_c = cell(1,Nphase);

for iph = 1:Nphase

    grains     = find(grain_to_phase == iph);
    ig0        = grains(1);
    Nc         = numel(c{ig0});
    c_c{iph}   = cell(1,Nc);
    den        = reshape(sum(p(:,:,grains),3),size(c{ig0}{1}));
    good       = den > eps;

    for ic = 1:Nc

        num = zeros(size(c{ig0}{ic}));

        for ig = grains
            num = num + reshape(p(:,:,ig),size(num)) .* c{ig}{ic};
        end

        tmp        = c{ig0}{ic};
        tmp(good)  = num(good) ./ den(good);
        c_c{iph}{ic} = tmp;

    end

end

end