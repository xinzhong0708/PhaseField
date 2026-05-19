function STATE = LE_Run(STATE,PARAM,MODEL)
%LE_RUN Local-equilibrium update using structured variables.

%Unpack structured variables
pars        = MODEL.pars;
phase_index = MODEL.phase_index(:).';

p           = STATE.p;
c           = STATE.c;
E           = STATE.E;
mu_e        = STATE.mu_e;
chi         = STATE.chi;
eta_vec     = PARAM.eta(:);
LE_state    = STATE.LE_state;

if isempty(LE_state)
    LE_state = struct();
end

%Reshape 2D into 1D
nx          = size(p,2);
ny          = size(p,1);
c           = Unpack_c(c);
E           = Unpack_E(E);
p           = Unpack_p(p);
mu_e        = Unpack_E(mu_e);
chi         = Unpack_Chi(chi);

%Keep grain-resolved copies
pars_full   = pars;
c_full      = c;
p_full      = p;

%Collapse repeated grains to thermodynamic phases
[pars,c,p,grain_to_phase] = Collapse_LE_Phases(pars_full,c_full,p_full,phase_index);

%Maximal number of allowed phase coexist
Pmax    = 4;

%Thermodynamic interpolation thresholds
p_tail  = 1e-4;
p_full  = 5e-3;

p_on    = 5e-4;
p_off   = 5e-5;

p_th    = Calc_Thermo_p(p,p_tail,p_full);

%Interface w damping
dp1     = 0.85;
dp2     = 0.99;

Np      = numel(c);
N       = numel(c{1}{1});

%Build raw p matrix and thermodynamic p matrix
pMat    = zeros(N,Np);
pThMat  = zeros(N,Np);

for ip = 1:Np

    p_ip          = squeeze(p(:,:,ip));
    pMat(:,ip)   = p_ip(:);

    p_th_ip       = squeeze(p_th(:,:,ip));
    pThMat(:,ip) = p_th_ip(:);

end

%Initialize hysteresis active set
if isempty(LE_state) || ~isfield(LE_state,'active') || isempty(LE_state.active) || ~isequal(size(LE_state.active),[N,Np])

    LE_state.active = pThMat > p_on;

    % Make sure every grid point has at least one active phase
    for i = 1:N

        if ~any(LE_state.active(i,:))

            [~,idmax] = max(pThMat(i,:));

            if pThMat(i,idmax) <= 0
                [~,idmax] = max(pMat(i,:));
            end

            LE_state.active(i,idmax) = true;

        end

    end

end

%Update hysteresis active set
active_old = LE_state.active;
active_new = false(N,Np);

for i = 1:N

    for ip = 1:Np

        if active_old(i,ip)

            % Active phase remains active until thermodynamic weight drops below p_off
            active_new(i,ip) = pThMat(i,ip) > p_off;

        else

            % Inactive phase becomes active only when thermodynamic weight rises above p_on
            active_new(i,ip) = pThMat(i,ip) > p_on;

        end

    end

    %Make sure every grid point has at least one active phase
    if ~any(active_new(i,:))

        [~,idmax] = max(pThMat(i,:));

        if pThMat(i,idmax) <= 0
            [~,idmax] = max(pMat(i,:));
        end

        active_new(i,idmax) = true;

    end

    %Limit maximum number of active phases
    if sum(active_new(i,:)) > Pmax

        score = pThMat(i,:);

        %Favor phases that were already active
        score(active_old(i,:)) = score(active_old(i,:)) + 0.5*p_on;

        [~,ord] = sort(score,'descend');

        tmp = false(1,Np);
        tmp(ord(1:Pmax)) = true;

        active_new(i,:) = tmp;

    end

end

LE_state.active = active_new;

%Store for diagnostics
LE_state.p_th = p_th;

%Handle all active subsets of size 1,2,3,4
for k = 1:min([Np,Pmax])

    combs = nchoosek(1:Np, k);

    for icomb = 1:size(combs,1)

        ph_act = combs(icomb,:);

        %Mask: active phases follow hysteresis state
        mask = true(1,N);

        for ip = 1:Np

            if ismember(ip, ph_act)

                mask = mask & LE_state.active(:,ip).';

            else

                mask = mask & ~LE_state.active(:,ip).';

            end

        end

        %Only calculate when mask exist
        if ~any(mask); continue; end

        %Slice local fields
        c_cur = Slice_c(c, ph_act, mask);
        E_cur = Slice_E(E, mask);
        eta   = eta_vec(mask);

        %One phase
        if k == 1
            ip                     = ph_act(1);
            [c_tmp,mu_tmp,chi_tmp] = LE_Calculator(pars(ip), ones(size(c_cur{1}{1})), c_cur(1), E_cur, eta, [0.8,100]);
            [c, mu_e, chi]         = Assign_LE_Back(c, mu_e, chi, c_tmp, mu_tmp, chi_tmp, ip, mask);
        end

        %Two phase
        if k == 2
            p_cur                  = p_th(:,mask,ph_act);
            pars_inter             = Apply_WScale_FromP(pars,p(:,mask,:),dp1,dp2,1);
            % pars_inter             = pars;
            [c_tmp,mu_tmp,chi_tmp] = LE_Calculator(pars_inter(ph_act), p_cur, c_cur, E_cur, eta, [0.6,1000]);
            [c, mu_e, chi]         = Assign_LE_Back(c, mu_e, chi, c_tmp, mu_tmp, chi_tmp, ph_act, mask);
        end

        %Three phase
        if k == 3
            p_cur                  = p_th(:,mask,ph_act);
            pars_inter             = Apply_WScale_FromP(pars,p(:,mask,:),dp1,dp2,1);
            % pars_inter             = pars;
            [c_tmp,mu_tmp,chi_tmp] = LE_Calculator(pars_inter(ph_act), p_cur, c_cur, E_cur, eta, [0.4,1000]);
            [c, mu_e, chi]         = Assign_LE_Back(c, mu_e, chi, c_tmp, mu_tmp, chi_tmp, ph_act, mask);
        end

        %Four phase
        if k == 4
            p_cur                  = p_th(:,mask,ph_act);
            pars_inter             = Apply_WScale_FromP(pars,p(:,mask,:),dp1,dp2,1);
            % pars_inter             = pars;
            [c_tmp,mu_tmp,chi_tmp] = LE_Calculator(pars_inter(ph_act), p_cur, c_cur, E_cur, eta, [0.2,1000]);
            [c, mu_e, chi]         = Assign_LE_Back(c, mu_e, chi, c_tmp, mu_tmp, chi_tmp, ph_act, mask);
        end

    end

end

%Pack up collapsed thermodynamic phases
g = cell(1,Np);

for ip = 1:Np
    g{ip} = reshape(PhaseG(pars{ip},c{ip}),ny,[]);
end

c_col = Pack_c(c,ny);
mu_e  = Pack_E(mu_e,ny);
chi   = Pack_chi(chi,ny);
e_col = Calc_e(pars,c_col);

%Calculate omega for collapsed phases
Ne      = length(E);
omg_col = zeros(ny,nx,Np);

for ip = 1:Np

    omg_col(:,:,ip) = g{ip};

    for ie = 1:Ne
        omg_col(:,:,ip) = omg_col(:,:,ip) - e_col{ip}{ie} .* mu_e{ie};
    end

end

%Copy collapsed result back to grain-resolved variables
Ngrain = numel(c_full);
c_out  = cell(1,Ngrain);
omg    = zeros(ny,nx,Ngrain);

for ig = 1:Ngrain

    iph = grain_to_phase(ig);

    c_out{ig}    = c_col{iph};
    omg(:,:,ig)  = omg_col(:,:,iph);

end

e = Calc_e(pars_full,c_out);

%Store collapsed active set information
LE_state.phase_index    = phase_index;
LE_state.grain_to_phase = grain_to_phase;

%Write back to STATE
STATE.c        = c_out;
STATE.e        = e;
STATE.mu_e     = mu_e;
STATE.chi      = chi;
STATE.omg      = omg;
STATE.LE_state = LE_state;

end


function p_th = Calc_Thermo_p(p,p_tail,p_full)
%CALC_THERMO_P Calculate thermodynamic phase weights from geometric p.
%
% Raw p is the geometric phase fraction.
% p_th is the thermodynamic phase weight.
%
% A tiny phase-field tail has small or zero thermodynamic weight.
% This avoids forcing local equilibrium on infinitesimal phase tails.

if p_full <= p_tail
    error('p_full must be larger than p_tail')
end

Np = size(p,3);

x = (p - p_tail) ./ max(p_full - p_tail,eps);
x = min(max(x,0),1);

%Smoothstep activation
a = x.^2 .* (3 - 2*x);

%Unnormalized thermodynamic weight
w = p .* a;

wsum = sum(w,3);

p_th = zeros(size(p));

%Normal case
good = wsum > eps;

for ip = 1:Np

    tmp = zeros(size(p(:,:,ip)));
    wi = w(:,:,ip);

    tmp(good) = wi(good) ./ wsum(good);

    p_th(:,:,ip) = tmp;

end

%Fallback: if all weights are zero, assign the dominant geometric phase
bad = ~good;

if any(bad(:))

    [~,idmax] = max(p,[],3);

    for ip = 1:Np

        tmp = p_th(:,:,ip);
        tmp(bad & idmax == ip) = 1;
        p_th(:,:,ip) = tmp;

    end

end

end


function [pars_c,c_c,p_c,grain_to_phase] = Collapse_LE_Phases(pars,c,p,phase_index)

phase_index = phase_index(:).';
phase_id    = unique(phase_index,'stable');

Ngrain      = numel(c);
Nphase      = numel(phase_id);
N           = size(p,2);

grain_to_phase = zeros(1,Ngrain);

for iph = 1:Nphase
    grain_to_phase(phase_index == phase_id(iph)) = iph;
end

%Representative thermodynamic data
pars_c = cell(1,Nphase);

for iph = 1:Nphase
    ig = find(grain_to_phase == iph,1,'first');
    pars_c{iph} = pars{ig};
end

%Collapse p by summing grains of the same phase
p_c = zeros(1,N,Nphase);

for ig = 1:Ngrain
    iph = grain_to_phase(ig);
    p_c(:,:,iph) = p_c(:,:,iph) + p(:,:,ig);
end

%Collapse c by p-weighted average
c_c = cell(1,Nphase);

for iph = 1:Nphase

    grains = find(grain_to_phase == iph);
    ig0    = grains(1);
    Nc     = numel(c{ig0});

    c_c{iph} = cell(1,Nc);

    den = zeros(size(c{ig0}{1}));

    for ig = grains
        den = den + reshape(p(:,:,ig),size(den));
    end

    for ic = 1:Nc

        num = zeros(size(c{ig0}{ic}));

        for ig = grains
            w   = reshape(p(:,:,ig),size(num));
            num = num + w .* c{ig}{ic};
        end

        tmp  = c{ig0}{ic};
        good = den > eps;

        tmp(good)  = num(good) ./ den(good);
        tmp(~good) = c{ig0}{ic}(~good);

        c_c{iph}{ic} = tmp;

    end

end

end