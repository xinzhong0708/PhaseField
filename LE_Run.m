function [c,mu_e,chi,omg,LE_state] = LE_Run(pars,p,c,E,mu_e,chi,eta_vec,level1,level2,LE_state)

if nargin < 10 || isempty(LE_state)
    LE_state = struct();
end
%Reshape 2D into 1D
nx         = size(p,2);
ny         = size(p,1);
c          = Unpack_c(c);
E          = Unpack_E(E);
p          = Unpack_p(p);
mu_e       = Unpack_E(mu_e);
chi        = Unpack_Chi(chi);
eta_vec    = eta_vec(:);

Pmax       = 4;
pmin       = 3e-2;

%Hysteresis thresholds
p_on       = pmin;
p_off      = 0.5*pmin;

Np         = numel(c);
N          = numel(c{1}{1});

%Build p matrix
pMat       = zeros(N,Np);
for ip = 1:Np
    p_ip          = squeeze(p(:,:,ip));
    pMat(:,ip)   = p_ip(:);
end

%Initialize hysteresis active set
if ~isfield(LE_state,'active') || isempty(LE_state.active) || ~isequal(size(LE_state.active),[N,Np])
    LE_state.active = pMat > p_on;
    % Make sure every grid point has at least one active phase
    for i = 1:N
        if ~any(LE_state.active(i,:))
            [~,idmax]              = max(pMat(i,:));
            LE_state.active(i,idmax) = true;
        end
    end
end

% Update hysteresis active set
active_old = LE_state.active;
active_new = false(N,Np);
for i = 1:N
    for ip = 1:Np
        if active_old(i,ip)
            % Active phase remains active until p drops below p_off
            active_new(i,ip) = pMat(i,ip) > p_off;
        else
            % Inactive phase becomes active only when p rises above p_on
            active_new(i,ip) = pMat(i,ip) > p_on;
        end
    end
    %Make sure every grid point has at least one active phase
    if ~any(active_new(i,:))
        [~,idmax]            = max(pMat(i,:));
        active_new(i,idmax)  = true;
    end

    %Limit maximum number of active phases
    if sum(active_new(i,:)) > Pmax
        score                = pMat(i,:);
        %Favor phases that were already active
        score(active_old(i,:)) = score(active_old(i,:)) + 0.5*p_on;
        [~,ord]              = sort(score,'descend');
        tmp                  = false(1,Np);
        tmp(ord(1:Pmax))     = true;
        active_new(i,:)      = tmp;
    end
end
LE_state.active = active_new;

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
        if ~any(mask);  continue;  end

        %Slice local fields
        c_cur = Slice_c(c, ph_act, mask);
        E_cur = Slice_E(E, mask);
        eta   = eta_vec(mask);

        %One phase
        if k == 1
            ip                     = ph_act(1);
            [c_tmp,mu_tmp,chi_tmp] = LE_Calculator(pars(ip), ones(size(c_cur{1}{1})), c_cur(1), E_cur, eta, level2);
            [c, mu_e, chi]         = Assign_LE_Back(c, mu_e, chi,  c_tmp, mu_tmp, chi_tmp, ip, mask);
        end

        %Two phase
        if k == 2
            p_cur                  = Norm_Phi(p(:,mask,ph_act));
            [c_tmp,mu_tmp,chi_tmp] = LE_Calculator(pars(ph_act), p_cur, c_cur, E_cur, eta, [0.3,1000]);
            [c, mu_e, chi]         = Assign_LE_Back(c, mu_e, chi, c_tmp, mu_tmp, chi_tmp, ph_act, mask);
        end

        %Three phase
        if k == 3
            p_cur                  = Norm_Phi(p(:,mask,ph_act));
            [c_tmp,mu_tmp,chi_tmp] = LE_Calculator(pars(ph_act), p_cur, c_cur, E_cur, eta, [0.2,1000]);
            [c, mu_e, chi]         = Assign_LE_Back(c, mu_e, chi, c_tmp, mu_tmp, chi_tmp, ph_act, mask);
        end

        %Four phase
        if k == 4
            p_cur                  = Norm_Phi(p(:,mask,ph_act));
            [c_tmp,mu_tmp,chi_tmp] = LE_Calculator(pars(ph_act), p_cur, c_cur, E_cur, eta, [0.1,1000]);
            [c, mu_e, chi]         = Assign_LE_Back(c, mu_e, chi, c_tmp, mu_tmp, chi_tmp, ph_act, mask);
        end
    end
end

%Pack up
g = cell(1,Np);
for ip = 1:Np
    g{ip} = reshape(PhaseG(pars{ip}, c{ip}), ny, []);
end

c     = Pack_c(c, ny);
mu_e  = Pack_E(mu_e, ny);
chi   = Pack_chi(chi, ny);
e     = Calc_e(pars, c);

%Calculate omega
Ne    = length(E);
omg   = ones(ny, nx, Np);
for ip = 1:Np
    omg(:,:,ip) = g{ip};
    for ie = 1:Ne
        omg(:,:,ip) = omg(:,:,ip) - e{ip}{ie} .* mu_e{ie};
    end
end

end

