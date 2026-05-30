function dc_all = LE_Mu_EProjected_Step(pars,p,c,mu_e,E_in,eta,lam,wmu)

Np     = length(c);
Ne     = length(mu_e);
N      = numel(mu_e{1});
mu_mat = cell2mat(mu_e(:));
E_mat  = cell2mat(E_in(:));

if isscalar(eta)
    eta_vec = eta*ones(1,N);
else
    eta_vec = eta(:).';
end

% current mixture composition
Emix = zeros(Ne,N);
Rall = cell(1,Np);

for ip = 1:Np
    Rall{ip} = PhaseThermo(pars{ip},c{ip});
    e_ip     = cell2mat(Rall{ip}.e(:));
    p_ip     = reshape(p(:,:,ip),1,N);
    Emix     = Emix + p_ip .* e_ip;
end

R_E_global = E_mat - Emix - mu_mat ./ reshape(eta_vec,1,N);

dc_all = c;
for ip = 1:Np
    for ic = 1:length(c{ip})
        dc_all{ip}{ic} = zeros(size(c{ip}{ic}));
    end
end

for ip = 1:Np

    R = Rall{ip};

    if isempty(R.mu_c) || isempty(R.H_c) || isempty(R.Jac)
        continue
    end

    mu_c = cell2mat(R.mu_c(:));
    H    = 0.5*(R.H_c + permute(R.H_c,[2 1 3]));
    J    = R.Jac;
    Nc   = size(mu_c,1);

    JTmu3 = pagemtimes(permute(J,[2 1 3]),reshape(mu_mat,Ne,1,N));
    R_mu  = reshape(JTmu3,Nc,N) - mu_c;

    p_ip = reshape(p(:,:,ip),1,N);
    R_E  = R_E_global .* p_ip;

    JTJ = pagemtimes(permute(J,[2 1 3]),J);
    HTH = pagemtimes(permute(H,[2 1 3]),H);

    A = JTJ + wmu*HTH + repmat(eye(Nc),1,1,N).*lam;

    JTRE3 = pagemtimes(permute(J,[2 1 3]),reshape(R_E,Ne,1,N));
    HTRM3 = pagemtimes(permute(H,[2 1 3]),reshape(R_mu,Nc,1,N));

    b  = JTRE3 + wmu*HTRM3;
    dc = reshape(pagemldivide(A,b),Nc,N);

    for ic = 1:length(c{ip})
        dc_all{ip}{ic} = dc(ic,:);
    end
end
end