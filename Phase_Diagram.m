clear
addpath ./ ./Thermo/Utilities/ ./Thermo/Solutions/ ./EOS ./glpkmex

T        = linspace(1000,1000,1) + 273.15;
P        = linspace(2.5,2.5,1) * 1e9;
solmod   = 'solution_models_PFM';

Cname    = {'Si' 'Fe' 'Mg' 'Ca' 'Al' 'Na' 'O'};

Nsys     = [0  0.110175735098627   0.152682123916501   0.057659754069281   0.189020427085893   0.034425893359538];

Nsys(1)  = 1-sum(Nsys);

Nsys     = [Nsys Nsys(1)*2 + Nsys(2) + Nsys(3) + Nsys(4) + Nsys(5)/2*3];

%Choose phases considered in Gibbs minimization
phs_name = {'Garnet','Feldspar','Clinopyroxene'};

td            = init_thermo(phs_name,Cname,solmod);
p             = props_generate(td);

%Minimization refinement
[g0,v0]       = tl_g0(T,P,td);
[g,Npc,pc_id] = tl_gibbs_energy(T,P,phs_name,td,p,g0,v0);

%Normalize on non-oxygen components
norm_fac      = sum(Npc(1:end-1,:),1);

g             = g ./ norm_fac.';
Npc           = Npc(1:end-1,:) ./ repmat(norm_fac,size(Npc,1)-1,1);

% %----------------------------------------------------------
% % GLPK Gibbs minimization
% %----------------------------------------------------------
% f             = real(g(:));
% Aeq           = Npc;
% beq           = Nsys(1:end-1).';
% 
% nvar          = numel(f);
% 
% LB            = zeros(nvar,1);
% UB            = ones(nvar,1);
% 
% ctype         = repmat('S',size(Aeq,1),1);   %Aeq*alph = beq
% vartype       = repmat('C',nvar,1);          %Continuous variables
% sense         = 1;                           %Minimization
% 
% param         = struct();
% param.msglev  = 1;
% param.presol  = 1;
% param.lpsolver = 1;                          %Important: revised simplex
% param.round   = 0;                           %Keep raw values for checking

% [alph,gmin,status,extra] = glpk(f,Aeq,beq,LB,UB,ctype);

%----------------------------------------------------------
% LINPROG Gibbs minimization
%----------------------------------------------------------
f             = real(g(:));
nvar          = numel(f);

LB            = zeros(nvar,1);
UB            = ones(nvar,1);

[alph,gmin] = linprog(f,[],[],Npc,Nsys(1:end-1).',LB,UB);

alph        = alph/sum(alph);
id          = find(alph > 1e-5);

disp('Active phase proportion')
disp(alph(id))

disp('Active phase ID')
disp(pc_id(id))

disp('Resulting phase proportion and endmember proportion')

for iph = 1:length(phs_name)

    idp      = find(pc_id == iph);
    ph_prop  = sum(alph(idp));

    if ph_prop > 1e-5

        ph_comp = alph(idp).' * p{iph} / ph_prop;

        fprintf('\n%s\n',phs_name{iph})
        fprintf('Phase proportion:\n')
        disp(ph_prop)

        fprintf('Endmember proportion:\n')
        disp(ph_comp)

    end

end