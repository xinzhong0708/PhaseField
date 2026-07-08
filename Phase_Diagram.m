clear
addpath ./ ./bin ./Thermo ./Thermo/Utilities ./Thermo/Solutions ./Thermo/EOS ./Thermo/glpkmex

% Use the same initial P-T source as Run_2D_Scaled.m.  The PF code reads
% Metadata.xlsx, so keeping this enabled prevents a silent phase-diagram/PF
% offset when the workbook is changed.
metadata_file = 'Metadata.xlsx';
sync_PT_from_metadata = true;

% load 1500
% phase         = [interp1(TIME,PHASE(:,1),NUM.t_phy) interp1(TIME,PHASE(:,2),NUM.t_phy) interp1(TIME,PHASE(:,3),NUM.t_phy) interp1(TIME,PHASE(:,4),NUM.t_phy)];
% [Tcur,Pcur]   =  PT_Path(NUM.t_phy,PARAM.PT);
T             =  673;
P             =  3e8;
solmod        = 'solution_models_PFM';

Cname         = {'Fe'      'Mg'       'Ca'       'Al'     'Si'    'O'};
Nsys          = [ 0.0897    0.0960    0.0574    0.1840];
Nsys          = [Nsys 1-sum(Nsys)];

% Choose phases considered in Gibbs minimization.  Keep the order consistent
% with the PFM map/MODEL.phs_name so printed phase proportions can be
% compared directly with PF phase columns.
phs_name      = {'Garnet','Clinopyroxene','Kyanite','Quartz'};

td            =  init_thermo(phs_name,Cname,solmod);
p             =  props_generate(td);     % generate endmember proportions

% Minimization refinement
[g0,v0]       =  tl_g0(T,P,td);
[g,Npc,pc_id] =  tl_gibbs_energy(T,P,phs_name,td,p,g0,v0);

%Normalize
g             =  g./sum(Npc(1:end-1,:))';
Npc           =  Npc(1:end-1,:)./sum(Npc(1:end-1,:));

%Normalize based on Npc
LB            =  zeros(1,length(g));
UB            =   ones(1,length(g));
[alph,gmin]   =  linprog(g,[],[],Npc,Nsys,LB,UB);



alph          = alph/sum(alph);
id            = find(alph > 1e-6);
disp(Npc(:,id))

disp('Active phase proportion')
disp(alph(id))

disp('Active phase ID')
disp(pc_id(id))

disp(Npc(:,id))

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

