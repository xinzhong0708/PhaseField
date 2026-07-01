clear; clc

script_dir = fileparts(mfilename('fullpath'));
repo_root  = fileparts(script_dir);
cd(repo_root)

addpath(fullfile(repo_root,'bin'),'-begin')
addpath(fullfile(repo_root,'Thermo'),'-begin')
addpath(fullfile(repo_root,'Thermo','Solutions'),'-begin')
addpath(fullfile(repo_root,'Maps'),'-begin')

state_file = fullfile(repo_root,'Maps','Map2d.mat');
metadata_file = fullfile(repo_root,'Metadata.xlsx');

S = load(state_file);
STATE = S.STATE;
MODEL = S.MODEL;
GRID  = S.GRID;
PARAM = S.PARAM;

[PHYS,NUM,PARAM,MODEL,GRID] = Read_PFM_Metadata(metadata_file,GRID,MODEL,STATE,PARAM.eta,PARAM);
[Tcur,Pcur] = PT_Path(NUM.t_phy,PARAM.PT);
PARAM.PT_update_count = 0;
[MODEL,PARAM] = Update_Model_PT(MODEL,PARAM,PHYS,Tcur,Pcur);

Ne = numel(STATE.E);
Ebulk = zeros(1,Ne);
for ie = 1:Ne
    Ebulk(ie) = mean(STATE.E{ie},'all');
end

fprintf('\n=== Equilibrium consistency diagnostic ===\n')
fprintf('State file : %s\n',state_file)
fprintf('P-T        : T = %.12g K, P = %.12g Pa\n',Tcur,Pcur)
fprintf('Phases     :')
fprintf(' %s',MODEL.phs_name{:})
fprintf('\n')
fprintf('PhaseG     : %s\n',which('PhaseG'))
fprintf('PhaseThermo: %s\n',which('PhaseThermo'))
fprintf('Ebulk      :')
fprintf(' %.12g',Ebulk)
fprintf('\n')

lp = PhaseDiagramLP(Tcur,Pcur,Ebulk,MODEL);
fprintf('\nPhase_Diagram LP phase proportions:\n')
PrintPhaseVector(MODEL.phs_name,lp.phase_prop)

phase_current = PhaseFractions(STATE.p,MODEL.phase_index,numel(MODEL.phs_name));
fprintf('\nCurrent PF phase proportions:\n')
PrintPhaseVector(MODEL.phs_name,phase_current)

pars_phase = PARAM.pars_phase_PT;
c0 = CollapseCByPhase(STATE.c,MODEL.phase_index,numel(MODEL.phs_name));
Ecell = num2cell(Ebulk);

eta0 = PHYS.eta;
if ~isscalar(eta0)
    eta0 = mean(eta0(:));
end

starts = [ ...
    max(lp.phase_prop,1e-8); ...
    max(phase_current,1e-8); ...
    ones(1,numel(MODEL.phs_name))/numel(MODEL.phs_name)];

best = struct('f',inf,'p',[],'omega',[],'mu',[],'c',[]);

opts = optimoptions('fmincon', ...
    'Display','none', ...
    'Algorithm','sqp', ...
    'MaxFunctionEvaluations',2000, ...
    'MaxIterations',300, ...
    'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-9);

Np = numel(MODEL.phs_name);
lb = 1e-8*ones(1,Np);
ub = ones(1,Np);
Aeq = ones(1,Np);
beq = 1;

for is = 1:size(starts,1)
    p0 = starts(is,:);
    p0 = p0/sum(p0);

    fun = @(pvec) KKSObjective(pvec,pars_phase,c0,Ecell,eta0);
    [pfit,ffit,exitflag] = fmincon(fun,p0,[],[],Aeq,beq,lb,ub,[],opts);

    [fval,detail] = KKSObjective(pfit,pars_phase,c0,Ecell,eta0);

    fprintf('\nKKS fmincon start %d: exitflag = %d, F = %.12g\n',is,exitflag,fval)
    PrintPhaseVector(MODEL.phs_name,pfit)
    fprintf('omega spread = %.12g\n',max(detail.omega)-min(detail.omega))

    if ffit < best.f
        best.f = ffit;
        best.p = pfit;
        best.omega = detail.omega;
        best.mu = detail.mu;
        best.c = detail.c;
    end
end

[~,lp_detail] = KKSObjective(max(lp.phase_prop,1e-8),pars_phase,c0,Ecell,eta0);
[~,cur_detail] = KKSObjective(max(phase_current,1e-8),pars_phase,c0,Ecell,eta0);

fprintf('\nBest KKS phase proportions:\n')
PrintPhaseVector(MODEL.phs_name,best.p)

fprintf('\nOmega at LP proportions:\n')
PrintOmega(MODEL.phs_name,lp_detail.omega)
fprintf('LP omega spread = %.12g\n',max(lp_detail.omega)-min(lp_detail.omega))

fprintf('\nOmega at current PF proportions:\n')
PrintOmega(MODEL.phs_name,cur_detail.omega)
fprintf('Current omega spread = %.12g\n',max(cur_detail.omega)-min(cur_detail.omega))

fprintf('\nOmega at best KKS proportions:\n')
PrintOmega(MODEL.phs_name,best.omega)
fprintf('Best omega spread = %.12g\n',max(best.omega)-min(best.omega))
fprintf('Best finite-eta composition residual mu/eta:\n')
for ie = 1:numel(best.mu)
    fprintf('  E%d: %.12g\n',ie,best.mu{ie}/eta0)
end


function out = PhaseDiagramLP(T,P,Ebulk,MODEL)

td = init_thermo(MODEL.phs_name,MODEL.Cname,MODEL.solmod);
p = props_generate(td);
[g0,v0] = tl_g0(T,P,td);
[g,Npc,pc_id] = tl_gibbs_energy(T,P,MODEL.phs_name,td,p,g0,v0);

g = g./sum(Npc(1:end-1,:))';
Npc = Npc(1:end-1,:)./sum(Npc(1:end-1,:));

Nsys = [Ebulk,1-sum(Ebulk)];

LB = zeros(1,length(g));
UB = ones(1,length(g));
opts = optimoptions('linprog','Display','none');
alph = linprog(g,[],[],Npc,Nsys,LB,UB,opts);
alph = alph/sum(alph);

Nphase = numel(MODEL.phs_name);
phase_prop = zeros(1,Nphase);
for ip = 1:Nphase
    phase_prop(ip) = sum(alph(pc_id == ip));
end

out.phase_prop = phase_prop;
out.alpha = alph;
out.pc_id = pc_id;
out.Npc = Npc;
out.g = g;

end


function [F,detail] = KKSObjective(pvec,pars_phase,c0,Ecell,eta0)

pvec = pvec(:).';
pvec = max(pvec,1e-12);
pvec = pvec/sum(pvec);

p = reshape(pvec,1,1,[]);

try
    [c,mu,chi] = LE_Calculator(pars_phase,p,c0,Ecell,eta0,[0.25,500,1e-10,12,1e-9]);
    e = Calc_e(pars_phase,c);

    Emix = zeros(1,numel(Ecell));
    for ip = 1:numel(pars_phase)
        for ie = 1:numel(Ecell)
            Emix(ie) = Emix(ie) + pvec(ip)*e{ip}{ie};
        end
    end

    Gmix = 0;
    omega = zeros(1,numel(pars_phase));
    for ip = 1:numel(pars_phase)
        gi = PhaseG(pars_phase{ip},c{ip});
        Gmix = Gmix + pvec(ip)*gi;
        omega(ip) = gi;
        for ie = 1:numel(Ecell)
            omega(ip) = omega(ip) - e{ip}{ie}*mu{ie};
        end
    end

    res = cell2mat(Ecell(:)).' - Emix;
    F = Gmix + 0.5*eta0*sum(res.^2);

    if ~isfinite(F)
        F = 1e30;
    end

catch ME
    warning('KKSObjective failed: %s',ME.message)
    F = 1e30;
    c = c0;
    mu = cell(size(Ecell));
    chi = [];
    omega = nan(1,numel(pars_phase));
end

if nargout > 1
    detail.c = c;
    detail.mu = mu;
    detail.chi = chi;
    detail.omega = omega;
end

end


function c_phase = CollapseCByPhase(c_grain,phase_index,Nphase)

c_phase = cell(1,Nphase);

for ip = 1:Nphase
    ig = find(phase_index == ip,1,'first');
    if isempty(ig)
        error('CollapseCByPhase: phase %d has no grain.',ip)
    end
    c_phase{ip} = c_grain{ig};
    for ic = 1:numel(c_phase{ip})
        c_phase{ip}{ic} = mean(c_phase{ip}{ic},'all');
    end
end

end


function phase_prop = PhaseFractions(p,phase_index,Nphase)

phase_prop = zeros(1,Nphase);
for ip = 1:Nphase
    grains = find(phase_index == ip);
    phase_prop(ip) = mean(sum(p(:,:,grains),3),'all');
end

end


function PrintPhaseVector(phs_name,p)

for ip = 1:numel(phs_name)
    fprintf('  %-16s %.12f\n',phs_name{ip},p(ip))
end

end


function PrintOmega(phs_name,omega)

for ip = 1:numel(phs_name)
    fprintf('  %-16s %.12g\n',phs_name{ip},omega(ip))
end

end
