%% Diagnose_PT_Zoning_Workflow.m
% Inspect the current P-T zoning workflow as MATLAB actually reads it.

clear

repo_root = fileparts(fileparts(mfilename('fullpath')));
if isempty(repo_root)
    repo_root = pwd;
end
cd(repo_root)

addpath('bin')
addpath('Thermo')
addpath('Thermo/Solutions')
addpath('Maps')

metadata_file = 'Metadata.xlsx';
map_file      = 'Map2d.mat';

fprintf('\n=== MATLAB paths ===\n')
fprintf('repo root   : %s\n',repo_root)
fprintf('PhaseG      : %s\n',which('PhaseG'))
fprintf('PhaseThermo : %s\n',which('PhaseThermo'))
fprintf('LE_Run_Mode : %s\n',which('LE_Run_Mode'))

fprintf('\n=== Map file ===\n')
map_resolved = which(map_file);
if isempty(map_resolved)
    error('Map file "%s" was not found.',map_file)
end

S = load(map_file);
fprintf('map file    : %s\n',map_resolved)
fprintf('MODEL.T/P   : %.12g K, %.12g Pa\n',S.MODEL.T,S.MODEL.P)
fprintf('phases      : %s\n',strjoin(S.MODEL.phs_name,', '))
fprintf('phase_index :')
fprintf(' %d',S.MODEL.phase_index)
fprintf('\n')

p_phase = collapse_p_by_phase(S.STATE.p,S.MODEL.phase_index,numel(S.MODEL.phs_name));
fprintf('phase fractions from STATE.p:\n')
for ip = 1:numel(S.MODEL.phs_name)
    fprintf('  %-18s %.12g\n',S.MODEL.phs_name{ip},mean(p_phase(:,:,ip),'all'))
end

fprintf('mean bulk E:\n')
for ie = 1:numel(S.STATE.E)
    fprintf('  E%-2d %.12g\n',ie,mean(S.STATE.E{ie},'all'))
end

fprintf('\n=== Metadata as read by Read_PFM_Metadata ===\n')
[PHYS,NUM,PARAM,MODEL,GRID] = Read_PFM_Metadata( ...
    metadata_file,S.GRID,S.MODEL,S.STATE,S.PARAM.eta,S.PARAM);

fprintf('initial dt / max dt / min dt code: %.12g %.12g %.12g\n', ...
    NUM.dt_phy,NUM.dt_max,NUM.dt_min)
fprintf('time scale t_sc                   : %.12g s\n',PHYS.t_sc)
fprintf('PT path code time range            : %.12g to %.12g\n', ...
    min(PARAM.PT.t_path),max(PARAM.PT.t_path))
fprintf('PT path SI time range              : %.12g to %.12g s\n', ...
    min(PARAM.PT.t_path)*PHYS.t_sc,max(PARAM.PT.t_path)*PHYS.t_sc)
fprintf('PT path first/last T               : %.12g to %.12g K\n', ...
    PARAM.PT.T_path(1),PARAM.PT.T_path(end))
fprintf('PT path first/last P               : %.12g to %.12g Pa\n', ...
    PARAM.PT.P_path(1),PARAM.PT.P_path(end))
fprintf('Update PT every                    : %d accepted calls\n',PARAM.update_PT_every)
fprintf('use_kappa_c / CS cap               : %d / %.12g\n', ...
    PARAM.use_kappa_c,PARAM.CS_chi_cap)
fprintf('M_L_floor_fac / interface_only     : %.12g / %d\n', ...
    PARAM.M_L_floor_fac,PARAM.M_L_interface_only)

fprintf('\nPHYS.M_phs_raw, rows follow MODEL.phs_name, columns MODEL.Cname without dependent species:\n')
fprintf('%-18s', 'phase')
for ie = 1:numel(MODEL.Cname)-2
    fprintf(' %14s',MODEL.Cname{ie})
end
fprintf('\n')
for ip = 1:numel(MODEL.phs_name)
    fprintf('%-18s',MODEL.phs_name{ip})
    for ie = 1:size(PHYS.M_phs_raw,2)
        fprintf(' %14.6e',PHYS.M_phs_raw(ip,ie))
    end
    fprintf('\n')
end

fprintf('\nPHYS.M_phs scaled:\n')
fprintf('%-18s', 'phase')
for ie = 1:numel(MODEL.Cname)-2
    fprintf(' %14s',MODEL.Cname{ie})
end
fprintf('\n')
for ip = 1:numel(MODEL.phs_name)
    fprintf('%-18s',MODEL.phs_name{ip})
    for ie = 1:size(PHYS.M_phs,2)
        fprintf(' %14.6e',PHYS.M_phs(ip,ie))
    end
    fprintf('\n')
end

fprintf('\n=== P-T update sampling ===\n')
sample_t = linspace(min(PARAM.PT.t_path),max(PARAM.PT.t_path),5);
for k = 1:numel(sample_t)
    [Tcur,Pcur] = PT_Path(sample_t(k),PARAM.PT);
    fprintf('  t_code %.12g  t_SI %.12g s  T %.12g K  P %.12g Pa\n', ...
        sample_t(k),sample_t(k)*PHYS.t_sc,Tcur,Pcur)
end

fprintf('\n=== Mobility field from initial map ===\n')
PARAM = Compute_M_And_L(S.STATE,PARAM,MODEL,PHYS);
for ie = 1:numel(S.STATE.E)
    Mie = PARAM.M{ie,ie};
    fprintf('  M%d min/median/max = %.6e %.6e %.6e\n', ...
        ie,min(Mie(:)),median(Mie(:)),max(Mie(:)))
end
fprintf('  L min/median/max  = %.6e %.6e %.6e\n', ...
    min(PARAM.L(:)),median(PARAM.L(:)),max(PARAM.L(:)))

fprintf('\n=== First LE at metadata first P-T ===\n')
[MODEL,PARAM] = Update_Model_PT(MODEL,PARAM,PHYS,PARAM.PT.T_path(1),PARAM.PT.P_path(1));
STATE_LE = LE_Run_Mode(S.STATE,PARAM,MODEL);
E_mix = Calc_E_Tot(STATE_LE.e,STATE_LE.p);
for ie = 1:numel(STATE_LE.E)
    dmu = max(STATE_LE.mu_e{ie}(:)) - min(STATE_LE.mu_e{ie}(:));
    res = max(abs(PARAM.eta(:).*(STATE_LE.E{ie}(:)-E_mix{ie}(:)) - STATE_LE.mu_e{ie}(:)));
    fprintf('  elem %d: mu jump %.6e, eta*(E-Emix)-mu max %.6e\n',ie,dmu,res)
end

fprintf('\n=== Newest numeric checkpoint compared with current metadata ===\n')
ckpt = newest_numeric_checkpoint(repo_root);
if isempty(ckpt)
    fprintf('  no numeric checkpoint files were found in the repo root\n')
else
    C = load(ckpt);
    fprintf('  checkpoint : %s\n',ckpt)
    if isfield(C,'MODEL') && isfield(C.MODEL,'T') && isfield(C.MODEL,'P')
        fprintf('  MODEL.T/P  : %.12g K, %.12g Pa\n',C.MODEL.T,C.MODEL.P)
    end
    if isfield(C,'NUM') && isfield(C.NUM,'t_phy')
        fprintf('  NUM.t_phy  : %.12g\n',C.NUM.t_phy)
    end
    if isfield(C,'PHYS') && isfield(C.PHYS,'M_phs_raw')
        fprintf('  checkpoint PHYS.M_phs_raw:\n')
        disp(C.PHYS.M_phs_raw)
        fprintf('  current metadata PHYS.M_phs_raw:\n')
        disp(PHYS.M_phs_raw)
    end
    if isfield(C,'PARAM') && isfield(C.PARAM,'PT')
        fprintf('  checkpoint PT path:\n')
        for k = 1:numel(C.PARAM.PT.t_path)
            fprintf('    %3d  t %.12g  T %.12g  P %.12g\n', ...
                k,C.PARAM.PT.t_path(k),C.PARAM.PT.T_path(k),C.PARAM.PT.P_path(k))
        end
    end
end

function p_phase = collapse_p_by_phase(p,phase_index,Nphase)

[ny,nx,~] = size(p);
p_phase = zeros(ny,nx,Nphase);

for ip = 1:Nphase
    grains = find(phase_index == ip);
    if ~isempty(grains)
        p_phase(:,:,ip) = sum(p(:,:,grains),3);
    end
end

end

function fname = newest_numeric_checkpoint(repo_root)

D = dir(fullfile(repo_root,'*.mat'));
best_num = -inf;
fname = '';

for i = 1:numel(D)
    [~,base] = fileparts(D(i).name);
    val = str2double(base);
    if ~isnan(val) && val > best_num
        best_num = val;
        fname = fullfile(D(i).folder,D(i).name);
    end
end

end
