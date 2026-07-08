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
T             =  573;
P             =  3e8;
solmod        = 'solution_models_PFM';

Cname         = {'Fe'      'Mg'       'Ca'       'Al'     'Si'    'O'};
%Molar fraction of elements
Nsys          = [0.0923    0.1054    0.0061    0.2643];
Nsys          = [Nsys 1-sum(Nsys)];

%Consider O
N_O           = [1    1     1   1.5   2];
Nsys          = [Nsys , sum(Nsys.*N_O)];
Nsys          =  Nsys/sum(Nsys);

% Choose phases considered in Gibbs minimization.  Keep the order consistent
% with the PFM map/MODEL.phs_name so printed phase proportions can be
% compared directly with PF phase columns.
phs_name      = {'Orthopyroxene','Garnet','Kyanite','Quartz'};
td            =  init_thermo(phs_name,Cname,solmod);
p             =  props_generate(td);     % generate endmember proportions

% Minimization refinement
[g0,v0]       =  tl_g0(T,P,td);
[g,Npc,pc_id] =  tl_gibbs_energy(T,P,phs_name,td,p,g0,v0);

%Normalize
% g             =  g./sum(Npc(1:end-1,:))';
% Npc           =  Npc(1:end-1,:)./sum(Npc(1:end-1,:));

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


function [T,P] = Read_First_Metadata_PT(metadata_file,T_default,P_default)

T = T_default;
P = P_default;

if exist(metadata_file,'file') ~= 2
    if exist('Metadata.xlsx','file') == 2
        metadata_file = 'Metadata.xlsx';
    elseif exist('MetaData.xlsx','file') == 2
        metadata_file = 'MetaData.xlsx';
    end
end

if exist(metadata_file,'file') ~= 2
    warning('Phase_Diagram: metadata file "%s" not found. Using hardcoded P-T.',metadata_file)
    return
end

try
    Cpt = readcell(metadata_file,'Sheet','PTt_Path');
catch ME
    warning('Phase_Diagram: could not read PTt_Path from "%s": %s. Using hardcoded P-T.', ...
        metadata_file,ME.message)
    return
end

header_row = [];

for i = 1:size(Cpt,1)

    row_txt = cell(1,size(Cpt,2));

    for j = 1:size(Cpt,2)
        row_txt{j} = strtrim(Cell_String_Local(Cpt{i,j}));
    end

    if any(strcmpi(row_txt,'P (Pa)')) && ...
       any(strcmpi(row_txt,'T (K)'))  && ...
       any(strcmpi(row_txt,'t (s)'))
        header_row = i;
        break
    end
end

if isempty(header_row)
    warning('Phase_Diagram: PTt_Path header row not found. Using hardcoded P-T.')
    return
end

colP = [];
colT = [];
colt = [];

for j = 1:size(Cpt,2)

    txt = strtrim(Cell_String_Local(Cpt{header_row,j}));

    if strcmpi(txt,'P (Pa)')
        colP = j;
    elseif strcmpi(txt,'T (K)')
        colT = j;
    elseif strcmpi(txt,'t (s)')
        colt = j;
    end
end

pt_rows = [];

for i = header_row+1:size(Cpt,1)

    if Is_Empty_Cell_Local(Cpt{i,colP}) || ...
       Is_Empty_Cell_Local(Cpt{i,colT}) || ...
       Is_Empty_Cell_Local(Cpt{i,colt})
        continue
    end

    pt_rows(end+1,:) = [ ...
        Cell_Number_Local(Cpt{i,colt}), ...
        Cell_Number_Local(Cpt{i,colT}), ...
        Cell_Number_Local(Cpt{i,colP})]; %#ok<AGROW>
end

if isempty(pt_rows)
    warning('Phase_Diagram: PTt_Path contains no valid rows. Using hardcoded P-T.')
    return
end

[~,ord] = sort(pt_rows(:,1));
pt_rows = pt_rows(ord,:);

T = pt_rows(1,2);
P = pt_rows(1,3);

fprintf('Phase_Diagram P-T from %s: T = %.12g K, P = %.12g Pa\n', ...
    metadata_file,T,P)

end


function val = Cell_Number_Local(x)

if isnumeric(x)
    val = x;
elseif ischar(x) || isstring(x)
    val = str2double(x);
else
    val = NaN;
end

if ~isfinite(val)
    error('Phase_Diagram: PTt_Path contains a nonnumeric P-T-time value.')
end

end


function s = Cell_String_Local(x)

if ismissing(x)
    s = '';
elseif ischar(x)
    s = x;
elseif isstring(x)
    s = char(x);
elseif isnumeric(x)
    s = num2str(x);
else
    s = '';
end

end


function tf = Is_Empty_Cell_Local(x)

tf = false;

if isempty(x)
    tf = true;
elseif ismissing(x)
    tf = true;
elseif ischar(x) || isstring(x)
    tf = strlength(string(x)) == 0;
elseif isnumeric(x)
    tf = isempty(x) || any(isnan(x),'all');
end

end

