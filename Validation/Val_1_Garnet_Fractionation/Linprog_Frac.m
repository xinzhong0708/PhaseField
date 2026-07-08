clear
addpath ./ ../.././bin ../../Thermo ../../Thermo/Utilities ../../Thermo/Solutions ../../Thermo/EOS ../../Thermo/glpkmex

% -------------------------------------------------------------------------
% PTt path
% Columns: P, T, t
% -------------------------------------------------------------------------
metadata_file = 'Metadata.xlsx';
sync_PT_from_metadata = true;
 
PTt           = [3.00E+08    773      0
                 3.00E+08    773      1.00E+08
                 1.00E+09    1073     2.00E+09
                 1.00E+09    1073     1.00E+13];

if sync_PT_from_metadata && exist(metadata_file,'file')
    try
        PTt0 = readmatrix(metadata_file,'Sheet','PTt');
        PTt0 = PTt0(all(isfinite(PTt0),2),:);
        if size(PTt0,2) >= 3
            PTt = PTt0(:,1:3);
        end
    catch
        % Keep hard-coded PTt if Metadata.xlsx does not contain PTt sheet
    end
end

% Time history
t_hist        = linspace(0,2e9,20);

% -------------------------------------------------------------------------
% Thermodynamic setup
% -------------------------------------------------------------------------
solmod        = 'solution_models_PFM';
Cname         = {'Fe','Mg','Ca','Al','Si','O'};

% Independent cation bulk composition.
% Order: Fe Mg Ca Al Si
Nsys0         = [0.083592923365060   0.159527052945144   0.006003759585895   0.171752528207798];
Nsys0         = [Nsys0 1-sum(Nsys0)];
Nsys0         = Nsys0(:);
Nsys0         = Nsys0./sum(Nsys0);

phs_name      = {'Orthopyroxene','Garnet','Kyanite','Quartz'};
igrt          = find(strcmpi(phs_name,'Garnet'));

if isempty(igrt)
    error('Garnet was not found in phs_name.')
end

% Fractionation control
remove_grt_frac = 0.99;     % 1 = remove all equilibrium garnet at each step

% Current remaining absolute cation bulk
Nabs          = Nsys0;
Nabs          = Nabs./sum(Nabs);

% -------------------------------------------------------------------------
% Refinement options
% -------------------------------------------------------------------------
options.nref       = 12;
options.dz         = 1/4;
options.z_window   = 0.085;
options.dz_fact    = 1.5;
options.ref_fact   = 1.25;
options.active_tol = 1e-8;
options.verbose    = 0;
options.solmod     = solmod;
options.Cname      = Cname;

% -------------------------------------------------------------------------
% Storage
% -------------------------------------------------------------------------
nstep         = numel(t_hist);
nph           = numel(phs_name);
nel           = numel(Nsys0);

P_hist        = zeros(nstep,1);
T_hist        = zeros(nstep,1);
gmin_hist     = zeros(nstep,1);

Mremain_hist  = zeros(nstep,1);
Mgrt_hist     = zeros(nstep,1);
Mgrt_cum_hist = zeros(nstep,1);
x_grt_mid     = nan(nstep,1);

Nsys_hist     = zeros(nstep,nel);
phase_hist    = zeros(nstep,nph);

% Garnet composition in independent cation basis: Fe Mg Ca Al Si
grt_elem_hist = nan(nstep,nel);

% Garnet divalent-site composition: Fe Mg Ca normalized by Fe+Mg+Ca
grt_X_hist    = nan(nstep,3);   % columns: XFe, XMg, XCa

alph_hist     = cell(nstep,1);
OUT_hist      = cell(nstep,1);

Mgrt_cum      = 0;

% Element indices in cation basis
Ccat = Cname(1:nel);

iFe = find(strcmpi(Ccat,'Fe'));
iMg = find(strcmpi(Ccat,'Mg'));
iCa = find(strcmpi(Ccat,'Ca'));

% -------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------
for it = 1:nstep

    % Current P and T
    P             = interp1(PTt(:,3),PTt(:,1),t_hist(it));
    T             = interp1(PTt(:,3),PTt(:,2),t_hist(it));

    P_hist(it)    = P;
    T_hist(it)    = T;

    % Current normalized remaining cation bulk
    Nsys               = Nabs./sum(Nabs);
    Nsys_hist(it,:)    = Nsys(:).';

    % ---------------------------------------------------------------------
    % Refined linprog equilibrium on cation basis
    % ---------------------------------------------------------------------
    OUT                = Equilibrium_Linprog_Refine(T,P,phs_name,Nsys,options);
    OUT_hist{it}       = OUT;

    alph               = OUT.alph(:);
    Npc                = OUT.Npc;
    pc_id              = OUT.pc_id(:);
    gmin               = OUT.gmin;

    alph_hist{it}      = alph;
    gmin_hist(it)      = gmin;

    % Phase proportions
    phase_hist(it,:)   = OUT.phase_prop(:).';

    % ---------------------------------------------------------------------
    % Garnet composition before removal
    % ---------------------------------------------------------------------
    id_grt = pc_id == igrt;

    if any(id_grt) && sum(alph(id_grt)) > 0

        % Garnet cation amount in the current normalized rock
        Ngrt = Npc(:,id_grt)*alph(id_grt);

        % Average garnet cation composition: Fe Mg Ca Al Si
        grt_elem = Ngrt./sum(Ngrt);
        grt_elem_hist(it,:) = grt_elem(:).';

        % Fe-Mg-Ca divalent-site composition
        div_sum = grt_elem(iFe) + grt_elem(iMg) + grt_elem(iCa);

        if div_sum > 0
            grt_X_hist(it,1) = grt_elem(iFe)./div_sum;
            grt_X_hist(it,2) = grt_elem(iMg)./div_sum;
            grt_X_hist(it,3) = grt_elem(iCa)./div_sum;
        end

    else

        Ngrt = zeros(nel,1);

    end

    % ---------------------------------------------------------------------
    % Fractionation: remove current equilibrium garnet from remaining bulk
    % ---------------------------------------------------------------------
    Mgrt_before = Mgrt_cum;

    if any(id_grt) && sum(alph(id_grt)) > 0

        % Absolute amount removed from current remaining cation bulk
        dNgrt = remove_grt_frac*sum(Nabs)*Ngrt;

        % Avoid tiny negative values from numerical roundoff
        dNgrt = min(dNgrt,Nabs);

        Nabs  = Nabs - dNgrt;
        Nabs(abs(Nabs) < 1e-14) = 0;

        Mgrt_step = sum(dNgrt);

    else

        Mgrt_step = 0;

    end

    Mgrt_hist(it)     = Mgrt_step;
    Mgrt_cum          = Mgrt_cum + Mgrt_step;
    Mgrt_cum_hist(it) = Mgrt_cum;
    Mremain_hist(it)  = sum(Nabs);

    % 1D linear garnet coordinate.
    % No cubic root is used here.
    if Mgrt_step > 0
        x_grt_mid(it) = Mgrt_before + 0.5*Mgrt_step;
    end

    if sum(Nabs) <= 0
        error('All material has been fractionated by step %d.',it)
    end

    fprintf('step %3d / %3d: T = %.1f K, P = %.3e Pa, grt = %.5f, removed = %.5e, remaining = %.5e\n', ...
        it,nstep,T,P,phase_hist(it,igrt),Mgrt_step,sum(Nabs))

end

% -------------------------------------------------------------------------
% Normalize 1D garnet coordinate
% -------------------------------------------------------------------------
good_grt = isfinite(x_grt_mid) & Mgrt_hist > 0;

if Mgrt_cum > 0
    x_grt_1D = x_grt_mid./Mgrt_cum;
else
    x_grt_1D = x_grt_mid;
end

% -------------------------------------------------------------------------
% Print final result
% -------------------------------------------------------------------------
fprintf('\nFinal remaining normalized bulk:\n')
for ie = 1:nel
    fprintf('  %-2s = %.8e\n',Ccat{ie},Nabs(ie)/sum(Nabs))
end

fprintf('\nFinal remaining mass = %.8e\n',sum(Nabs))
fprintf('Cumulative removed garnet mass = %.8e\n',Mgrt_cum)

fprintf('\nFinal phase proportions before last fractionation:\n')
for iph = 1:nph
    fprintf('  %-15s %.8f\n',phs_name{iph},phase_hist(end,iph))
end

% -------------------------------------------------------------------------
% Plot phase proportions
% -------------------------------------------------------------------------
figure(1); clf
plot(t_hist,phase_hist,'LineWidth',1.5)
xlabel('time')
ylabel('phase proportion')
legend(phs_name,'Location','best')
box on

% -------------------------------------------------------------------------
% Plot cumulative removed garnet
% -------------------------------------------------------------------------
figure(2); clf
plot(t_hist,Mgrt_cum_hist,'LineWidth',1.5)
xlabel('time')
ylabel('cumulative removed garnet')
box on

% -------------------------------------------------------------------------
% Plot remaining bulk composition
% -------------------------------------------------------------------------
figure(3); clf
plot(t_hist,Nsys_hist,'LineWidth',1.5)
xlabel('time')
ylabel('remaining normalized cation bulk composition')
legend(Ccat,'Location','best')
box on

% -------------------------------------------------------------------------
% Plot garnet Fe-Mg-Ca composition profile
% x = 0 is core, x = 1 is rim.
% This assumes 1D linear garnet growth, so no cubic root is used.
% -------------------------------------------------------------------------
figure(4); clf

id_first = find(good_grt,1,'first');

if ~isempty(id_first)

    xplot  = [0; x_grt_1D(good_grt)];
    XFe    = [grt_X_hist(id_first,1); grt_X_hist(good_grt,1)];
    XMg    = [grt_X_hist(id_first,2); grt_X_hist(good_grt,2)];
    XCa    = [grt_X_hist(id_first,3); grt_X_hist(good_grt,3)];

    plot(xplot,XFe,'-','LineWidth',1.5)
    hold on
    plot(xplot,XMg,'-','LineWidth',1.5)
    plot(xplot,XCa,'-','LineWidth',1.5)

end

xlabel('normalized 1D garnet position, core to rim')
ylabel('garnet divalent-site fraction')
legend({'XFe','XMg','XCa'},'Location','best')
box on

% -------------------------------------------------------------------------
% Plot full garnet cation composition profile
% -------------------------------------------------------------------------
figure(5); clf

if any(good_grt)
    plot([x_grt_1D(good_grt)],[grt_elem_hist(good_grt,1:3)],'LineWidth',1.5)
    % plot([0 ; x_grt_1D(good_grt)],[grt_elem_hist(good_grt(1),1:3); grt_elem_hist(good_grt,1:3)],'LineWidth',1.5)
end

xlabel('normalized 1D garnet position, core to rim')
ylabel('garnet normalized cation composition')
legend(Ccat,'Location','best')
box on