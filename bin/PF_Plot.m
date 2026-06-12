function PF_Plot(pos,what,STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%PF_PLOT Phase-index aware plotting helper.
%
% Examples:
%   PF_Plot([3,3,1],'E1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%   PF_Plot([3,3,2],'mu_e1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%   PF_Plot([3,3,3],'dt',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%   PF_Plot([3,3,4],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%   PF_Plot([3,3,5],'omg12',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%   PF_Plot([3,3,6],'p2',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%   PF_Plot([3,3,8],'c21',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%   PF_Plot([3,3,9],'Phase%',STATE,GRID,MODEL,TIME,DTPHY,PHASE)

subplot(pos(1),pos(2),pos(3)); cla

what = strtrim(what);
w    = lower(what);

phase_index = MODEL.phase_index(:).';
phase_ids   = unique(phase_index,'stable');
Nphase      = numel(phase_ids);

if isfield(MODEL,'phs_name') && numel(MODEL.phs_name) >= max(phase_ids)
    phs_name = MODEL.phs_name(phase_ids);
elseif isfield(MODEL,'phs_name') && numel(MODEL.phs_name) == Nphase
    phs_name = MODEL.phs_name;
else
    phs_name = arrayfun(@(i) sprintf('Phase%d',i),phase_ids,'UniformOutput',false);
end

% Collapse grain p to thermodynamic phase p
p_phase = zeros(size(STATE.p,1),size(STATE.p,2),Nphase);

for iph = 1:Nphase
    grains = find(phase_index == phase_ids(iph));
    p_phase(:,:,iph) = sum(STATE.p(:,:,grains),3);
end

% ------------------------------------------------------------
% dt plot
% ------------------------------------------------------------
if strcmp(w,'dt')

    n = find(DTPHY(:) ~= 0,1,'last');
    if isempty(n); n = numel(DTPHY); end

    plot(DTPHY(1:n),'b.')
    title('dt')
    return

end

% ------------------------------------------------------------
% Phase proportion history
% ------------------------------------------------------------
if strcmp(w,'phase%')

    n = find(any(PHASE ~= 0,2) | TIME(:) ~= 0,1,'last');
    if isempty(n); n = size(PHASE,1); end

    plot(TIME(1:n),PHASE(1:n,1:Nphase),'.-')
    legend(phs_name,'Location','best')
    title('Phase proportion')
    xlabel('time')
    ylabel('phase fraction')
    return

end

% ------------------------------------------------------------
% 2D thermodynamic phase map
% ------------------------------------------------------------
if strcmp(w,'phase2d')

    [~,phase_ID] = max(p_phase,[],3);

    pcolor(GRID.x,GRID.y,phase_ID)
    shading flat
    axis equal tight
    colormap(gca,jet(Nphase))
    caxis([0.5,Nphase+0.5])

    cb = colorbar;
    cb.Ticks = 1:Nphase;
    cb.TickLabels = phs_name;

    title('2D phase map')
    xlabel('x \mum')
    ylabel('y \mum')
    return

end

% ------------------------------------------------------------
% E component, e.g. E1, E2
% ------------------------------------------------------------
tok = regexp(w,'^e(\d+)$','tokens');

if ~isempty(tok)

    ie = str2double(tok{1}{1});
    pcolor(GRID.x,GRID.y,STATE.E{ie})
    shading interp; axis equal tight; colorbar
    title(sprintf('E%d',ie))
    xlabel('x \mum')
    ylabel('y \mum')
    return

end

% ------------------------------------------------------------
% mu_e component, e.g. mu1, mu_e1
% ------------------------------------------------------------
tok = regexp(w,'^mu_?e?(\d+)$','tokens');

if ~isempty(tok)

    ie = str2double(tok{1}{1});
    pcolor(GRID.x,GRID.y,STATE.mu_e{ie})
    shading interp; axis equal tight; colorbar
    title(sprintf('\\mu_e%d',ie))
    xlabel('x \mum')
    ylabel('y \mum')
    return

end

% ------------------------------------------------------------
% p by thermodynamic phase, e.g. p2
% ------------------------------------------------------------
tok = regexp(w,'^p(\d+)$','tokens');

if ~isempty(tok)

    iph = str2double(tok{1}{1});
    pcolor(GRID.x,GRID.y,p_phase(:,:,iph))
    shading interp; axis equal tight; colorbar
    title(sprintf('p %s',phs_name{iph}))
    xlabel('x \mum')
    ylabel('y \mum')
    return

end

% ------------------------------------------------------------
% phi summed by thermodynamic phase, e.g. phi2
% ------------------------------------------------------------
tok = regexp(w,'^phi(\d+)$','tokens');

if ~isempty(tok)

    iph = str2double(tok{1}{1});
    grains = find(phase_index == phase_ids(iph));
    phi_phase = sum(STATE.phi(:,:,grains),3);

    pcolor(GRID.x,GRID.y,phi_phase)
    shading interp; axis equal tight; colorbar
    title(sprintf('\\phi %s',phs_name{iph}))
    xlabel('x \mum')
    ylabel('y \mum')
    return

end

% ------------------------------------------------------------
% c phase/endmember, e.g.
%   c21      = phase 2, endmember 1, masked by phase 2
%   c2_1     = same, useful for multi-digit indices
%   c21all   = plot everywhere without phase mask
%   c2_1all  = plot everywhere without phase mask
% ------------------------------------------------------------
plot_all = false;

tok = regexp(w,'^c(\d+)_(\d+)all$','tokens');

if isempty(tok)
    tok = regexp(w,'^c(\d)(\d)all$','tokens');
end

if ~isempty(tok)
    plot_all = true;
else
    tok = regexp(w,'^c(\d+)_(\d+)$','tokens');

    if isempty(tok)
        tok = regexp(w,'^c(\d)(\d)$','tokens');
    end
end

if ~isempty(tok)

    iph = str2double(tok{1}{1});
    ic  = str2double(tok{1}{2});

    grains = find(phase_index == phase_ids(iph));

    % Phase-weighted c for this thermodynamic phase.
    % This is safer than just taking the first grain.
    c_plot = zeros(size(STATE.p(:,:,1)));
    wsum   = zeros(size(STATE.p(:,:,1)));

    for ig = grains
        wloc   = STATE.p(:,:,ig);
        c_plot = c_plot + wloc .* STATE.c{ig}{ic};
        wsum   = wsum   + wloc;
    end

    good = wsum > eps;
    c_plot(good) = c_plot(good) ./ wsum(good);

    % Mask out regions where this phase is not present
    if ~plot_all
        p_cut = 5e-1;
        c_plot(wsum < p_cut) = NaN;
    end

    pcolor(GRID.x,GRID.y,c_plot)
    shading interp
    axis equal tight
    colorbar

    if plot_all
        title(sprintf('c%d%d %s all',iph,ic,phs_name{iph}))
    else
        title(sprintf('c%d%d %s masked',iph,ic,phs_name{iph}))
    end

    xlabel('x \mum')
    ylabel('y \mum')

    return

end

% ------------------------------------------------------------
% omega, e.g.
%   omg1  = omega of phase 1
%   omg12 = omega phase 1 - phase 2
% also supports omg1_2
% ------------------------------------------------------------
tok = regexp(w,'^omg(\d+)_(\d+)$','tokens');

if isempty(tok)
    tok = regexp(w,'^omg(\d)(\d)$','tokens');
end

if ~isempty(tok)

    iph1 = str2double(tok{1}{1});
    iph2 = str2double(tok{1}{2});

    grains1 = find(phase_index == phase_ids(iph1));
    grains2 = find(phase_index == phase_ids(iph2));

    omg1 = mean(STATE.omg(:,:,grains1),3);
    omg2 = mean(STATE.omg(:,:,grains2),3);

    pcolor(GRID.x,GRID.y,omg1 - omg2)
    shading interp; axis equal tight; colorbar
    title(sprintf('\\omega_{%s} - \\omega_{%s}',phs_name{iph1},phs_name{iph2}))
    xlabel('x \mum')
    ylabel('y \mum')
    return

end

tok = regexp(w,'^omg(\d+)$','tokens');

if ~isempty(tok)

    iph = str2double(tok{1}{1});
    grains = find(phase_index == phase_ids(iph));
    omg_phase = mean(STATE.omg(:,:,grains),3);

    pcolor(GRID.x,GRID.y,omg_phase)
    shading interp; axis equal tight; colorbar
    title(sprintf('\\omega %s',phs_name{iph}))
    xlabel('x \mum')
    ylabel('y \mum')
    return

end


% ------------------------------------------------------------
% ------------------------------------------------------------
% Phase proportion as stacked color blocks with black boundaries
%   x axis: time
%   y axis: cumulative phase fraction from 0 to 1
%
% Call:
%   PF_Plot([3,3,9],'PhaseStack',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%   PF_Plot([3,3,9],'PhaseBlock',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
% ------------------------------------------------------------
if strcmp(w,'phasestack') || strcmp(w,'phaseblock')

    n = find(any(PHASE ~= 0,2) | TIME(:) ~= 0,1,'last');

    if isempty(n)
        n = size(PHASE,1);
    end

    Y = PHASE(1:n,1:Nphase);
    X = TIME(1:n);

    % If TIME is still empty/zero, use iteration number
    if all(X == 0)
        X = (1:n).';
        xlabel_text = 'iteration';
    else
        xlabel_text = 'time';
    end

    % Normalize each row to exactly sum to 1
    s = sum(Y,2);
    good = s > eps;

    for ip = 1:Nphase
        Y(good,ip) = Y(good,ip) ./ s(good);
    end

    % Colored stacked blocks
    h = area(X,Y);
    set(h,'LineStyle','none')

    cmap = jet(Nphase);

    for ip = 1:Nphase
        h(ip).FaceColor = cmap(ip,:);
    end

    hold on

    % Black boundaries between cumulative phase blocks
    Ycum = cumsum(Y,2);

    for ip = 1:Nphase-1
        plot(X,Ycum(:,ip),'k-','LineWidth',1.0)
    end

    % Optional black top and bottom frame
    plot(X,zeros(size(X)),'k-','LineWidth',0.8)
    plot(X,ones(size(X)),'k-','LineWidth',0.8)

    hold off

    ylim([0 1])
    xlim([min(X) max(X)])
    xlabel(xlabel_text)
    ylabel('phase fraction')
    title('Phase proportion')
    legend(h,phs_name,'Location','eastoutside')

    return

end

% ------------------------------------------------------------
% Phase proportion as stacked color blocks with black boundaries
%   x axis: time
%   y axis: cumulative phase fraction from 0 to 1
%
% Call:
%   PF_Plot([3,3,9],'PhaseStack',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%   PF_Plot([3,3,9],'PhaseBlock',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
% ------------------------------------------------------------
if strcmp(w,'phasestack') || strcmp(w,'phaseblock')

    n = find(any(PHASE ~= 0,2) | TIME(:) ~= 0,1,'last');

    if isempty(n)
        n = size(PHASE,1);
    end

    Y = PHASE(1:n,1:Nphase);
    X = TIME(1:n);

    % If TIME is still empty/zero, use iteration number
    if all(X == 0)
        X = (1:n).';
        xlabel_text = 'iteration';
    else
        xlabel_text = 'time';
    end

    % Normalize each row to exactly sum to 1
    s = sum(Y,2);
    good = s > eps;

    for ip = 1:Nphase
        Y(good,ip) = Y(good,ip) ./ s(good);
    end

    % Colored stacked blocks
    h = area(X,Y);
    set(h,'LineStyle','none')

    cmap = jet(Nphase);

    for ip = 1:Nphase
        h(ip).FaceColor = cmap(ip,:);
    end

    hold on

    % Black boundaries between cumulative phase blocks
    Ycum = cumsum(Y,2);

    for ip = 1:Nphase-1
        plot(X,Ycum(:,ip),'k-','LineWidth',1.0)
    end

    % Optional black top and bottom frame
    plot(X,zeros(size(X)),'k-','LineWidth',0.8)
    plot(X,ones(size(X)),'k-','LineWidth',0.8)

    hold off

    ylim([0 1])
    xlim([min(X) max(X)])
    xlabel(xlabel_text)
    ylabel('phase fraction')
    title('Phase proportion')
    legend(h,phs_name,'Location','eastoutside')

    return

end

error('Unknown plot request: %s',what)

end