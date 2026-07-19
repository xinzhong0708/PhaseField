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
%   PF_Plot([3,3,7],'c51',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%   PF_Plot([3,3,8],'e51',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%   PF_Plot([3,3,9],'Phase%',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
%
% Notes:
%   c51 means thermodynamic phase 5, composition variable 1.
%   e51 means thermodynamic phase 5, element variable 1.
%   For multi-digit indices use c5_10 or e5_10.
%   c51/e51 are phase-weighted over all grains belonging to phase 5.
%   normal c plots hide pixels with phase p < 1e-2; c51all plots everywhere.
%   omega plots are unmasked and shown on all grids.

ax = subplot(pos(1),pos(2),pos(3));

% Store the original subplot slot once, before a colorbar can resize it.
% Reuse this position during later plot updates.
tile_pos = getappdata(ax,'PFPlotTilePosition');

if isempty(tile_pos)
    tile_pos = ax.Position;
    setappdata(ax,'PFPlotTilePosition',tile_pos)
end

cla(ax,'reset')
setappdata(ax,'PFPlotTilePosition',tile_pos)

% Finalize_Plot restores the same plotting area for every panel and places
% colorbars in a reserved strip to the right of the axes.
cb_tag = sprintf('PFPlotColorbar_%d_%d_%d',pos(1),pos(2),pos(3));

old_cb = findall(gcf,'Tag',cb_tag);
if ~isempty(old_cb)
    delete(old_cb)
end

setappdata(ax,'PFPlotColorbar',[]);
setappdata(ax,'PFPlotColorbarTag',cb_tag);

cleanup_plot = onCleanup(@() Finalize_Plot(ax,tile_pos));

what0 = strtrim(what);
w     = lower(what0);

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

%c plotting cutoff. Normal c plots hide regions where phase amount is small.
%c...all still plots everywhere.
c_plot_p_cut    = 1e-1;
omg_local_p_cut = 1e-3;

% ------------------------------------------------------------
% dt plot
% ------------------------------------------------------------
if strcmp(w,'dt')

    n = find(DTPHY(:) ~= 0,1,'last');

    if isempty(n)
        n = numel(DTPHY);
    end

    plot(ax,1:n,DTPHY(1:n),'.-','LineWidth',0.9,'MarkerSize',7)
    grid(ax,'on')
    xlabel(ax,'iteration')
    ylabel(ax,'dt')
    title(ax,'Time step')
    return

end

% ------------------------------------------------------------
% Phase proportion history
% ------------------------------------------------------------
if strcmp(w,'phase%')

    n = find(any(PHASE ~= 0,2) | TIME(:) ~= 0,1,'last');

    if isempty(n)
        n = size(PHASE,1);
    end

    plot(ax,TIME(1:n),PHASE(1:n,1:Nphase),'.-','LineWidth',0.9,'MarkerSize',6)
    grid(ax,'on')
    legend(ax,phs_name,'Location','best','Box','off')
    title(ax,'Phase proportion')
    xlabel(ax,'time')
    ylabel(ax,'phase fraction')
    return

end

% ------------------------------------------------------------
% Phase proportion as stacked color blocks with black boundaries
% ------------------------------------------------------------
if strcmp(w,'phasestack') || strcmp(w,'phaseblock')

    n = find(any(PHASE ~= 0,2) | TIME(:) ~= 0,1,'last');

    if isempty(n)
        n = size(PHASE,1);
    end

    Y = PHASE(1:n,1:Nphase);
    X = TIME(1:n);

    if all(X == 0)
        X = (1:n).';
        xlabel_text = 'iteration';
    else
        xlabel_text = 'time';
    end

    s = sum(Y,2);
    good = s > eps;

    for ip = 1:Nphase
        Y(good,ip) = Y(good,ip)./s(good);
    end

    h = area(ax,X,Y);
    set(h,'LineStyle','none')

    cmap = jet(Nphase);

    for ip = 1:Nphase
        h(ip).FaceColor = cmap(ip,:);
    end

    hold(ax,'on')

    Ycum = cumsum(Y,2);

    for ip = 1:Nphase-1
        plot(ax,X,Ycum(:,ip),'k-','LineWidth',0.8)
    end

    plot(ax,X,zeros(size(X)),'k-','LineWidth',0.7)
    plot(ax,X,ones(size(X)),'k-','LineWidth',0.7)

    hold(ax,'off')
    grid(ax,'on')

    ylim(ax,[0 1])
    if numel(X) > 1 && max(X) > min(X)
        xlim(ax,[min(X) max(X)])
    end
    xlabel(ax,xlabel_text)
    ylabel(ax,'phase fraction')
    title(ax,'Phase proportion')
    % legend(h,phs_name,'Location','eastoutside')
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

    cb = Add_Colorbar(ax);
    cb.Ticks = 1:Nphase;
    cb.TickLabels = phs_name;

    title('2D phase map')
    xlabel('x \mum')
    ylabel('y \mum')
    return

end

% ------------------------------------------------------------
% Conserved E component, e.g. E1
% Use uppercase E to avoid conflict with phase-resolved e51.
% Lowercase e1 is kept as a backward-compatible alias.
% ------------------------------------------------------------
tok = regexp(what0,'^E(\d+)$','tokens');

if isempty(tok)
    tok = regexp(w,'^e(\d+)$','tokens');
end

if ~isempty(tok)

    ie = str2double(tok{1}{1});

    pcolor(GRID.x,GRID.y,STATE.E{ie})
    shading interp
    axis equal tight
    Add_Colorbar(ax);
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
    shading interp
    axis equal tight
    Add_Colorbar(ax);
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
    Check_Phase_Index(iph,Nphase,what0)

    pcolor(GRID.x,GRID.y,p_phase(:,:,iph))
    shading interp
    axis equal tight
    Add_Colorbar(ax);
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
    Check_Phase_Index(iph,Nphase,what0)

    grains    = find(phase_index == phase_ids(iph));
    phi_phase = sum(STATE.phi(:,:,grains),3);

    pcolor(GRID.x,GRID.y,phi_phase)
    shading interp
    axis equal tight
    Add_Colorbar(ax);
    title(sprintf('\\phi %s',phs_name{iph}))
    xlabel('x \mum')
    ylabel('y \mum')
    return

end

% ------------------------------------------------------------
% c phase/component, e.g.
%   c51      = phase 5, component 1, all grains of phase 5
%   c5_1     = same, useful for multi-digit indices
%   c51all   = no phase mask
%   c5_1all  = no phase mask
% ------------------------------------------------------------
[ok,iph,ic,plot_all] = Parse_Phase_Component(w,'c');

if ok

    Check_Phase_Index(iph,Nphase,what0)

    grains = find(phase_index == phase_ids(iph));

    if isempty(grains)
        error('PF_Plot: phase %d has no grains.',iph)
    end

    c_plot = PhaseWeighted_CellField(STATE.c,STATE.p,grains,ic);

    if ~plot_all
        c_plot(p_phase(:,:,iph) < c_plot_p_cut) = NaN;
    end

    Plot_Field(ax,GRID,c_plot)

    if plot_all
        title(sprintf('c%d%d %s all grains',iph,ic,phs_name{iph}))
    else
        title(sprintf('c%d%d %s',iph,ic,phs_name{iph}))
    end

    return

end

% ------------------------------------------------------------
% e phase/component, e.g.
%   e51      = phase 5, element/component 1, all grains of phase 5
%   e5_1     = same, useful for multi-digit indices
%   e51all   = no phase mask
%   e5_1all  = no phase mask
% ------------------------------------------------------------
[ok,iph,ie,plot_all] = Parse_Phase_Component(w,'e');

if ok

    Check_Phase_Index(iph,Nphase,what0)

    if ~isfield(STATE,'e') || isempty(STATE.e)
        error('PF_Plot: STATE.e is empty. Run Calc_e or LE first.')
    end

    grains = find(phase_index == phase_ids(iph));

    if isempty(grains)
        error('PF_Plot: phase %d has no grains.',iph)
    end

    e_plot = PhaseWeighted_CellField(STATE.e,STATE.p,grains,ie);

    if ~plot_all
        e_plot(p_phase(:,:,iph) <= eps) = NaN;
    end

    Plot_Field(ax,GRID,e_plot)

    if plot_all
        title(sprintf('e%d%d %s all grains',iph,ie,phs_name{iph}))
    else
        title(sprintf('e%d%d %s',iph,ie,phs_name{iph}))
    end

    return

end

% ------------------------------------------------------------
% ------------------------------------------------------------
% omega difference on coexisting grids only, e.g.
%   omg12   = local omega phase 1 - phase 2, only where both coexist
%   omg1_2  = same
% ------------------------------------------------------------
tok = regexp(w,'^omg(\d+)_(\d+)$','tokens');

if isempty(tok)
    tok = regexp(w,'^omg(\d)(\d)$','tokens');
end

if ~isempty(tok)

    iph1 = str2double(tok{1}{1});
    iph2 = str2double(tok{1}{2});

    Check_Phase_Index(iph1,Nphase,what0)
    Check_Phase_Index(iph2,Nphase,what0)

    grains1 = find(phase_index == phase_ids(iph1));
    grains2 = find(phase_index == phase_ids(iph2));

    if isempty(grains1) || isempty(grains2)
        error('PF_Plot: one of the requested phases has no grains.')
    end

    % Local phase-weighted omega, not global phase average
    omg1 = PhaseWeighted_ArrayField(STATE.omg,STATE.p,grains1);
    omg2 = PhaseWeighted_ArrayField(STATE.omg,STATE.p,grains2);

    % Only show true local coexistence grids
    omg_p_cut = 1e-3;
    coexist = p_phase(:,:,iph1) > omg_p_cut & ...
              p_phase(:,:,iph2) > omg_p_cut;

    domg = omg1 - omg2;
    domg(~coexist) = NaN;

    Plot_Field(ax,GRID,domg)
    title(sprintf('\\omega_{%s} - \\omega_{%s}', ...
        phs_name{iph1},phs_name{iph2}))
    return

end


tok = regexp(w,'^omg(\d+)(avg)?$','tokens');

if ~isempty(tok)

    iph = str2double(tok{1}{1});
    opt = '';

    if numel(tok{1}) >= 2 && ~isempty(tok{1}{2})
        opt = tok{1}{2};
    end

    Check_Phase_Index(iph,Nphase,what0)

    grains = find(phase_index == phase_ids(iph));

    if isempty(grains)
        error('PF_Plot: phase %d has no grains.',iph)
    end

    if strcmp(opt,'avg')
        omg_phase = PhaseMean_ArrayField(STATE.omg,grains);
        Plot_Field(ax,GRID,omg_phase)
        title(sprintf('\\omega %s avg grains',phs_name{iph}))
    else
        omg_phase = LocalGrain_ArrayField(STATE.omg,STATE.p,grains,omg_local_p_cut);
        Plot_Field(ax,GRID,omg_phase)
        title(sprintf('\\omega %s local grains',phs_name{iph}))
    end

    return

end

error('Unknown plot request: %s',what)

end

% =========================================================================
% Local helpers
% =========================================================================
function [ok,iph,ic,plot_all] = Parse_Phase_Component(w,prefix)

ok       = false;
iph      = [];
ic       = [];
plot_all = false;

expr = ['^' prefix '(\d+)_(\d+)(all)?$'];
tok  = regexp(w,expr,'tokens');

if isempty(tok)
    expr = ['^' prefix '(\d)(\d)(all)?$'];
    tok  = regexp(w,expr,'tokens');
end

if isempty(tok)
    return
end

ok  = true;
iph = str2double(tok{1}{1});
ic  = str2double(tok{1}{2});

if numel(tok{1}) >= 3 && ~isempty(tok{1}{3})
    plot_all = true;
end

end


function A = PhaseWeighted_CellField(F,p,grains,ic)

[ny,nx,~] = size(p);

A    = zeros(ny,nx);
wsum = zeros(ny,nx);

for ig = grains

    if ic > numel(F{ig})
        error('PF_Plot: component %d does not exist for grain %d.',ic,ig)
    end

    wloc = p(:,:,ig);
    A    = A    + wloc.*F{ig}{ic};
    wsum = wsum + wloc;

end

good = wsum > eps;

A(good)  = A(good)./wsum(good);
A(~good) = NaN;

end


function A = PhaseWeighted_ArrayField(F,p,grains)

[ny,nx,~] = size(p);

A    = zeros(ny,nx);
wsum = zeros(ny,nx);

for ig = grains

    wloc = p(:,:,ig);
    A    = A    + wloc.*F(:,:,ig);
    wsum = wsum + wloc;

end

good = wsum > eps;

A(good)  = A(good)./wsum(good);
A(~good) = NaN;

end



function A = PhaseMean_ArrayField(F,grains)

if isempty(grains)
    error('PF_Plot: empty grain list.')
end

A = mean(F(:,:,grains),3);

end


function A = LocalGrain_ArrayField(F,p,grains,p_cut)

[ny,nx,~] = size(p);

idmap = LocalGrainID_ByP_Local(p,grains,p_cut);
idmap = Fill_GrainID_Nearest_Local(idmap);

A = NaN(ny,nx);

for ig = grains

    mask = idmap == ig;

    if any(mask(:))
        tmp = F(:,:,ig);
        A(mask) = tmp(mask);
    end
end

end


function idmap = LocalGrainID_ByP_Local(p,grains,p_cut)

[ny,nx,~] = size(p);

if isempty(grains)
    idmap = zeros(ny,nx);
    return
end

[pmax,iloc] = max(p(:,:,grains),[],3);

idmap = zeros(ny,nx);

for k = 1:numel(grains)
    idmap(iloc == k & pmax > p_cut) = grains(k);
end

end


function idmap = Fill_GrainID_Nearest_Local(idmap)

[ny,nx] = size(idmap);

if ~any(idmap(:) > 0)
    return
end

for it = 1:(ny+nx)

    missing = idmap == 0;

    if ~any(missing(:))
        break
    end

    old = idmap;

    L = idmap(:,[1,1:nx-1]);
    R = idmap(:,[2:nx,nx]);
    U = idmap([1,1:ny-1],:);
    D = idmap([2:ny,ny],:);

    take = missing & L > 0;
    idmap(take) = L(take);

    missing = idmap == 0;
    take = missing & R > 0;
    idmap(take) = R(take);

    missing = idmap == 0;
    take = missing & U > 0;
    idmap(take) = U(take);

    missing = idmap == 0;
    take = missing & D > 0;
    idmap(take) = D(take);

    if isequal(old,idmap)
        break
    end
end

end


function mask_pair = PairNearMask_Local(p,grains1,grains2,p_cut,thick)

p1 = sum(p(:,:,grains1),3);
p2 = sum(p(:,:,grains2),3);

m1 = p1 > p_cut;
m2 = p2 > p_cut;

m1 = DilateMask_Local(m1,thick);
m2 = DilateMask_Local(m2,thick);

mask_pair = m1 & m2;

end


function mask = DilateMask_Local(mask,thick)

if thick <= 0
    return
end

ker  = ones(2*thick+1,2*thick+1);
mask = conv2(double(mask),ker,'same') > 0;

end


function Plot_Field(ax,GRID,A)

pcolor(ax,GRID.x,GRID.y,A)
shading(ax,'interp')
axis(ax,'equal')
axis(ax,'tight')
Add_Colorbar(ax);
xlabel(ax,'x \mum')
ylabel(ax,'y \mum')

end



function cb = Add_Colorbar(ax)
%ADD_COLORBAR Create one colorbar for the current PF_Plot panel.
%
% The axes size is restored later by Finalize_Plot, so the colorbar does
% not make field plots smaller than history or phase-stack plots.

old_cb = getappdata(ax,'PFPlotColorbar');

if isgraphics(old_cb)
    delete(old_cb)
end

cb = colorbar(ax);

cb_tag = getappdata(ax,'PFPlotColorbarTag');
if ~isempty(cb_tag)
    cb.Tag = cb_tag;
end

setappdata(ax,'PFPlotColorbar',cb);

end


function Finalize_Plot(ax,tile_pos)
%FINALIZE_PLOT Apply one consistent panel size and visual style.
%
% All panels reserve the same narrow strip on the right. Field plots use
% this strip for a colorbar; line and area plots leave it empty. Therefore,
% phasestack, phase-history and field panels have the same axes size.

if ~isgraphics(ax)
    return
end

fig = ancestor(ax,'figure');
set(fig,'Color','w')

set(ax, ...
    'Units','normalized', ...
    'PositionConstraint','innerposition', ...
    'Box','on', ...
    'Layer','top', ...
    'FontName','Arial', ...
    'FontSize',9, ...
    'LineWidth',0.8, ...
    'TickDir','out', ...
    'TickLength',[0.015 0.015])

% Reserve identical colorbar space for every subplot.
reserve_fraction = 0.18;

ax_pos    = tile_pos;
ax_pos(3) = tile_pos(3)*(1-reserve_fraction);
ax.Position = ax_pos;

% Consistent title and axis-label typography.
if isgraphics(ax.Title)
    set(ax.Title,'FontName','Arial','FontSize',10,'FontWeight','normal')
end

if isgraphics(ax.XLabel)
    set(ax.XLabel,'FontName','Arial','FontSize',9)
end

if isgraphics(ax.YLabel)
    set(ax.YLabel,'FontName','Arial','FontSize',9)
end

% Place the colorbar manually in the reserved strip. Resetting ax.Position
% above prevents MATLAB from shrinking the plotted panel.
cb = getappdata(ax,'PFPlotColorbar');

if isgraphics(cb)

    cb.Units = 'normalized';

    gap      = 0.025*tile_pos(3);
    cb_width = 0.045*tile_pos(3);
    cb_y     = ax_pos(2) + 0.04*ax_pos(4);
    cb_h     = 0.92*ax_pos(4);

    cb.Position = [ ...
        ax_pos(1) + ax_pos(3) + gap, ...
        cb_y, ...
        cb_width, ...
        cb_h];

    set(cb, ...
        'FontName','Arial', ...
        'FontSize',8, ...
        'LineWidth',0.7, ...
        'TickDirection','out')
end

end


function Check_Phase_Index(iph,Nphase,what)

if iph < 1 || iph > Nphase
    error('PF_Plot: phase index %d in "%s" is outside 1:%d.',iph,what,Nphase)
end

end
