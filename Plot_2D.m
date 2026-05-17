clear; figure(1); clf; colormap jet

phase_name = {'Olv','Cpx','Grt','Qtz'};
nphase     = numel(phase_name);

% Discrete colormap
cmap = jet(nphase);
colormap(cmap)

interval = 100;
plot_no  = 14000;

for i = interval:interval:plot_no

    load(num2str(i))

    [~,phase_ID] = max(STATE.p,[],3);

    pcolor(GRID.x*1e6,GRID.y*1e6,phase_ID);
    axis equal;
    xlim([min(GRID.x*1e6) max(GRID.x*1e6)])
    ylim([min(GRID.y*1e6) max(GRID.y*1e6)])
    shading interp

    % Force phase IDs to map cleanly to discrete colors
    clim([0.5 nphase+0.5])

    xlabel('x \mum')
    ylabel('y \mum')
    title(sprintf('Step %d',i))

    % Legend on top
    hold on
    hleg = gobjects(1,nphase);
    for ip = 1:nphase
        hleg(ip) = plot(nan,nan,'s', ...
            'MarkerFaceColor',cmap(ip,:), ...
            'MarkerEdgeColor',cmap(ip,:), ...
            'MarkerSize',8);
    end
    hold off

    legend(hleg,phase_name, 'Location','northoutside', 'Orientation','horizontal');

    drawnow

    save_the_image('tiff',100,[9,9],num2str(1e6+i))

end