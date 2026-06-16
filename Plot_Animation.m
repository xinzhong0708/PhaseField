clear;figure(1);clf;colormap jet

for iplot = 1400:200:3200
    load(num2str(iplot))
    PF_Plot([3,3,1],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    PF_Plot([3,3,2],'E1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    title('Fe')
    PF_Plot([3,3,3],'E2',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    title('Mg')
    PF_Plot([3,3,4],'E3',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    title('Ca')
    PF_Plot([3,3,5],'E4',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    title('Al')
    % PF_Plot([3,3,6],'E5',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    % title('Na')
    subplot(3,3,6),plot(Tcur-273,Pcur/1e9,'ko');hold on;xlabel('T ^oC');ylabel('P GPa')
    xlim([350,900]);ylim([0.2,1.5]);title('P-T-path')

    PF_Plot([3,3,7],'c51',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    title('Pyrope')
    PF_Plot([3,3,8],'c52',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    title('Almandine')
    PF_Plot([3,3,9],'PhaseStack',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    drawnow
    save_the_image('tiff',100,[25,25],num2str(100000+iplot))
end



% PF_Plot([2,2,1],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
% PF_Plot([2,2,2],'E1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
% title('Fe')
% PF_Plot([2,2,3],'E2',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
% title('Mg')
% PF_Plot([2,2,4],'PhaseStack',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
% drawnow
