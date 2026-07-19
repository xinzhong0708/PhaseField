clear;figure(1);clf;colormap jet
for ii = 100:100:14600
    %Load
    load(num2str(ii))
    %Plot
    PF_Plot([2,4,1],'E1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    title('FeO')
    PF_Plot([2,4,2],'E2',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    title('MgO')
    PF_Plot([2,4,3],'E3',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    title('CaO')
    PF_Plot([2,4,4],'E4',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    title('Al2O3')
    PF_Plot([2,4,5],'c11',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    PF_Plot([2,4,6],'c41',STATE,GRID,MODEL,TIME,DTPHY,PHASE)


    PF_Plot([2,4,7],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    title('Phase proportion')
    PF_Plot([2,4,8],'phasestack',STATE,GRID,MODEL,TIME,DTPHY,PHASE)

    save_the_image('tiff',100,[22,11],num2str(100000+ii))
    drawnow


end