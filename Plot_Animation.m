clear;figure(1);clf

for ii = 100:100:5000
    load(num2str(ii))
    PF_Plot([2,2,1],'Phase2d',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    PF_Plot([2,2,2],'E1',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    PF_Plot([2,2,3],'c61',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    PF_Plot([2,2,4],'PhaseStack',STATE,GRID,MODEL,TIME,DTPHY,PHASE)
    drawnow
    save_the_image('tiff',100,[25,25],num2str(100000+ii))
end