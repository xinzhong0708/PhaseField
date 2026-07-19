clear;figure(1);clf
for iload = 4000:200:10000
    load(num2str(iload))


    E_1D = [];
    id = find(STATE.p(2,:,4)>0.999);
    for ie = 1:3
        E_1D(ie,:) = STATE.E{ie}(2,id);
    end
    E_1D = E_1D./sum(E_1D);
    subplot(211);plot(GRID.x(id),E_1D')
    xlim([150,350])
    ylim([0,0.8])
    subplot(212);plot(GRID.x,log10(PARAM.M{1}(2,:)))
    xlim([150,350])
    pause(0.2)


end