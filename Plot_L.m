
clf


name = {'run_L7.mat','run_L8.mat','run_L9.mat','run_L10.mat','run_L11.mat'};
subplot(121);hold on
for ii = 1:length(name)
    load(name{ii})
    plot(TIME,PHASE1,'-')

    p1 = polyfit(TIME(200:end),PHASE1(200:end),1);

    L(ii) = PHYS.L;
    S(ii) = p1(1);

end
set(gca,'xscale','lin')

subplot(122)

plot(L,S,'o')