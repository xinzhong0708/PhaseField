clear;hold on
load 10000

id = find(STATE.p(2,:,4)>0.99);

for ie = 1:3
    E_1D(ie,:) = STATE.E{ie}(2,id);
end

E_1D = E_1D./sum(E_1D);

x    = GRID.x(id);
plot(-1+2*(x-x(1))/(x(end)-x(1)),E_1D','k--')
xlim([0 1])