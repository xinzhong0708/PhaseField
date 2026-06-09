clear
load Run_L1

%Analytical
xl_eq  = 0.5;
x0     = 0.3;
xs_eq  = 0;
fun    = @(lam) (xl_eq - xs_eq).*lam.*sqrt(pi).*exp(lam.^2).*erfc(lam) - (x0 - xl_eq);

lambda = fzero(fun,[-10,10]);
t      = linspace(0,200000,1000);
s      = 2*lambda*sqrt(PHYS.D_esti*t);
clf
plot(t,abs(s),'k','LineWidth',2)

load Run_L1
hold on
plot(TIME,PHASE(:,1)-0.1,'--','LineWidth',1.5)

load Run_L2
hold on
plot(TIME,PHASE(:,1)-0.1,'--','LineWidth',1.5)

load run_L0_5
hold on
plot(TIME,PHASE(:,1)-0.1,'--','LineWidth',1.5)

load run_L0_2
hold on
plot(TIME,PHASE(:,1)-0.1,'--','LineWidth',1.5)

load run_L0_1
hold on
plot(TIME,PHASE(:,1)-0.1,'--','LineWidth',1.5)

load run_L00_1
hold on
plot(TIME,PHASE(:,1)-0.1,'--','LineWidth',1.5)


legend('Ana','L=1','L=2','L=0.5','L=0.2','L=0.1','L=0.01','location','best');grid on;box on
% set(gca,'xscale','log')
% set(gca,'yscale','log')
xlim([0 2e5])