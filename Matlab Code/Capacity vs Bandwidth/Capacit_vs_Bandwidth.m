clc
clear all
close all


B = 20e6: 10e6: 2e9;
h=1e-14;
P=8;
T = 290;
k = physconst('Boltzmann');
No=k*T;
C = B.*log2(1+(P.*h)./(B.*No));
figure
plot(B/(1e6),C/(1e6),LineWidth=2,Color='red');
xlabel('Bandwidth [MHz]','Interpreter','Latex')
ylabel('Capacity [Mb/s]','Interpreter','Latex')
set(gca,'fontsize',18);
set(gca,"YScale","log")
ylim([18, 30]);
grid on