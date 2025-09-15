syms s
G = abs((fresnelc(s) + 1i*fresnels(s))/s) ;
figure
fplot(G,[0 10])
grid on
hold on
plot([1 1]*1.57, ylim, '--k', Color="r")               
hold off
xlabel('Gamma $(\gamma)$','Interpreter','Latex')
ylabel('Error Function $(\gamma)$','Interpreter','Latex')
grid on
set(gca,'fontsize',18);