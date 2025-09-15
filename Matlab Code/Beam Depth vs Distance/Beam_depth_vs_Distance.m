clc
clear all
close all

fc = 28e9; % carrier frequency
N = 256;
c = 3e8;
lambda_c = c/fc;
d = lambda_c / 2;
A = (N-1).*d;
RD = (2*A.^2)/lambda_c;

relativeRange = linspace(0.02,0.1,5);
F = relativeRange*RD;

thetas = linspace(-90,90,128);
 
r = zeros(length(thetas),length(F));


for k =1:length(thetas) 
    theta = thetas(k);
for i = 1:length(F)
f = F(i);
%r(i) =  RD * f * (1/(RD-10*f) - 1/(RD+10*f))
if f >= RD*(cosd(theta))^2 / 10
r(k,i) =1;
else 
r(k,i) =  20*(f^2)*RD*(cosd(theta))^2 / (RD^2 * (cosd(theta))^4 - 100*f^2 );  
end
end
end
% plot(F, r,LineWidth=2)
% xlabel('Focus distance (f)','Interpreter','Latex')
% ylabel('Depth of focus ($r_{BD}$)','Interpreter','Latex')
% xline(RD/10,'--b',{'Rayeligh distance/10'});
% set(gca,'fontsize',18);
% set(gca,"YScale","log")
% grid on
% box on
% ylim([0 max(r)+100])
% ax = gca;
% exportgraphics(ax,'Beamdepth_vs_focus.eps','ContentType','vector')
% exportgraphics(ax,'Beamdepth_vs_focus.jpg')
%max_BD = max(max(r));
%r1 = changem( r , max_BD , 0 );
%r1 = (r1./max_BD);

% figure
% ab = bar3(thetas,log10(r),1.6, "white");
% xlabel('Focus distance (f)','Interpreter','Latex')
% zlabel('Depth of focus ($r_{BD}$)','Interpreter','Latex')
% ylabel('Angle ($\theta$)','Interpreter','Latex')
% colorbar

set(groot,'defaultAxesTickLabelInterpreter','latex');
figure
p = ribbon(thetas,r,0.05,"white");
xlabel('Focus distance (f)','Interpreter','Latex')
zlabel('Depth of focus ($r_{BD}$)','Interpreter','Latex')
ylabel('Angle ($\theta$)','Interpreter','Latex')

%set(gca,'XScale','log');
set(gca,'ZScale','log');
colormap("parula")
p(1).MeshStyle = "column";
p(2).MeshStyle = "column";
p(3).MeshStyle = "column";
p(4).MeshStyle = "column";
p(5).MeshStyle = "column";

p(1).LineStyle ="none";
p(2).LineStyle ="none";
p(3).LineStyle ="none";
p(4).LineStyle ="none";
p(5).LineStyle ="none";

pbaspect([1.5  1 0.75])


yticks(linspace(-90,90,3))
yticklabels({'-$\pi/2$','0','$\pi/2$'});

xticks([1 2 3 4 5])
xticklabels({ '0.02$r_{RD}$','0.04$r_{RD}$','0.06$r_{RD}$','0.08$r_{RD}$','0.10$r_{RD}$'});
xlim([1 5])
set(gca,'fontsize',14);