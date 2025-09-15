clc
clear all
close all

fc = 28e9; % carrier frequency
N = 256;
c = 3e8;
lambda = c/fc;
d = lambda/2;
D = (N-1).*d;
RD = (2*D.^2)/lambda;


fc= 3e8/lambda;
theta0 = linspace(-pi/2,pi/2, N);

 % D = 1 * N;
 % theta0 = -1 + 2/N : 2/N : 1;


label = [];
for  i =1:length(theta0)
F = 2*D;
theta = theta0(i);
rr= [];

while F <= (RD/10)*(cos(theta).^2)
rr = [F rr];

% r =  RD*  F * (1/(RD-10*F) - 1/(RD+10*F));
r =  RD*cos(theta).^2 * F * (1/(RD*cos(theta).^2-10*F) - 1/(RD*cos(theta).^2+10*F));

F = F+r/2;

end
t =(RD/10)*(cos(theta).^2);
if t>2*D
    v = t ;
else
    v =[];
end
rr =[rr   v];
label1 = [ones(1,length(rr))*(theta) ; rr ];
label = [label label1];
end
dict = zeros(N,size(label,2));
for i= 1:size(label,2)
Nt = N;
theta = label(1,i);
r = label(2,i);
% at = near_field_manifold( Nt, d, fc, r, theta );
at = polar_domain_manifold( Nt, d, fc, r, theta );
dict(:,i) =at;
end

set(groot,'defaultAxesTickLabelInterpreter','latex');
figure
c  = sin(linspace(0,pi,length(label)));
scatter([label(1,:)],[label(2,:)],15,c,'filled');
xlabel('Angle','Interpreter','Latex');
ylabel('Distance','Interpreter','Latex');
yline(130,'--',{'Rayleigh distance'},'fontsize',14,'Interpreter','Latex');
yline(RD/10,'--',{'Rayleigh distance /10'},'fontsize',14,'Interpreter','Latex');
box on

ax = gca;
colormap(gca,"lines")


yticks([0 floor(RD/10) floor(130)])
yticklabels({'$2D$','$\frac{r_{RD}}{10}$','$r_{RD}$'})
yscale log
ylim([2*D 180])

xlim([-pi/2.4 pi/2.4])
xticks([-pi/2.4 0 pi/2.4])
xticklabels({'$-\pi/2$','$0$','$\pi/2$'})

set(gca,'fontsize',18);


