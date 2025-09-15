clc
clear all
close all
c = 3e8;
fc = 6e9;
lambda = c/fc;
N = 257;
D = (N-1)*(lambda/2);
RD = 2*D.^2/lambda;
d = lambda/2;
nn = linspace(-N/2, N/2,N);
r00 = linspace(2*D/2,100*RD,1);
%% user is inside the increase decrease arena
r00 =5*D;
theta00 = linspace(-pi/2,pi/2,N);
y1= r00*sin(theta00);
r = zeros(length(theta00),length(r00),length(nn));
for i =1:length(theta00)
    theta0 = theta00(i);
    for j= 1 :length(r00)
      r0 =r00(j);
      r(i,j,:) = sqrt(r0^2 + ((nn-1)*d).^2 - 2*r0*(nn-1)*d*sin(theta0)) - r0;
    end
end
figure
r1 = transpose(squeeze(r));
mesh(theta00,nn,r1);
xlabel('Angle of Arrival')
ylabel('Antenna Index')
zlabel('Distance Difference Function')

%% user is inside the increase  arena
r00 =10*D;
theta00 = linspace(.0500,pi/2,N);
y2= r00*sin(theta00) ;
r = zeros(length(theta00),length(r00),length(nn));
for i =1:length(theta00)
    theta0 = theta00(i);
    for j= 1 :length(r00)
      r0 =r00(j);
      r(i,j,:) = sqrt(r0^2 + ((nn-1)*d).^2 - 2*r0*(nn-1)*d*sin(theta0)) - r0;
    end
end
figure
r2 = transpose(squeeze(r));
mesh(theta00,nn,r2);
xlabel('Angle of Arrival')
ylabel('Antenna Index')
zlabel('Distance Difference Function')

%% user is inside the decrease arena
r00 =10*D;
theta00 = linspace(-pi/2,-.0500,N);
y3= r00*sin(theta00);
r = zeros(length(theta00),length(r00),length(nn));
for i =1:length(theta00)
    theta0 = theta00(i);
    for j= 1 :length(r00)
      r0 =r00(j);
      r(i,j,:) = sqrt(r0^2 + ((nn-1)*d).^2 - 2*r0*(nn-1)*d*sin(theta0)) - r0;
    end
end
figure
r3 = transpose(squeeze(r));
mesh(theta00,nn,r3);
xlabel('Angle of Arrival')
ylabel('Antenna Index')
zlabel('Distance Difference Function')