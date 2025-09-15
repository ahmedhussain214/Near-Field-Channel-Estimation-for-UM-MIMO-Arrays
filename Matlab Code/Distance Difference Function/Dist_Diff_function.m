clear all
clc
close all
fc = 28e9;
c = 3e8;
lambda = c/fc;
array = phased.ULA('NumElements',256,'ElementSpacing',lambda/2);
D = (array.NumElements*array.ElementSpacing);%Antenna Array Dimension
RD = (2*D.^2)/lambda; %Rayleigh Distance
sensorpos = getElementPosition(array);
% viewArray(array)
numele = array.NumElements;
index = linspace(-numele/2, numele/2,numele);

%% When user is located at Dist > D/2
Userpos1 = [25 2 0].';
theta1 = atand((Userpos1(2)-sensorpos(2,129))/(Userpos1(1)-sensorpos(1,128)))
Userpos1_dist = sqrt(sum(Userpos1(1:end).^2)); 
DDF = ((sensorpos(:,1:end) - Userpos1)) ;
Dist1 = sqrt(sum(DDF(:,1:end).^2))- Userpos1_dist;
figure;plot(index, Dist1)
xlabel('Antenna index')
ylabel('Distance Difference function')
title('When user is located at d > D/2')
%% When user is located at Dist < D/2
Userpos2 = [10 -5 0].';
Userpos2_dist = sqrt(sum(Userpos2(1:end).^2)); 
DDF = ((sensorpos(:,1:end) - Userpos2)) ;
Dist2 = sqrt(sum(DDF(:,1:end).^2)) - Userpos2_dist;
figure;
plot(index, Dist2)
xlabel('Antenna index')
ylabel('Distance Difference function')
title('When user is located at d < D/2')
%% When user is located at mod of Dist < D/2
Userpos3 = [10 0 0].';
theta3 = atand((Userpos3(2)-sensorpos(2,129))/(Userpos3(1)-sensorpos(1,128)))
RD3 = RD/10 *cosd(theta3).^2
Userpos3_dist = sqrt(sum(Userpos3(1:end).^2)); 
DDF = ((sensorpos(:,1:end) - Userpos3));
Dist3 = sqrt(sum(DDF(:,1:end).^2))-Userpos3_dist;
figure;
plot(index,Dist3)
xlabel('Antenna index')
ylabel('Distance Difference function')
title('When user is located at mod of Dist < D/2')

