clear all; 
clc
close all

%% system parameters
N = 256; % number of beams (transmit antennas)
fc = 28e9; % carrier frequency
c = physconst('LightSpeed');
lambda = c/fc;
d = lambda / 2;
D = (N-1)*(lambda/2);
RD = 2*D.^2/lambda;
nn = -(N-1)/2:1:(N-1)/2;



%%  Energy Spread Effect in angle Domain

 theta0 = 30;
 r0 = RD/100; 
%Channel Model Near Field
 r = sqrt(r0^2 + (nn*d).^2 - 2*r0*nn*d*sind(theta0));
 at = exp(-1j*2*pi*fc*(r - r0)/c)/sqrt(N);
 %FFT Based Angle Estimation
 K = N;
 Y = fftshift(abs(fft(at,K)));
 angle_scale = -90:180/K:90 -180/K;
 [~, max_index] = max(Y);
 theta_axis=asind(angle_scale/180*lambda/d);
 angle_estimate = (theta_axis(max_index));
 figure
 plot(theta_axis,pow2db(Y.^2),'Linewidth',2,'markersize',5,'color','[0.8500 0.3250 0.0980]','MarkerFaceColor','w')
 axis([-60 60 -35 15])
 xlabel('DOA (degree)','Interpreter','Latex')
 ylabel('Power Spectrum [dB]','Interpreter','Latex')
 grid on
 set(gca,'fontsize',18);
 a_1 = gcf;
 exportgraphics(a_1,'energy_spread.eps','ContentType','vector')
 exportgraphics(a_1,'energy_spread.jpg')
%% 1D Angle Estimae error vs Distance

r01 = linspace(1,RD,1000);
Error_1D = zeros(length(r01),1);
theta0 = 10;
error_pred_dist = 0.5*N*d*cotd(theta0)*cosd(theta0);
for i = 1: length(r01)
    r0 = r01(i); 
%Channel Model Near Field
 r = sqrt(r0^2 + (nn*d).^2 - 2*r0*nn*d*sind(theta0));
 at = exp(-1j*2*pi*fc*(r - r0)/c)/sqrt(N);
 % FFT Based Angle Estimation
 K = N;
 Y = fftshift(abs(fft(at,K)));
 angle_scale = -90:180/K:90 -180/K;
 [~, max_index] = max(Y);
 theta_axis=asind(angle_scale/180*lambda/d);
 angle_estimate = (theta_axis(max_index));
 Error_1D(i) = sqrt((angle_estimate - theta0).^2);
end
figure
plot(r01,pow2db(Error_1D),'Linewidth',2,'markersize',5,'color','[0.6350 0.0780 0.1840]','MarkerFaceColor','w')
xlabel('Distance (meter)','Interpreter','Latex')
xticks([  RD/10^2  RD/10  RD ])
xticklabels({'.01r_{RD}','.1r_{RD}','0.2r_{RD}'})
ylabel('Angle Estimation Error [dB]','Interpreter','Latex')
ylim([-10, 15])
set(gca,'fontsize',18);
set(gca,"XScale","log")
grid on

 a_2 = gcf;
 exportgraphics(a_2,'angle_estimation_error_1D.eps','ContentType','vector')
 exportgraphics(a_2,'angle_estimation_error_1D.jpg')
%% Angle Estimate Error in 2D (Angle and Distance) 
relativeRange = logspace(-2.5,1,256);

theta01 = linspace(-80,80,256);
%r01 = logspace(1,RD/5,256);
r01 = relativeRange.*RD;
Error = zeros(length(theta01),length(r01));
for j=1:length(theta01)
    theta0 = theta01(j);
for i = 1: length(r01)
    r0 = r01(i); 
%% Channel Model Far Field
bt = 1/sqrt(N) * exp(1j*pi*sind(theta0)*nn);
%% Channel Model Near Field
 r = sqrt(r0^2 + (nn*d).^2 - 2*r0*nn*d*sind(theta0));
 at = exp(-1j*2*pi*fc*(r - r0)/c)/sqrt(N);
 %% FFT Based Angle Estimation
 K = N;
 Y = fftshift(abs(fft(at,K)));
 angle_scale = -90:180/K:90 -180/K;
 [~, max_index] = max(Y);
 theta_axis=asind(angle_scale/180*lambda/d);
 angle_estimate = (theta_axis(max_index));
 Error(j,i) = sqrt((angle_estimate - theta0).^2);
 % figure
 % plot(theta_axis,pow2db(Y.^2))
 % axis([-90 90 -20 1.1*pow2db(N)])
 % xlabel('DOA \theta/degree')
 % ylabel('Power Spectrum/dB')
 % title('FFT Based Angle estimation')
 % grid on
end
% figure
% plot(r01,pow2db(Error),LineWidth=2)
% xlabel('Distance in Meters')
% ylabel('Angle Estimate Error in dB ')
end
Error = pow2db(Error);
set(groot,'defaultAxesTickLabelInterpreter','latex');
figure
imagesc(relativeRange,theta01,Error);
set(gca,'XScale','log');
colormap("winter")
a = colorbar;
a.Label.String = ' Angular Spread [dB]';
a.Label.Interpreter = 'Latex';
xlabel('Propagation Distance','Interpreter','Latex')
xticks([ .01 0.1 1 10   ])
xticklabels({'.01$r_{RD}$','0.1$r_{RD}$','1$r_{RD}$','10$r_{RD}$'})
%yticks([   -80 -60 -30  0 30 60 80  ])
yticks([   -80   0  80  ])
yticklabels ({'$-\pi/2$', '-$\pi/3$','$-\pi/6$','0','$\pi/6$','$\pi/3$','$\pi/2$' });
yticklabels ({'$-\pi/2$','0','$\pi/2$' });

ylabel('Angle','Interpreter','Latex')
set(gca,'fontsize',18);
 % a_3 = gcf;
 % exportgraphics(a_3,'angle_estimation_error_2D.eps','ContentType','vector')
 % exportgraphics(a_3,'angle_estimation_error_2D.jpg')
