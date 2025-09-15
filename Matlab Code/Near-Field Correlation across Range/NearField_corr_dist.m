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

 theta0 = 0;
 theta1 = 1;
 relativeRange = logspace(-2,1,8192);
 Fs =[RD/25 RD/5 10*RD];
 
 r2s= logspace(log10(2*D),log10(1e2*RD),1000);
 r2s = relativeRange.*RD;
 %r2s= linspace(2*D,100*RD,8192);
 for L = 1:length(Fs)
 r1 = Fs(L);
 for k = 1:length(r2s)
 r2 = r2s(k);
 at1 = near_field_manifold( N, d, fc, r1, deg2rad(theta0) );
 at2 = near_field_manifold( N, d, fc, r2, deg2rad(theta0) );
 corr(L,k) = abs(at1*at2');
 end
 end

set(groot,'defaultAxesTickLabelInterpreter','latex');
hold on; box on; grid on;
C = linspecer(4);
p1 = plot(relativeRange,corr(1,:), 'color',C(1, :) , 'Linewidth', 2);
plot(0.0286,0.5,'o','color',C(1, :) , 'Linewidth', 2)
plot(0.0667,0.5,'o','color',C(1, :) , 'Linewidth', 2)
annotation('doublearrow',[0.26 0.335],[0.5 0.5])
annotation('textbox',[0.27 0.28 0.2 0.3],'String',{'$r_{BD}$'},'FitBoxToText','on','Interpreter','Latex','fontsize', 16);
annotation('textbox',[0.15 0.25 0.2 0.3],'String',{'$f_{min} $'},'FitBoxToText','on','Interpreter','Latex','fontsize', 16);
annotation('textbox',[0.45 0.25 0.2 0.3],'String',{'$f_{max}$'},'FitBoxToText','on','Interpreter','Latex','fontsize', 16);

p2 = plot(relativeRange,corr(2,:), 'color',C(2, :) , 'Linewidth', 2);
p3 = plot(relativeRange,corr(3,:), 'color',C(3, :) , 'Linewidth', 2);
% p4 = plot(relativeRange,corr(4,:), 'color',C(4, :) , 'Linewidth', 2);
set(gca,'XScale','log');

xlabel('Range','Interpreter','Latex','fontsize', 12)
ylabel('Correlation','Interpreter','Latex','fontsize', 12)
legend([p1,p2,p3],{'f =$r_{RD}$/25', 'f = $r_{RD}$/5' 'f=10$r_{RD}$' }, 'interpreter', 'latex', 'fontsize', 12,Location='best');

xticks([.01 .1 1 10])
xticklabels({'.01$r_{RD}$','0.1$r_{RD}$','1$r_{RD}$','10$r_{RD}$'})
xlim([.01 10])
ylim([0.1 1])

set(gca,'fontsize',18);


%% Code to Verify Beam Depth
%F = Fs;
% % Find the maximum gain of the signal
% max_gain = max(corr);
% 
% % Calculate the 3dB point (assuming signal is in linear scale, not dB)
% dB_point = max_gain / ((2));
% 
% % Find the indices where the signal is closest to the 3dB point
% [~, idx1] = min(abs(corr(1:find(corr == max_gain, 1)) - dB_point));
% [~, idx2] = min(abs(corr(find(corr == max_gain, 1):end) - dB_point));
% idx2 = idx2 + find(corr == max_gain, 1) - 1;
% 
% % Get the values at these points
% val1 = r2s(idx1);
% val2 = r2s(idx2);
% 
% % Subtract the values at these points
% r_BD_A = val2 - val1
% BD =  RD*cosd(theta0).^2 * F * (1/(RD*cosd(theta0).^2-10*F) - 1/(RD*cosd(theta0).^2+10*F))

