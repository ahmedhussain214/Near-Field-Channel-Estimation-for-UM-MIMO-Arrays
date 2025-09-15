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
 r0s = [ 10*RD RD RD /5 RD/9 RD/20 RD/40 RD/80] ; 
 r0s = flip( r0s) 
C = linspecer(7);
%Channel Model Near Field
  for m = 1:length(r0s)
    r0 = r0s(m);
    r = sqrt(r0^2 + (nn*d).^2 - 2*r0*nn*d*sind(theta0));
    b = exp(-1j*2*pi*fc*(r - r0)/c)/sqrt(N);


 %FFT Based Angle Estimation
 K = 2*N;
 Y = fftshift(abs(fft(b,K)));
 angle_scale = -90:180/K:90 -180/K;
 [~, max_index] = max(Y);
 theta_axis=asind(angle_scale/180*lambda/d);
 Y = Y./max(Y);
  % Y = Y/max(Y); 
  % Y1  = Y1/max(Y1);
 h1(m,:) = pow2db(Y.^2);
 

  end

p= ribbon(angle_scale,h1.',.025)
set(groot,'defaultAxesTickLabelInterpreter','latex');

% p.EdgeColor = 'b';
p(1).MeshStyle = "column";
p(2).MeshStyle = "column";
p(3).MeshStyle = "column";
p(4).MeshStyle = "column";
p(5).MeshStyle = "column";
p(6).MeshStyle = "column";
p(7).MeshStyle = "column";

p(1).LineStyle ="none";
p(2).LineStyle ="none";
p(3).LineStyle ="none";
p(4).LineStyle ="none";
p(5).LineStyle ="none";
p(6).LineStyle ="none";
p(7).LineStyle ="none";
pbaspect([2  1 .8])
% p.LineJoin = "miter";

zlim([-40 0])
ylim([-60 60])
xlabel('Range','Interpreter','Latex')
ylabel('Angle','Interpreter','Latex')
zlabel('Normalized Gain (dB)','Interpreter','Latex')

ax = gca;
set(gca,'XTick',linspace(1,7,7))
ax.XTickLabel = {'0.01$r_{RD}$','0.02$r_{RD}$', '0.05$r_{RD}$','0.1$r_{RD}$','$0.2r_{RD}$','$r_{RD}$', '10$r_{RD}$'};
ax.YTickLabel = {'-$\pi/2$','0','$\pi/2$' };
set(gca,'fontsize',12);
colormap("winter")

