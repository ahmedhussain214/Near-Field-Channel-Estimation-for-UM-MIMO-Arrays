clc
clear all; 
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
SNR_dB = 10;  SNR_linear=10.^(SNR_dB/10.);
sigma2=1/SNR_linear;
L=3;
fs = 100e6;
M=1;
sector = deg2rad(60);
T_Rmin = [34 600];


for iii = 1:length(T_Rmin)
Rmin = T_Rmin(iii);   
Rmax = Rmin;
%% LOS Angle/Range Estimation
 
% Channel Model Near Field
%  theta0 = 35;  
%  r0_max = (cosd(theta0)).^2 * RD/10;
%  r0_min = 2*D;
%  r0 = r0_min + (r0_max - r0_min)*rand(1)
%  r_user = sqrt(r0^2 + (nn*d).^2 - 2*r0*nn*d*sind(theta0));
%  x = exp(-1j*2*pi*fc*(r_user  - r0)/c)/sqrt(N);
%  x = channel_norm(x);
%  r_maxt = RD/100 :1:RD/10;
%  A = (2*randi(2,128,N) - 3)/sqrt(N);
%  noise = sqrt(sigma2)*(randn(128,1)+1i*randn(128,1))/sqrt(2); 
%  Error_vector = [];
%  Angle_Estimates = [];
%  for k=1:length(r_maxt)
%  % Fixed Distance Angle Codebook
% theta_n = linspace(-90,90,512);
% dict = [];
% r_max = r_maxt(k); 
% for ii=1:length(theta_n)
% theta_l = theta_n(ii);    
% r = sqrt(r_max^2 + (nn*d).^2 - 2*r_max*nn*d*sind(theta_l));
% at = exp(-1j*2*pi*fc*(r - r_max)/c)/sqrt(N);
% dict = [dict ; at];
% end
% Y1 = A*x.'+ noise;
% phii =  A*dict.';
% [X, support] = OMP( Y1,phii, 2, 512, 1 );
% angle_estimate = mean((theta_n(support)));
% bt = near_field_manifold( N, d, fc, r_max, deg2rad(mean(angle_estimate)));
% error =A*bt.';
% error = abs(error'*Y1);
% Error_vector =  [Error_vector error];
% Angle_Estimates = [Angle_Estimates angle_estimate];  
%  end
% figure
% plot(r_maxt,Error_vector,'-o','Linewidth',1.2,'markersize',3,'color','blue','MarkerFaceColor','w')
% xlabel('Range','Interpreter','Latex')
% ylabel('Correlation','Interpreter','Latex')
% hold on
% xline(r0,'--r',{'True Range'});
% hold off
% grid on
% set(gca,'fontsize',18);
% box off
% [a, b] = max(Error_vector);
% corr_angle = Angle_Estimates(b)
% corr_range = r_maxt(b)

%% MultiChannel Angle / Range Estimation
 
% Channel Model Near Field
 [H, hc, r0, theta, G] = near_field_channel(N, 1, L, d, fc, fs, M, Rmin, Rmax,sector,1);
 disp(rad2deg(theta))
 H = channel_norm(H);
 r_maxt = [RD/100 RD/50 RD/10 RD/2.5 RD 2*RD] ;
 A = (2*randi(2,128,N) - 3)/sqrt(N);
 noise = sqrt(sigma2)*(randn(128,1)+1i*randn(128,1))/sqrt(2); 
 Error_vector = [];
 Angle_Estimates = [];
 Y1 = A*H.'+ noise;

    % Quacode
    rho = 3;
    beta = 1.2;
    rho_max = 64;
    [Polarcodebook,label] = PolarCodeBook(N, 2, d, lambda, beta, rho,rho_max);
    S = size(Polarcodebook, 2);
    [Haf_hat3,sup3] = SOMP(Y1, A*Polarcodebook,L, S, M);
    

 for k=1:length(r_maxt)
 % Fixed Distance Angle Codebook
theta_n = linspace(-90,90,512);
dict = [];
r_max = r_maxt(k); 
for ii=1:length(theta_n)
theta_l = theta_n(ii);    
at = near_field_manifold( N, d, fc, r_max, deg2rad(theta_l));
dict = [dict ; at];
end
phii =  A*dict.';
[X, support] = OMP( Y1,phii, L, 512, 1 );
angle_estimate_t = theta_n(support);
bt = [];
for bb = 1:length(angle_estimate_t)
angle_estimate = angle_estimate_t(bb);    
btt = near_field_manifold( N, d, fc, r_max, deg2rad(angle_estimate));
bt = [bt  Y1'*(A*btt.')];
end
%error =A*bt.';
error = sum(abs(bt));
Error_vector = [ Error_vector error] ;
Angle_Estimates = [Angle_Estimates ; angle_estimate_t];  
 end
Error_vector_plot(iii,:) = Error_vector;
end

figure
plot(r_maxt,Error_vector_plot(1,:),'-o','Linewidth',2,'markersize',3,'color','magenta','MarkerFaceColor','w')
hold on
plot(r_maxt,Error_vector_plot(2,:),'-v','Linewidth',2,'markersize',3,'color','cyan','MarkerFaceColor','w')
xlabel('Range [m]','Interpreter','Latex')
ylabel('Correlation','Interpreter','Latex')
% xline(r0,'--r',{'True Range'});
xline(RD,'--b',{'Rayleigh distance'});
xline(RD/10,'--b',{'Rayleigh distance /10'});
grid on
set(gca,'fontsize',18);
box on
legend('Near-Field User at 34 m','Far-Field User at 600 m', 'FontSize',8,'Location','best');
[a, b] = max(Error_vector);
corr_angle = Angle_Estimates(b,:);
Dai_Angle = rad2deg(label(1,sup3));
corr_range = r_maxt(b);
Dai_Range = label(2,sup3);
ax = gca;
exportgraphics(ax,'Segregation_of_users.eps','ContentType','vector')
savefig('Segregation_of_users.fig')




% %% Angle Estimation Wrt to Angles
%  theta0_tot = linspace(0,30,5);
%  r0_min = 2*D;
%  Num_of_it = 1;
%  RMSE_Ranges = [];
%  RMSE_Angles = [];
% % Channel Model Near Field
%  for j = 1:length(theta0_tot)
%  theta0 = theta0_tot(j);  
%  r0_max = (cosd(theta0)).^2 * RD/10; 
%  r0_tot = r0_min + (r0_max - r0_min)*rand(Num_of_it);
%  R_E = [];
%  A_E = [];
%  for it = 1: length(r0_tot)
%  r0 =  r0_tot(it);   
%  r_user = sqrt(r0^2 + (nn*d).^2 - 2*r0*nn*d*sind(theta0));
%  x = exp(-1j*2*pi*fc*(r_user  - r0)/c)/sqrt(N);
%  r_max = RD/10;
%  % Fixed Distance Angle Codebook
%  for k= 1:3
% theta_n = linspace(-90,90,512);
% dict = [];
% for ii=1:length(theta_n)
% theta_l = theta_n(ii);    
% r = sqrt(r_max^2 + (nn*d).^2 - 2*r_max*nn*d*sind(theta_l));
% at = exp(-1j*2*pi*fc*(r - r_max)/c)/sqrt(N);
% dict = [dict ; at];
% end
%  y = x*dict';
%  % figure
%  % plot(theta_n,abs(y));
% A = (2*randi(2,128,N) - 3)/sqrt(N);
% Y1 = A*x.';
% phii =  A*dict.';
% [X, support] = OMP( Y1,phii, 1, 512, 1 );
% angle_estimate = (theta_n(support));
% [dict2, label2] = PCodeBook2(N, 2*D, RD, deg2rad(angle_estimate),lambda);
% S2 = size(label2, 2);
% [X2, A2] = OMP( Y1,  A*dict2, 1, S2, 1 );
% Range_Estimate = label2(2,A2);
% r_max = Range_Estimate;
%  end
% Range_Error = (r0 - Range_Estimate);
% Angle_Error = angle_estimate -  theta0;
% R_E = [R_E Range_Error ];
% A_E = [A_E Angle_Error ];
%  end
% Average_R_E = mean(R_E);
% Average_A_E = mean(A_E);
% RMSE_Ranges = [RMSE_Ranges Average_R_E ];
% RMSE_Angles = [RMSE_Angles Average_A_E ];
%  end
%  figure
%  plot(theta0_tot,RMSE_Angles,LineWidth=2)
%  xlabel('Angle (deg)','Interpreter','Latex')
%  ylabel('Angle Estimation error (deg)','Interpreter','Latex')
%  grid on
%  set(gca,'fontsize',18);
% 
%  figure
%  plot(theta0_tot,RMSE_Ranges,LineWidth=2)
%   xlabel('Angle (deg)','Interpreter','Latex')
%  ylabel('Range Estimation error (meters)','Interpreter','Latex')
%  grid on
%  set(gca,'fontsize',18);





