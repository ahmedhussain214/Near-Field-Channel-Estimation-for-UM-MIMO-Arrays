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
T_Rmin = [ RD/80 RD/60 RD/40 RD/20  RD/5  100*RD]

F = RD/100;
rr =[];
BD=1;
while (F < 1000*RD) && (BD>0)
rr = [rr F ];
theta0 = 0;
BD =  RD*cosd(theta0).^2 * F * (1/(RD*cosd(theta0).^2-10*F) - 1/(RD*cosd(theta0).^2+10*F));
F = F+5*BD;
end


for iii = 1:length(T_Rmin)
Rmin = T_Rmin(iii);   
Rmax = Rmin;
 
%% MultiChannel Angle / Range Estimation
 
% Channel Model Near Field
 [H, hc, r0, theta, G] = near_field_channel(N, 1, L, d, fc, fs, M, Rmin, Rmax,sector,1);
 disp(rad2deg(theta))
 H = channel_norm(H);
 %r_maxt = [RD/100 RD/50 RD/20 RD/9 RD/5 10*RD ] ;
 r_maxt = [rr 10*RD ];
 %r_maxt = linspace(RD/80, RD, 1024);
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
Error_vector_plot(iii,:) = Error_vector/max(Error_vector);
end
set(groot,'defaultAxesTickLabelInterpreter','latex');
C = linspecer(6);
figure
plot(r_maxt,Error_vector_plot(1,:),'-o','Linewidth',2,'markersize',3,'color',C(1, :),'MarkerFaceColor','w')
hold on
plot(r_maxt,Error_vector_plot(2,:),'-o','Linewidth',2,'markersize',3,'color',C(2, :),'MarkerFaceColor','w')
hold on
plot(r_maxt,Error_vector_plot(3,:),'-o','Linewidth',2,'markersize',3,'color',C(3, :),'MarkerFaceColor','w')
hold on
plot(r_maxt,Error_vector_plot(4,:),'-o','Linewidth',2,'markersize',3,'color',C(4, :),'MarkerFaceColor','w')
hold on
plot(r_maxt,Error_vector_plot(5,:),'-o','Linewidth',2,'markersize',3,'color',C(5, :),'MarkerFaceColor','w')
hold on
plot(r_maxt,Error_vector_plot(6,:),'-o','Linewidth',2,'markersize',3,'color',C(6, :),'MarkerFaceColor','w')
set(gca,'XScale','log');
xlabel('Distance (meters)','Interpreter','Latex')
ylabel('Normalized Correlation','Interpreter','Latex')
% xline(r0,'--r',{'True Range'});
xline(RD,'--b',{'Rayleigh distance'});
xline(RD/10,'--b',{'Rayleigh distance /10'});
grid on
box on
legend('Near-Field User at .012$r_{RD}$','Near-Field User at 0.016$r_{RD}$','Near-Field User at .02$r_{RD}$', ...
    'Near-Field User at .05$r_{RD}$','Near-Field User at 0.2$r_{RD}$','Far-Field User at 10$r_{RD}$', ...
    'FontSize',8,'Location','best','Interpreter','Latex');
xticks([0.01*RD 0.1*RD 1*RD 10*RD ])
%xticklabels(r_maxt)
xticklabels({'.01$r_{RD}$','0.1$r_{RD}$','1$r_{RD}$','10$r_{RD}$'})
xlim([2*D 10*RD])
set(gca,'fontsize',18);





