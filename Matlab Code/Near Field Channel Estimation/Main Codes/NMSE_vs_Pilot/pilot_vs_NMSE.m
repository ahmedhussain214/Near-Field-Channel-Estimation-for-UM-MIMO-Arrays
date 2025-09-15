clear all
close all
clc

N_iter = 250; 
%% system parameters
N = 256; % number of beams (transmit antennas)
K = 4; % number of users
N_RF = 4; % number of RF chains
M = 64; % number of subcarriers
L = 3; % number of paths per user
fc = 28e9; % carrier frequency
sector = pi/6; % User Coverage/2
fs = 100e6; % bandwidth
Q1 = 16:8:72;  % number of pilot blocks
c = 3e8;  % Speed of light 
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % element spacing
A = (N-1).*d; % Array Aperture
RD = (2*A.^2)/lambda_c; % Rayleigh Distance


SNR_dB = 10;  SNR_linear=10.^(SNR_dB/10.);
sigma2=1./SNR_linear;

NMSE = zeros(1,length(Q1));
NMSE0 = zeros(1,length(Q1));
NMSE1 = zeros(1,length(Q1));
NMSE2 = zeros(1,length(Q1));
NMSE3 = zeros(1,length(Q1));
NMSE4 = zeros(1,length(Q1));
NMSE5 = zeros(1,length(Q1));
NMSE6 = zeros(1,length(Q1));

%% DFT    
s = 2;
D = s*N; 
row = (-(N - 1)/2:(N - 1)/2)' ;
col = -1 + 2/D : 2/D : 1 ;
DFT  =  exp( 1j*  pi * row * col ) / sqrt(N);

%% Linlong codebook
rho = 2*A;
beta = 1.2;
rho_max = 348;
[Polarcodebook,label] = PolarCodeBook(N, s, d, lambda_c, beta, rho,rho_max);
S = size(Polarcodebook, 2);

%% BF-OMP codebook
[dict1, label1] = PCodeBook(N, A, RD, lambda_c);
S1 = size(label1, 2);
%% Refined BF-OMP codebook
[dict3, label3] = PCodeBooka(N, A, RD, lambda_c);
S3 = size(label3, 2);


R_len = length(Q1);
t0 = clock;
for i_r = 1:R_len
    Q =Q1(i_r);
     error0 = 0; error1 = 0; error2 = 0; error3 =0;error4 =0;error5 =0; error6 =0;
    parfor iter = 1:N_iter
        fprintf('Q = %.4f[%d/%d] | iteration:[%d/%d] | run %.4f s\n', Q, i_r, R_len, iter, N_iter, etime(clock, t0)); 
    Rmin = 34+ rand(1)*(215);
    Rmax = Rmin;
        % Wideband spatial channel
        [H, hc, r, theta, G] = near_field_channel(N, K, L, d, fc, fs, M, Rmin, Rmax,sector,1);
        true_angle = rad2deg (theta);
        H = channel_norm(H);
        Phi = (2*randi(2,Q*N_RF,N) - 3)/sqrt(N);
        noise = sqrt(sigma2)*(randn(Q*N_RF,M)+1i*randn(Q*N_RF,M))/sqrt(2);
       
        for k = 1 : K
            Hsf =  reshape(H(k, :, :), [N, M]);    % Hsf = F*Haf  
            % adaptive selecting matrix
            Z = Phi*Hsf + noise;
            Znorm = norm(Z, 'fro');
            %% Oracle LS
             Hsf_hat0 = Oracle_LS(Z, Phi, N, M, r(k,:), theta(k,:), d, fc);
             %% Far field on-grid
            [Haf_hat1, sup1] = SOMP(Z, Phi*DFT, 2*L, D, M);          
            %% Near field on-grid
            [Haf_hat3,sup3] = SOMP(Z, Phi*Polarcodebook,2*L, S, M);
           %%    Near Field Beam Focus
           [Haf_hat2, sup2] = SOMP_N( Z,  Phi*dict1, 2*L, S1, M,sigma2);
           %%    Refined  Beam Focus
           [Haf_hat4,sup4,] = SOMP_N( Z,  Phi*dict3, 2*L, S3, M,sigma2);
            %% SIGW P-SOMP
            [A4, G4, r4, theta4, Reslist2] =P_SIGW( Z/Znorm, Phi'/Znorm, Polarcodebook(:, sup3), Haf_hat3(sup3,:), label(2, sup3), label(1, sup3), ...
            fc, d, L, 3, 20, 10, 1, 1e-8);
            Hsf_hat4 = A4*G4;
              %% SIGW BF-OMP
            [A6, G6, r6, theta6, Reslist1] =P_SIGW_N( Z/Znorm, Phi'/Znorm, dict1(:, sup2), Haf_hat2(sup2,:), label1(2, sup2), label1(1, sup2), ...
            fc, d,4,1000/Q, 30, 10, 1, 1e-8);
            Hsf_hat5 = A6*G6;

            
            error0 = error0 + norm(Hsf - Hsf_hat0,'fro')^2/norm(Hsf,'fro')^2;
            error1 = error1 + norm(Hsf - Polarcodebook*Haf_hat3,'fro')^2/norm(Hsf,'fro')^2;
            error2 = error2 + norm(Hsf - dict1*Haf_hat2,'fro')^2/norm(Hsf,'fro')^2;
            error3 = error3 + norm(Hsf - dict3*Haf_hat4,'fro')^2/norm(Hsf,'fro')^2;
            error4 = error4 + norm(Hsf - DFT*Haf_hat1,'fro')^2/norm(Hsf,'fro')^2;
            error5 = error5 + norm(Hsf - Hsf_hat4,'fro')^2/norm(Hsf,'fro')^2;
            error6 = error6 + norm(Hsf - Hsf_hat5,'fro')^2/norm(Hsf,'fro')^2;
        end     
    end
    NMSE0(i_r) = error0/K/N_iter;
    NMSE1(i_r) = error1/K/N_iter;
    NMSE2(i_r) = error2/K/N_iter;
    NMSE3(i_r) = error3/K/N_iter;
    NMSE4(i_r) = error4/K/N_iter;
    NMSE5(i_r) = error5/K/N_iter;
    NMSE6(i_r) = error6/K/N_iter;
end

x = 1:1:size(Q1,2);
figure; hold on; box on; grid on;

plot(Q1(x),10*log10(NMSE4(x)),'-v','Linewidth',2,'markersize',5,'color','[0.9290 0.6940 0.1250]','MarkerFaceColor','w');
plot(Q1(x),10*log10(NMSE1(x)),'-s','Linewidth',2,'markersize',5,'color','[0 0.4470 0.7410]','MarkerFaceColor','w');
plot(Q1(x),10*log10(NMSE2(x)),'-v','Linewidth',2,'markersize',5,'color','[0.6350 0.0780 0.1840]','MarkerFaceColor','w');
plot(Q1(x),10*log10(NMSE3(x)),'-o','Linewidth',2,'markersize',5,'color','[0.4660 0.6740 0.1880]','MarkerFaceColor','w');
plot(Q1(x),10*log10(NMSE5(x)),'--','Linewidth',2,'markersize',5,'color','blue','MarkerFaceColor','w');
plot(Q1(x),10*log10(NMSE6(x)),'--','Linewidth',2,'markersize',5,'color','green','MarkerFaceColor','w');
plot(Q1(x),10*log10(NMSE0(x)),'--','Linewidth',2,'markersize',5,'color','black','MarkerFaceColor','w');

xlabel('Pilot length ($Q$)','Interpreter','Latex');
ylabel('NMSE (dB)','Interpreter','Latex');
% xlim([Q1(1),Q1(end)])
grid on
legend('FF-SOMP[13]','P-SOMP[16]','BF-SOMP','RF-BF-SOMP',' SIGW P-SOMP[16]','SIGW BF-SOMP','Oracle LS', 'FontSize',8,'Location','best');
set(gca,'fontsize',18);
% print('NMSE12.pdf','-dpdf')
ax = gca;
exportgraphics(ax,'Pilot_vs_NMSE_FF.eps','ContentType','vector')
exportgraphics(ax,'Pilot_vs_NMSE_FF   .jpg')
exportgraphics(ax,'Pilot_vs_NMSE_fF.eps','ContentType','vector')
savefig('Pilot_vs_NMSE_FF.fig')
