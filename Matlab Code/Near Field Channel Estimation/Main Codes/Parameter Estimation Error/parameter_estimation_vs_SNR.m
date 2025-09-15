clear all
close all
clc

N_iter = 1000; 
%% system parameters
N = 256; % number of beams (transmit antennas)
K = 1; % number of users
N_RF = 4; % number of RF chains
M = 64; % number of subcarriers
L = 1; % number of paths per user
fc = 28e9; % carrier frequency
sector = pi/6; % User Coverage/2
fs = 100e6; % bandwidth
Q = 32;  % number of pilot blocks
c = 3e8;  % Speed of light 
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % element spacing
A = (N-1).*d; % Array Aperture
RD = (2*A.^2)/lambda_c; % Rayleigh Distance


SNR_dB = [-10 -5 0 5 10 ];  SNR_linear=10.^(SNR_dB/10.);
sigma1=1./SNR_linear;

NMSE0 = zeros(1,length(sigma1));
NMSE1 = zeros(1,length(sigma1));
NMSE2 = zeros(1,length(sigma1));

NMSE3 = zeros(1,length(sigma1));
NMSE4 = zeros(1,length(sigma1));
NMSE5 = zeros(1,length(sigma1));

%DFT    
s = 2;
D = s*N; 
row = (-(N - 1)/2:(N - 1)/2)' ;
col = -1 + 2/D : 2/D : 1 ;
DFT  =  exp( 1j*  pi * row * col ) / sqrt(N);

% Quacode
rho = 3;
beta = 1.2;
rho_max = 34;
[Polarcodebook,label] = PolarCodeBook(N, s, d, lambda_c, beta, rho,rho_max);
S = size(Polarcodebook, 2);

%mycodebook
[dict1, label1] = PCodeBook(N, A, RD, lambda_c);
S1 = size(label1, 2);



R_len = length(sigma1);
t0 = clock;
for i_r = 1:R_len

    sigma2 =sigma1(i_r);
     error0 = 0; error1 = 0; error2 = 0; error3 =0;error4 =0; error5 =0;
    for iter = 1:N_iter
        fprintf('R = %.4f[%d/%d] | iteration:[%d/%d] | run %.4f s\n', 0, i_r, R_len, iter, N_iter, etime(clock, t0)); 
        
    Rmin = 3+5*rand(1);
    Rmax = Rmin;

        % Wideband spatial channel
        [H, hc, r, theta, G] = near_field_channel(N, K, L, d, fc, fs, M, Rmin, Rmax,sector,0);
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
                      
            %% Near field on-grid
            [Haf_hat3,sup3] = SOMP(Z, Phi*Polarcodebook,L, S, M);
            Dai_Range = label(2,sup3);
            Dai_angle = (label(1,sup3));
           %%    Near Field Beam Focus
           [X1, A1] = SOMP( Z,  Phi*dict1, L, S1, M );

           angle_estimate1 = label1(1,A1);
           rad2deg (angle_estimate1);
           range_estimate= mean(label1(2,A1));

           if length(angle_estimate1) > 1

               angle_estimate2 = max(abs(angle_estimate1));
           else
               angle_estimate2 = abs(angle_estimate1);
           end

           % FD =  RD * range_estimate * (1/(RD-10*range_estimate) - 1/(RD+10*range_estimate));
           if range_estimate > RD/(2*10)   
           r_maxt = range_estimate-5:0.5:range_estimate+5;
           end

           if (range_estimate < RD/(2*10)) && (range_estimate > RD/(50))    
           r_maxt = range_estimate-2.5:0.25:range_estimate+2.5;
           end

           if range_estimate < RD/(50)   
           r_maxt = range_estimate-1:0.1:range_estimate+1;
           end

           Error_vector = [];
           Angle_Estimates = [];
           Supports = [];
           Gs = [];
            for kk=1:length(r_maxt)
             % Fixed Distance Angle Codebook
            r_max = r_maxt(kk);  
            [dict2, label2] = fd_angle_codebook(N, r_max,lambda_c, angle_estimate2);
            S2 = size(label2, 2);
            [X, support] = SOMP( Z,Phi*dict2, L, S2, M );
            angle_estimate_t = label2(support);
            Supports = [Supports;support];
            Gs = [Gs ; transpose(X(support)) ];  
            bt = [];
            for bb = 1:length(angle_estimate_t)
            angle_estimate = angle_estimate_t(bb);    
            btt = near_field_manifold( N, d, fc, r_max, deg2rad(angle_estimate));
            bt = [bt  Z'*(Phi*btt.')];
            end
            %error =A*bt.';
            error = sum(abs(bt),'all');
            Error_vector = [ Error_vector error] ;
            Angle_Estimates = [Angle_Estimates ; angle_estimate_t];
            end

            [a, b] = max(Error_vector);
            corr_angle = Angle_Estimates(b,:);
            corr_range = r_maxt(b);
            % sup = Supports(b,:);
            % g = Gs(b,:);


            BF_OMP = Oracle_LS(Z, Phi, N, M,corr_range.*ones(1,L), deg2rad(corr_angle), d, fc);
            % Dai = Oracle_LS2(Z, Phi, N, M,corr_range.*ones(1,L), deg2rad(corr_angle), d, fc,g,sup);
            % Dai = Oracle_LS2(Z, Phi, N, M,Dai_Range, asin(Dai_angle), d, fc,Haf_hat3,sup3);
            % Dai = Oracle_LS(Z, Phi, N, M,Dai_Range, asin(Dai_angle), d, fc);

           % [dict2, label2] = PCodeBook2(N, A, RD, angle_estimate,lambda_c, range_estimate);
           % S2 = size(label2, 2);
           % [X2, A2] = SOMP( Z,  Phi*dict2, L, S2, M );
           corr_angle = deg2rad(corr_angle);
            mean(Dai_Range);

            error0 = error0 + 1/r(1,1).*(abs(r(1,1) - mean(Dai_Range)));
            error1 = error1 + 1/r(1,1).*(abs(r(1,1) - range_estimate));
            error2 = error2 + 1/r(1,1).*(abs(r(1,1) - corr_range));


            error3 = error3 + sum(abs(theta(k,:) - flip(Dai_angle)));
            error4 = error4 + sum(abs(theta(k,:) - flip(angle_estimate1)));
            error5 = error5 + sum(abs(theta(k,:) - flip(corr_angle)));

        end     
    end
    NMSE0(i_r) = error0/K/N_iter;
    NMSE1(i_r) = error1/K/N_iter;
    NMSE2(i_r) = error2/K/N_iter;

    NMSE3(i_r) = error3/K/N_iter;
    NMSE4(i_r) = error4/K/N_iter;
    NMSE5(i_r) = error5/K/N_iter; 
end

x = 1:1:size(SNR_dB,2);
figure; hold on; box on; grid on;
plot(SNR_dB(x),(NMSE3(x)),'-o','Linewidth',1.2,'markersize',5,'color','red','MarkerFaceColor','w');
plot(SNR_dB(x),(NMSE4(x)),'-s','Linewidth',1.2,'markersize',5,'color','blue','MarkerFaceColor','w');
plot(SNR_dB(x),(NMSE5(x)),'-v','Linewidth',1.2,'markersize',5,'color','green','MarkerFaceColor','w');

xlabel('SNR (dB)','Interpreter','Latex');
ylabel('Mean Angle Estimation error','Interpreter','Latex');
set(gca,'fontsize',18);
grid on
legend( 'P-OMP','BF-OMP','Refined BF-OMP','FontSize',12,'Location','best'  );

x = 1:1:size(SNR_dB,2);
figure; hold on; box on; grid on;
plot(SNR_dB(x),(NMSE0(x)),'-o','Linewidth',1.2,'markersize',5,'color','red','MarkerFaceColor','w');
plot(SNR_dB(x),(NMSE1(x)),'-s','Linewidth',1.2,'markersize',5,'color','blue','MarkerFaceColor','w');
plot(SNR_dB(x),(NMSE2(x)),'-v','Linewidth',1.2,'markersize',5,'color','green','MarkerFaceColor','w');
xlabel('SNR (dB)','Interpreter','Latex');
ylabel('Mean Distance Estimation error','Interpreter','Latex');
set(gca,'fontsize',18);
grid on
legend( 'P-OMP','BF-OMP','Refined BF-OMP','FontSize',12,'Location','best'  );
