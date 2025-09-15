clear all
close all
clc

N_iter = 10; 
%% system parameters
N1 = 128:32:256; % number of beams (transmit antennas)
K = 4; % number of users
N_RF = 4; % number of RF chains
M = 64; % number of subcarriers
L = 3; % number of paths per user
fc = 28e9; % carrier frequency
sector = pi/3; % User Coverage/2
fs = 100e6; % bandwidth
Q = 32;  % number of pilot blocks
c = 3e8;  % Speed of light 
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % element spacing
R_len = length(N1);

NMSE = zeros(1,length(N1));
NMSE0 = zeros(1,length(N1));
NMSE1 = zeros(1,length(N1));
NMSE2 = zeros(1,length(N1));
NMSE3 = zeros(1,length(N1));
NMSE4 = zeros(1,length(N1));

for i_r = 1:R_len
N = N1(i_r);
A = (N-1).*d; % Array Aperture
RD = (2*A.^2)/lambda_c; % Rayleigh Distance


SNR_dB = 10;  SNR_linear=10.^(SNR_dB/10.);
sigma2=1./SNR_linear;



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

t0 = clock;

    Rmin = 5;
    Rmax = Rmin;
    
     error0 = 0; error1 = 0; error2 = 0; error3 =0;error4 =0;
    for iter = 1:N_iter
        fprintf('R = %.4f[%d/%d] | iteration:[%d/%d] | run %.4f s\n', Rmin, i_r, R_len, iter, N_iter, etime(clock, t0)); 
        
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
               angle_estimate2 = (angle_estimate1);
           end

           % FD =  RD * range_estimate * (1/(RD-10*range_estimate) - 1/(RD+10*range_estimate));
           if range_estimate > RD/(2*10)   
           r_maxt = range_estimate-5:0.5:range_estimate+5;
           end

           if (range_estimate < RD/(2*10)) && (range_estimate > RD/(50))    
           r_maxt = range_estimate-2.5:0.25:range_estimate+2.5;
           end

           if range_estimate < RD/(50)   
           r_maxt = range_estimate-1:0.15:range_estimate+1;
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

            
            error0 = error0 + norm(Hsf - Hsf_hat0,'fro')^2/norm(Hsf,'fro')^2;
            error1 = error1 + norm(Hsf - Polarcodebook*Haf_hat3,'fro')^2/norm(Hsf,'fro')^2;
            error2 = error2 + norm(Hsf - BF_OMP,'fro')^2/norm(Hsf,'fro')^2;
            error3 = error3 + norm(Hsf - dict1*X1,'fro')^2/norm(Hsf,'fro')^2;
            % error4 = error4 + norm(Hsf - Dai,'fro')^2/norm(Hsf,'fro')^2;
        end     
    end
    NMSE0(i_r) = error0/K/N_iter;
    NMSE1(i_r) = error1/K/N_iter;
    NMSE2(i_r) = error2/K/N_iter;
    NMSE3(i_r) = error3/K/N_iter;
    % NMSE4(i_r) = error4/K/N_iter;
end

x = 1:1:size(N1,2);
figure; hold on; box on; grid on;
plot(N1(x),10*log10(NMSE0(x)),'-o','Linewidth',1.2,'markersize',5,'color','red','MarkerFaceColor','w');
plot(N1(x),10*log10(NMSE1(x)),'-s','Linewidth',1.2,'markersize',5,'color','blue','MarkerFaceColor','w');
plot(N1(x),10*log10(NMSE2(x)),'-v','Linewidth',1.2,'markersize',5,'color','green','MarkerFaceColor','w');
plot(N1(x),10*log10(NMSE3(x)),'-v','Linewidth',1.2,'markersize',5,'color','black','MarkerFaceColor','w');
% plot(R(x),10*log10(NMSE4(x)),'-v','Linewidth',1.2,'markersize',5,'color','cyan','MarkerFaceColor','w');
xlabel('Num of Antenna Element','Interpreter','Latex');
ylabel('NMSE (dB)','Interpreter','Latex');
% ylim([-40, 0]);
set(gca,'fontsize',18);
grid on
legend( 'Oracle LS','P-OMP','Refined BF-OMP','BF-OMP' );
