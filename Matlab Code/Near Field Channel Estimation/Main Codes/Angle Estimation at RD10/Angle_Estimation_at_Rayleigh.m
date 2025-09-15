close all
clc
SNR_dB = 0;  SNR_linear=10.^(SNR_dB/10.);
sigma2=1/SNR_linear;
N_iter = 100; 
%%% system parameters
N = 256; % number of beams (transmit antennas)
K = 1; % number of users
N_RF = 4; % number of RF chains
M = 1; % number of subcarriers
L = 1; % number of paths per user


fc = 28e9; % carrier frequency
sector = pi/3;
fs = 100e6; % bandwidth
Q = 32;  % number of pilot blocks
tmax = 20e-9; % maximum of the path delay
f = zeros(1,M);
c = 3e8;
lambda_c = c/fc;
d = lambda_c / 2;
A = (N-1).*d;
RD = (2*A.^2)/lambda_c;
R = 3:30:500;
NMSE = zeros(1,length(R));
NMSE0 = zeros(1,length(R));
NMSE1 = zeros(1,length(R));
NMSE2 = zeros(1,length(R));
NMSE3 = zeros(1,length(R));
NMSE4 = zeros(1,length(R));


%DFT    
s = 2;
D = s*N; %字典规模
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


R_len = length(R);
t0 = clock;
for i_r = 1:R_len
    Rmin = R(i_r);
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
            %% Far field on-grid
            [Haf_hat1, sup1] = SOMP(Z, Phi*DFT, L, D, M);
            angle_estimate_farfield = asin(col(sup1));       
            %% Near field on-grid
            [Haf_hat3,sup3] = SOMP(Z, Phi*Polarcodebook,L, S, M);
            Dai_angle = (label(1,sup3));
           %%    Near Field Beam Focus
           [X1, A1] = SOMP( Z,  Phi*dict1, L, S1, M );
           angle_estimateBF = label1(1,A1);


            error0 = error0 + sum(abs(theta(k,:) - flip(angle_estimate_farfield)));
            error1 = error1 + sum(abs(theta(k,:) - flip(Dai_angle)));
            error2 = error2 + sum(abs(theta(k,:) - flip(angle_estimateBF)));
        end     
    end
    NMSE0(i_r) = error0/K/N_iter;
    NMSE1(i_r) = error1/K/N_iter;
    NMSE2(i_r) = error2/K/N_iter;
    
end

x = 1:3:size(R,2);
figure; hold on; box on; grid on;
plot(R(x),pow2db(NMSE0(x)),'-o','Linewidth',1.2,'markersize',5,'color','cyan','MarkerFaceColor','w');
plot(R(x),pow2db(NMSE1(x)),'-s','Linewidth',1.2,'markersize',5,'color','blue','MarkerFaceColor','w');
plot(R(x),pow2db(NMSE2(x)),'-v','Linewidth',1.2,'markersize',5,'color','green','MarkerFaceColor','w');
xlabel('Distance (meters)','Interpreter','Latex');
ylabel('Angle Error (dB)','Interpreter','Latex');
set(gca,'fontsize',18);
grid on
xline(RD/10,'--r',{'Rayleigh Distance/10'});
xline(RD,'--r',{'Rayleigh Distance'});
legend( 'Far Field','Dai','BF-OMP','FontSize',12,'Location','best');