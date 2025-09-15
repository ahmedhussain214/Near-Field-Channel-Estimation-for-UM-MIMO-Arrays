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
sector = pi/3; % User Coverage/2
fs = 100e6; % bandwidth
Q = 32;  % number of pilot blocks
c = 3e8;  % Speed of light 
lambda_c = c/fc; % wavelength
d = lambda_c / 2; % element spacing
A = (N-1).*d; % Array Aperture
RD = (2*A.^2)/lambda_c; % Rayleigh Distance


SNR_dB = [ -20 -15 -10 -5 0 5 10 15 20 ];  SNR_linear=10.^(SNR_dB/10.);
sigma1=1./SNR_linear;

ss = 1:1:5*L;
NMSE0 = zeros(length(sigma1),length(ss));

%mycodebook
[dict1, label1] = PCodeBook(N, A, RD, lambda_c);
S1 = size(label1, 2);
[dict3, label3] = PCodeBooka(N, A, RD, lambda_c);
S3 = size(label3, 2);


R_len = length(sigma1);
t0 = clock;
for i_r = 1:R_len

    sigma2 =sigma1(i_r);
    
     c = 1/sigma2*ones(1,Q*N_RF);
     C = diag(c);
      Error = zeros(1,length(ss));
    for st = 1:length(ss)
        s=ss(st);
        error = 0; 
    parfor iter = 1:N_iter
        fprintf('SNR = %.4f[%d/%d] | iteration:[%d/%d] | run %.4f s\n', sigma2, i_r, R_len, iter, N_iter, etime(clock, t0)); 
    Rmin = 3+ rand(1)*(30);
    Rmax = Rmin;
        % Wideband spatial channel
        [H, hc, r, theta, G] = near_field_channel(N, K, L, d, fc, fs, M, Rmin, Rmax,sector,1);
        true_angle = rad2deg (theta);
        H = channel_norm(H);
        Phi = (2*randi(2,Q*N_RF,N) - 3)/sqrt(N);
        noise = sqrt(sigma2)*(randn(Q*N_RF,M)+1i*randn(Q*N_RF,M))/sqrt(2); 
        
            Hsf =  reshape(H(1, :, :), [N, M]);    % Hsf = F*Haf
            
            % adaptive selecting matrix
            Z = Phi*Hsf + noise;
            
            Znorm = norm(Z, 'fro');
         
           %%    Near Field Beam Focus

           [Haf_hat5,sup5] = SOMP( Z,  Phi*dict1, s, S1, M);

            error = error + norm(Hsf - dict1*Haf_hat5,'fro')^2/norm(Hsf,'fro')^2;
            
        
    end
Error(st) = error/N_iter;
    end
    NMSE0(i_r,:) = Error;
   
end
[~ , indices] = min(NMSE0,[],2);

