clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Codebook Size vs Number of antenna elements %
%% system parameters
NN1 = 64:128:512; % number of beams (transmit antennas)
S_Dai = [];
S_BF =[];
for i = 1:length(NN1)
N = NN1(i);    
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
%% Polar codebook (Dai)
s = 1;
rho = 2*A;
beta = 1.2;
rho_max = RD;
[Polarcodebook,label] = PolarCodeBook(N, s, d, lambda_c, beta, rho,rho_max);
S = size(Polarcodebook, 2);
%% BF Polar Codebook
[dict1, label1] = PCodeBook(N, A, RD, lambda_c);
S1 = size(label1, 2);
S2 = S1/3;
%%
S_Dai = [S_Dai S]; 
S_BF = [S_BF S2];
end
figure; hold on; box on; grid on;
plot(NN1,S_Dai,'-o','Linewidth',2,'markersize',5,'color','blue','MarkerFaceColor','w');
plot(NN1,S_BF,'-^','Linewidth',2,'markersize',5,'color','red','MarkerFaceColor','w');
xlabel('Number of Antenna Elements','Interpreter','Latex');
ylabel('Codebook Size','Interpreter','Latex');
set(gca,'fontsize',18);
grid on
legend( 'Polar Codebook [16]','Proposed Polar Codebook','FontSize',12,'Location','Northwest');
% ax = gca;
% exportgraphics(ax,'Polarcodebooksize_vs_Num_antenna_elements.eps','ContentType','vector')
% exportgraphics(ax,'Polarcodebooksize_vs_Num_antenna_elements.jpg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Codebook Size vs Frequency %
%% system parameters
NN2 = 20e9:10e9:100e9; % number of beams (transmit antennas)
S_Dai_1 = [];
S_BF_1 =[];
for i = 1:length(NN2)
N = 256;    
N_RF = 4; % number of RF chains
M = 64; % number of subcarriers
L = 3; % number of paths per user
fc = NN2(i); % carrier frequency
sector = pi/3; % User Coverage/2
fs = 100e6; % bandwidth
Q = 32;  % number of pilot blocks
c = 3e8;  % Speed of light 
lambda_c = c/fc; % wavelength
A = 1; % Fix Array Aperture
d = A/(N-1);% element spacing
RD = (2*A.^2)/lambda_c; % Rayleigh Distance
%% Polar codebook (Dai)
s = 1;
rho = 2*A;
beta = 1.2;
rho_max = RD;
[Polarcodebook,label] = PolarCodeBook(N, s, d, lambda_c, beta, rho,rho_max);
S = size(Polarcodebook, 2);
%% BF Polar Codebook
[dict1, label1] = PCodeBook(N, A, RD, lambda_c);
S1 = size(label1, 2);
S2 = S1/3;
%%
S_Dai_1 = [S_Dai_1 S]; 
S_BF_1 = [S_BF_1 S2];
end

figure; hold on; box on; grid on;
plot(NN2/1e9,S_Dai_1,'-o','Linewidth',2,'markersize',5,'color','blue','MarkerFaceColor','w');
plot(NN2/1e9,S_BF_1,'-^','Linewidth',2,'markersize',5,'color','red','MarkerFaceColor','w');
xlabel('Frequency (GHz)','Interpreter','Latex');
ylabel('Codebook Size','Interpreter','Latex');
set(gca,'fontsize',18);
grid on
legend( 'Polar Codebook [16]','Proposed Polar Codebook','FontSize',12,'Location','Northwest');
% ax = gca;
% exportgraphics(ax,'Polarcodebooksize_vs_frequency.eps','ContentType','vector')
% exportgraphics(ax,'Polarcodebooksize_vs_frequency.jpg')

%%%% Single plot






% Create the first axis
figure; hold on; box on; grid on;
yyaxis left;

plot(S_Dai, NN1 , '--square', 'LineWidth', 1.2, 'MarkerSize', 5, 'MarkerFaceColor', 'w');
hold on
plot( S_BF,NN1, '-*', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'w');
xlabel( 'Codebook Size', 'Interpreter', 'Latex');
ylabel( 'Number of Antenna elements', 'Interpreter', 'Latex');



yyaxis right;

plot(S_Dai_1,NN2/1e9,  '--square', 'LineWidth', 1.2, 'MarkerSize', 5, 'MarkerFaceColor', 'w');
hold on
plot(S_BF_1,NN2/1e9,  '-*', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'w');
xlabel( 'Codebook Size', 'Interpreter', 'Latex');
ylabel( 'Frequency (GHz)', 'Interpreter', 'Latex');
legend('Polar Codebook [16]', 'Proposed Polar Codebook', 'Polar Codebook [16]', 'Proposed Polar Codebook','FontSize', 12, 'Location', 'Southeast','Interpreter', 'Latex');
xlim([64 12000])
set(gca,'fontsize',18);



