clear all;close all;clc;
Nt = 256;
Nr = 1;
fc = 28e9;
c = 3e8;
lambda_c = c/fc;
d = lambda_c / 2;
eps = 1e-3;

D = 0.5 * Nt^2 * lambda_c;

num_r = 1000;
num_theta = 1000;
r_start = 1.5;
r_end = 348;


r_list = linspace(r_start, r_end, num_r);
theta_list = linspace(-pi/2, pi/2, num_theta);
x_list = r_list.' * cos(theta_list);
y_list = r_list.' * sin(theta_list);

num_subfig = 1;
g = zeros(num_subfig, num_r, num_theta);

r_center = [20];

for i_subfig = 1:1
    r = r_center(i_subfig);
    theta = -pi/6;
    H = near_field_manifold(Nt, d, fc, r, theta);
    w = narrow_focus(H,Nt);
    for idx_r = 1:num_r
        for idx_theta = 1:num_theta
            H = near_field_manifold(Nt, d, fc, r_list(idx_r), theta_list(idx_theta));
            g(i_subfig, idx_r, idx_theta)=  abs(H*w);
        end
    end
end




%% fig 1
figure;
hold on
meshc(x_list,y_list,(squeeze(g(1,:,:))));
grid('off');
colorbar
colormap('parula')
xlabel('x-axis [m]','Interpreter','Latex')
ylabel('y axis [m]','Interpreter','Latex')
set(gca,'fontsize',18);
