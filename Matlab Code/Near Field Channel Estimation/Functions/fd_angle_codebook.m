function [dict, label] = fd_angle_codebook(N, r_max,lambda, angle_estimate)

c = 3e8;
fc= 3e8/lambda;
d = lambda/2;
nn = -(N-1)/2:1:(N-1)/2;
angle_estimate = rad2deg(angle_estimate);
theta_n = linspace(-(angle_estimate+5),angle_estimate+5,512);
dict = [];
for ii=1:length(theta_n)
theta_l = theta_n(ii);    
r = sqrt(r_max^2 + (nn*d).^2 - 2*r_max*nn*d*sind(theta_l));
at = exp(-1j*2*pi*fc*(r - r_max)/c)/sqrt(N);
dict = [dict ; at];
end
label = theta_n;
dict = dict.';