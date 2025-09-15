function [H, hc, r0, theta0, G] = near_field_channel(Nt, K, L, d, fc, B, M, Rmin, Rmax, sector, random)
% K: user num
% L: path num
% Nt: antenna num

r0 = rand(K, L) * (Rmax - Rmin) + Rmin;
% sector = 1*pi/2;
if random
    
    theta0 = rand(K, L) * sector - sector/2 ;
else 
    theta_list = [pi/18, pi/8, -pi/8, pi/12, -pi/16, pi/6];
    theta0 = ones(K, L).*theta_list(1:L);
end
% 
% loc = zeros(K, L, 2);
% loc(:,:,1) = r0 .* cos(theta0);
% loc(:,:,2) = r0 .* sin(theta0);


H = zeros( K,  Nt, M+1);
nn = -(Nt-1)/2:1:(Nt-1)/2;
ssf = (randn(K,L) + 1j*randn(K,L))/sqrt(2);
G = zeros(K, L, M);

c = 3e8;

for k = 1:K
    for m = 1:M+1
       if m == M+1
            f = fc;
       else
            f=fc+B/(M)*(m-1-(M-1)/2);
       end
       
       for l = 1:L
           at = near_field_manifold( Nt, d, fc, r0(k,l), theta0(k,l) );
           g = ssf(k,l) / f;
           H(k, :, m) = squeeze(H(k, :, m)) + g * exp(-1j*2*pi*f*r0(k,l)/c) * at;
           if m <= M
                G(k, l, m) = g * exp(-1j*2*pi*f*r0(k,l)/c);
           end
       end 
    end
end
hc = H(:,:,M+1);
H = H(:,:,1:M);
end

