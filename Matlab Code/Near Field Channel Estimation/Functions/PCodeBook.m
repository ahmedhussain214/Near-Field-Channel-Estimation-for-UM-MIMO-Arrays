function [dict, label] = PCodeBook(N, D, RD, lambda)
d = lambda/2;
fc= 3e8/lambda;
theta0 = linspace(-pi/2,pi/2, 3*N);

 % D = 1 * N;
 % theta0 = -1 + 2/N : 2/N : 1;


label = [];
for  i =1:length(theta0)
F = 2*D;
theta = theta0(i);
rr= [];

while F <= (RD/10)*(cos(theta).^2)
rr = [F rr];

% r =  RD*  F * (1/(RD-10*F) - 1/(RD+10*F));
r =  RD*cos(theta).^2 * F * (1/(RD*cos(theta).^2-10*F) - 1/(RD*cos(theta).^2+10*F));

F = F+r/2;

end
t =(RD/10)*(cos(theta).^2);
if t>2*D
    v = [t RD];
else
    v =RD;
end
rr =[rr   v];
label1 = [ones(1,length(rr))*(theta) ; rr ];
label = [label label1];
end
dict = zeros(N,size(label,2));
for i= 1:size(label,2)
Nt = N;
theta = label(1,i);
r = label(2,i);
% at = near_field_manifold( Nt, d, fc, r, theta );
at = polar_domain_manifold( Nt, d, fc, r, theta );
dict(:,i) =at;
end
end
% 















% function [dict, label] = PCodeBook(N, D, RD, lambda)
% d = lambda/2;
% fc= 3e8/lambda;
% % theta0 = linspace(-pi/2,pi/2, N);
% 
%   S = 3 * N;
%  theta0 = -1 + 2/S : 2/S : 1;
% 
% 
% label = [];
% for  i =1:length(theta0)
% F = 2*D;
% theta = asin(theta0(i));
% rr= [];
% 
% while F <= (RD/10)*(cos(theta).^2)
% rr = [F rr];
% 
% % r =  RD*  F * (1/(RD-10*F) - 1/(RD+10*F));
% r =  RD*cos(theta).^2 * F * (1/(RD*cos(theta).^2-10*F) - 1/(RD*cos(theta).^2+10*F));
% 
% F = F+r/2;
% 
% end
% t =(RD/10)*(cos(theta).^2);
% if t>2*D
%     v = [t RD];
% else
%     v =RD;
% end
% rr =[rr   v];
% label1 = [ones(1,length(rr))*(sin(theta)) ; rr ];
% label = [label label1];
% end
% dict = zeros(N,size(label,2));
% for i= 1:size(label,2)
% Nt = N;
% theta = label(1,i);
% r = label(2,i);
% % at = near_field_manifold( Nt, d, fc, r, theta );
% at = polar_domain_manifold( Nt, d, fc, r, theta );
% dict(:,i) =at;
% end
% end


















% function [dict, label] = PCodeBook(N, D, RD, lambda)
% d = lambda/2;
% fc= 3e8/lambda;
% theta0 = linspace(-pi/2,pi/2, 2*N);
% 
% label = [];
% for  i =1:length(theta0)
% r_min = 3;    
% F = (RD/10)*cos(theta).^2;
% theta = theta0(i);
% rr= [];
% 
% while F >= r_min 
% rr = [F rr];
% r =  RD*cos(theta).^2 * F * (1/(RD*cos(theta).^2-10*F) - 1/(RD*cos(theta).^2+10*F));
% F = F+r/1.5;
% 
% end
% rr =[rr  RD];
% label1 = [ones(1,length(rr))*theta ; rr ];
% label = [label label1];
% end
% dict = zeros(N,size(label,2));
% for i= 1:size(label,2)
% Nt = N;
% theta = label(1,i);
% r = label(2,i);
% at = polar_domain_manifold( Nt, d, fc, r, theta );
% dict(:,i) =at;
% end
% end