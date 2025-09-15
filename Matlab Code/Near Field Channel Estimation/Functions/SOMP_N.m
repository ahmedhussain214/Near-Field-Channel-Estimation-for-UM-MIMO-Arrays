function [X, support,i] = SOMP_N(Y, Phi, s, N, M,sigma)
%%% initialization
R = Y;
support = [];
norm_Phi = sqrt(sum(abs(Phi).^2, 1)).';
% n= length(Y)*(3.1);%% Q_32
% n = length(Y)*(2.5); %% Q_64
n = length(Y)* (1.5+288/length(Y));
e = sqrt(sigma)*(sqrt(n + sqrt(n*log(n))));
i=1;
r = e+1;
%%% SOMP
while (r>e)&&(i<3*s)
%     for n = 1 : N
%         t(n) = norm(Phi(:,n)'* R,'fro')^2;
%     end
    t = Phi'*R./norm_Phi;
    t = sum(abs(t).^2, 2);
    [~,order] = max(t);
    support = [support order];
    Phi_s = Phi(:,support);
    X = zeros(N,M);
    X(support,:) = pinv(Phi_s'*Phi_s)*Phi_s'*Y;
    R = Y - Phi*X;
    r = norm(R,2);
    i =i+1;
end