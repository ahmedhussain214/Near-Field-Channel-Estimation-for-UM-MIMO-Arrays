function [X, support,R] = SOMP(Y, Phi, s, N, M)
%%% initialization
R = Y;
support = [];
norm_Phi = sqrt(sum(abs(Phi).^2, 1)).';

%%% SOMP
for i = 1 : s
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
end