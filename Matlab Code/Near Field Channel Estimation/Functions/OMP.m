function [X, support] = OMP( Y, Phi, s, N, M )
%%% initialization
R = Y;
%%% OMP
X = zeros(N, M);
A = zeros(s, M);
norm_Phi = sqrt(sum(abs(Phi).^2, 1)).';
for m = 1 : M
    support = [];
    for i = 1 : s 
        [~,order] = max(abs(Phi' * R(:,m))./norm_Phi);
        support = [support order];
        Phi_s = Phi(:,support);
        Temp = zeros(N,1);
        Temp(support,:) = pinv( Phi_s' * Phi_s )*Phi_s'*Y(:,m);
        R(:,m) = Y(:,m) - Phi*Temp;
    end
    A(:,m) = unique(support);
    X(:,m) = Temp;
end