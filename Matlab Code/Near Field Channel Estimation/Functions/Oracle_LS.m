function H_sf = Oracle_LS(Y, Phi, N, M, r, theta, d, f)
    L = size(r, 2);
    A = zeros(N, L);
    for i = 1:L
%        A(:, i) = polar_domain_manifold( N, d, f, r(i), theta(i) );
        A(:, i) = near_field_manifold( N, d, f, r(i), theta(i) );
    end
    % Y = Phi * A * X
    B = Phi * A;
    X = pinv(B'*B)*B' * Y;
%     H_sf = zeros(N, M);
    H_sf = A * X;
end

