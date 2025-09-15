function [A, G, r, theta, Res_list] = P_SIGW_N(Y, W, A0, G0, r0, theta0,  f, d, L_min, pruning, N_iter, lr_alpha, lr_th, th)
%IR_Near Y = WH * A * G + n

% [A4, G4, r4, theta4] =P_SIGW( Z/Znorm, Phi'/Znorm, Polarcodebook(:, sup3), Haf_hat3(sup3,:), label(2, sup3), label(1, sup3), ...
%             fc, d, L, 3, 20, 10, 1, 1e-8


[Q, M] = size(Y);
[N, L] = size(A0);
% Ini
A = A0;
G = G0;
G_prev = G;
theta0 = sin(theta0);
theta = theta0;
r = r0;
alpha = (1 - theta.^2) ./ 2 ./ r;

c = 3e8;
row = (-(N - 1)/2:(N - 1)/2)' ;
% lambda = 10;

Res = Y - W'*A*G;
Res_list = [];
Rnorm = norm(Res, 'fro')^2;
Res_list = [Res_list, Rnorm];

c1 = 0.1;
cl_alpha = 0.01;
% calculate gradient


for iter = 1:N_iter
   
    %% update theta
    dtheta0 = df_dtheta(Y, W, A, r, theta, f,d);
    direction_theta = dtheta0;
    lr_th_iter = lr_th;
    
    while lr_th_iter > 1e-6  
        theta_iter = theta - lr_th_iter * direction_theta;
        theta_iter(theta_iter > 1) = theta_iter(theta_iter > 1) - 2;
        theta_iter(theta_iter < -1) = theta_iter(theta_iter < -1) + 2;
        
        
        
        % update A
        rn_iter = - ( row * d ) .* theta_iter + ( row * d ).^2 .* alpha;
       
        A_iter = exp( - 1j*2*pi*f*( rn_iter )/c ) / sqrt(N);
        R_iter = W' * A_iter;
        G_iter = pinv( R_iter'*R_iter + 1e-15*eye(size(R_iter, 2)) ) * R_iter' * Y;
        
        Res2 = Y - W'*A_iter*G_iter;
        Rnorm2 = norm(Res2, 'fro')^2;
        if Rnorm2 < Rnorm - c1 * lr_th_iter * norm(direction_theta, 'fro')^2
            break;
        else
            lr_th_iter = lr_th_iter * 0.5;
        end
     end

    % update theta(A) and gains(G)
    theta = theta - lr_th_iter * direction_theta;
    theta(theta > 1) = theta(theta > 1) - 2;
    theta(theta < -1) = theta(theta < -1) + 2;

    
    % update A
    rn = - ( row * d ) .* theta + ( row * d ).^2 .* alpha;
    A = exp( - 1j*2*pi*f*( rn )/c) / sqrt(N);

    R = W' * A;
    G = pinv( R'*R + 1e-15 * eye(size(R, 2)) ) * R' * Y;
    
    %% update distance
%     dr0 = df_dr(Y, W, A, r, theta,f,d);
    dalpha0 = df_dalpha(Y, W, A, r, theta,f,d);
    direction_alpha = dalpha0;
    lr_alpha_iter = lr_alpha;
    while lr_alpha_iter > 1e-6  
        alpha_iter = alpha - lr_alpha_iter * direction_alpha;
        alpha_iter( alpha_iter < 0 ) = 1e-10;

        % update A
        rn_iter =  - ( row * d ) .* theta + ( row * d ).^2 .* alpha_iter;
        A_iter = exp( - 1j*2*pi*f*( rn_iter  )/c)/sqrt(N);

        R_iter = W' * A_iter;
        G_iter = pinv( R_iter'*R_iter + 1e-15*eye(size(R_iter, 2))) * R_iter' * Y;
        
        Res2 = Y - W'*A_iter*G_iter;
        Rnorm2 = norm(Res2, 'fro')^2;
        if Rnorm2 < Rnorm - cl_alpha * lr_alpha_iter * norm(direction_alpha, 'fro')^2
            break;
        else
            lr_alpha_iter = lr_alpha_iter * 0.5;
        end
    end
     
    % update distanace(r) and gains(G)
    alpha = alpha - lr_alpha_iter * direction_alpha;
    alpha( alpha <= 0 ) = 1e-10;

    % update A
    rn = r - ( row * d ) .* theta + ( row * d ).^2 .* alpha;
    A = exp( - 1j*2*pi*f*( rn )/c) / sqrt(N);

    R = W' * A;
    G = pinv( R'*R + 1e-15 * eye(size(R, 2)) ) * R' * Y;

    % obtain res and res norm
    Res = Y - W'*A*G;
    Rnorm = norm(Res, 'fro')^2;
    Res_list = [Res_list, Rnorm];
    
    % update epsilon
    gamma = norm(G - G_prev, 'fro');
  
    % pruning
    %if L > L_min && iter > 10
    if L > L_min
        Gnorm = sum(abs(G).^2, 2);
        index = find(Gnorm > pruning);
%         index = find( Gnorm > pruning * max(Gnorm));
        r = r(index);
        theta = theta(index);
        alpha = alpha(index);
        A = A(:, index);
        G = G(index, :);
        dtheta0 = dtheta0(index);
        dalpha0 = dalpha0(index);
        direction_theta = direction_theta(index);
        direction_alpha = direction_alpha(index);
        L = numel(index);
    end


    % early stopping
    if gamma < th
       break; 
    end
    % update previous G
    G_prev = G;  
end

end


function df = df_dtheta(Y, W, A, r, theta, f,d)
% f = - tr{ Y'R(R'R + D/lambda)^(-1)R'Y}
% R = W'A
    dA = dA_dtheta(A, r, theta, f,d);
    R = W'*A;
%     size(R'*R)
%     size(D)
    E = R'*R ;
    EI = pinv(E);
    [N, L] = size(A);
    M = size(Y, 2);
    df = zeros(1, L);
    for l = 1:L
        dA_dl = zeros(N, L);
        dA_dl(:, l) = dA(:, l);
        dR_dl = W'*dA_dl;
        dE_dl = dR_dl' * R + R' * dR_dl;
        dEI_dl = - EI * dE_dl * EI;
        dP_dl = dR_dl * EI * R' + R * dEI_dl * R' + R * EI * dR_dl';
        df_dl = - trace( Y' * dP_dl * Y);
        df(l) = real(df_dl)/M;
    end
end


function df = df_dr(Y, W, A, r, theta, f,d)
% f = - tr{ Y'R(R'R + D/lambda)^(-1)R'Y}
% R = W'A
    dA = dA_dr(A, r, theta, f,d);
    R = W'*A;
    E = R'*R;
    EI = pinv(E);
    [N, L] = size(A);
    M = size(Y, 2);
    df = zeros(1, L);
    for l = 1:L
        dA_dl = zeros(N, L);
        dA_dl(:, l) = dA(:, l);
        dR_dl = W'*dA_dl;
        dE_dl = dR_dl' * R + R' * dR_dl;
        dEI_dl = - EI * dE_dl * EI;
        dP_dl = dR_dl * EI * R' + R * dEI_dl * R' + R * EI * dR_dl';
        df_dl = - trace( Y' * dP_dl * Y);
        df(l) = real(df_dl)/M;
    end
end

function df = df_dalpha(Y, W, A, r, theta, f,d)
% f = - tr{ Y'R(R'R + D/lambda)^(-1)R'Y}
% R = W'A
    dA = dA_dalpha(A, r, theta, f,d);
    R = W'*A;
    E = R'*R;
    EI = pinv(E);
    [N, L] = size(A);
    M = size(Y, 2);
    df = zeros(1, L);
    for l = 1:L
        dA_dl = zeros(N, L);
        dA_dl(:, l) = dA(:, l);
        dR_dl = W'*dA_dl;
        dE_dl = dR_dl' * R + R' * dR_dl;
        dEI_dl = - EI * dE_dl * EI;
        dP_dl = dR_dl * EI * R' + R * dEI_dl * R' + R * EI * dR_dl';
        df_dl = - trace( Y' * dP_dl * Y);
        df(l) = real(df_dl)/M;
    end
end

function dA = dA_dtheta(A, r, theta, f,d)
    N = size(A, 1);
    c = 3e8;
    row = (-(N - 1)/2:(N - 1)/2)' ;
    rn = sqrt( r.^2 + ( row * d ).^2 - 2 * row * d * ( r.*theta ));
%     dA = (-1j*2*pi*f/c)*(  -  (row * d * r)./rn  ).* A;
%     dA =(1j*2*pi*f/c).*( row * d  + (row * d).^2 * (theta./ r))  .* A;
    dA =(1j*2*pi*f/c).*( row * d )  .* A;
end


function dA = dA_dr(A, r, theta, f, d)
    N = size(A, 1);
    c = 3e8;
    row = (-(N - 1)/2:(N - 1)/2)' ;
    rn = sqrt( r.^2 + ( row * d ).^2 - 2 * row * d * ( r.*theta ));
%     dA = (-1j*2*pi*f/c).*( r./ rn - row * d * theta ./rn - 1) .* A;
    dA = (1j*2*pi*f/c).*( (row * d).^2 * (1 - theta.^2)./ (2 * r.^2) ) .* A;
end


function dA = dA_dalpha(A, r, theta, f, d)
    N = size(A, 1);
    c = 3e8;
    row = (-(N - 1)/2:(N - 1)/2)' ;
    rn = sqrt( r.^2 + ( row * d ).^2 - 2 * row * d * ( r.*theta ));
%     dA = (-1j*2*pi*f/c).*( r./ rn - row * d * theta ./rn - 1) .* A;
    dA = (1j*2*pi*f/c).*( - (row * d).^2  ) .* A;
end
