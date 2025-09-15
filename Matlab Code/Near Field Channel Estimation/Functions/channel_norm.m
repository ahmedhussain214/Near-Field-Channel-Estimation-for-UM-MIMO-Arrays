function hn = channel_norm(h)
[K, N, M] = size(h);

hn = zeros(size(h));

for k = 1:K
    hk = h(k, :, :);
    hn(k, :, :) = hk / sqrt(sum(abs(hk(:)).^2)) * sqrt(N * M);
end


end

