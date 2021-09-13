function W = Wasserstein2(m, K, Y)
% 2-Wasserstein distances for GP

mest = mean(Y);
Kest = cov(Y);
sqrtK = sqrtm(K);
W2 = sum((m - mest(:)).^2) + trace(K + Kest -2*sqrtm(sqrtK*Kest*sqrtK));
W= sqrt(W2);

