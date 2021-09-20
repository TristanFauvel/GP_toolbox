function W = Wasserstein2(m, K, input3, input4)
% 2-Wasserstein distances for GP

if nargin == 3
    Y = input3;
    mest = mean(Y);
    Kest = cov(Y);
elseif nargin ==4
    sqrtK = sqrtm(K);
    mest = input3;
    Kest = input4;
end
sqrtK = sqrtm(K);
W2 = sum((m - mest(:)).^2) + trace(K + Kest -2*sqrtm(sqrtK*Kest*sqrtK));

W= sqrt(W2);

