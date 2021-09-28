function [L, dLdx] = log_mvngauss(fx, mu,Sigma, dmudx, dSigmadx)

L = chol(Sigma)'; %in Rasmussen and Williams L = cholesky(K) is a lower triangular matrix, whereas in matlab it is an upper one
alpha= L'\(L\((fx-mu))); %see Algorithm 2.25 in Williams and Rasmussen
invK = inv(L')*inv(L);

m = numel(mu);
[D, nfx] = size(fx);
nx = size(dmudx,2);

negL= 0.5*dot((fx-mu),alpha) + sum(log(diag(L))) + 0.5*m*log(2*pi); %eq 2.30, algo 2.1

palpha = zeros(m,m,nfx);

dLdx = zeros(nfx, nx);
 for i = 1:nfx
    palpha(:,:,i) = alpha(:,i)*alpha(:,i)'-invK;
    for j = 1:nx        
        dLdx(i,j) = 0.5*trace(palpha(:,:,i)*dSigmadx(:,:,j)) + dmudx(:,j)'*alpha(:,i);
    end
end
 

L = -negL;
 