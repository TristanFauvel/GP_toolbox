function [L, dLdx] = log_mvngauss(fx, mu,Sigma, dmudx, dSigmadx)
m = numel(mu);
[D, nfx] = size(fx);

L = chol(Sigma)'; %in Rasmussen and Williams L = cholesky(K) is a lower triangular matrix, whereas in matlab it is an upper one
alpha= L'\(L\((fx-mu))); %see Algorithm 2.25 in Williams and Rasmussen

invK = inv(L')*inv(L);

nx = numel(mu); % dimension of x

negL= 0.5*dot((fx-mu),alpha) + sum(log(diag(L))) + 0.5*m*log(2*pi); %eq 2.30, algo 2.1
L = -negL;

dLdx = zeros(nfx, nx);
for i = 1:nfx
    % dLdx(i,:) = 0.5*trace(mtimesx((alpha(:,i)*alpha(:,i)'-invK),dSigmadx)) + dmudx'*alpha(:,i);
    A = (alpha(:,i)*alpha(:,i)'- invK);
    for j=1:nx
        dLdx(i,j) = 0.5*trace(A*dSigmadx(:,:,j)) + dmudx(:,j)'*alpha(:,i);
    end
end



