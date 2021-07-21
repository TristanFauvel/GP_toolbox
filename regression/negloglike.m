function [negL, dnegL] = negloglike(theta, x_tr, y_tr, kernelfun, meanfun, varargin)
 

N = numel(y_tr);
% covariance and derivative
if nargout>1
    [K, dK] = kernelfun(theta.cov, x_tr,x_tr, true, 'nugget');
else
    K = kernelfun(theta.cov, x_tr,x_tr, true, 'nugget');
end
% 
% t=0;
% ill_conditioned =1;
% while ill_conditioned == 1
%     try chol(K);
%         ill_conditioned =0;
%     catch ME
%         disp('Matrix is not symmetric positive definite, added jitter of 1e-08 to the diagonal')
%         K= K + 1E-8*eye(size(K));
%     end
% end
%  
if nargout>1
     [prior_mean, dprior_mean_dtheta]=meanfun(x_tr, theta.mean);
else
    prior_mean=meanfun(x_tr, theta.mean);
end

L = chol(K)'; %in Rasmussen and Williams L = cholesky(K) is a lower triangular matrix, whereas in matlab it is an upper one
alpha= L'\(L\((y_tr(:)-prior_mean(:)))); %see Algorithm 2.25 in Williams and Rasmussen
negL= 0.5*(y_tr(:)-prior_mean(:))'*alpha + sum(log(diag(L))) + 0.5*N*log(2*pi); %eq 2.30, algo 2.1

nhyp=numel(theta.cov)+numel(theta.mean); %number of hyperparameters

invK = inv(L')*inv(L);
if nargout>1
    dnegL = zeros(nhyp, 1);
    for j = 1:numel(theta.cov)
          dnegL(j) = -0.5*trace((alpha*alpha'-invK)*dK(:, :, j));
    end
    dnegL(numel(theta.cov)+1:nhyp) = -dprior_mean_dtheta*alpha;
end

