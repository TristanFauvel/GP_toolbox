function B = nugget_regularization(A)
% Add a diagonal term to the covariance matrix A to reach a condition
% number lower than kmax, following Mohammadi et al 2017
kmax= 10^8;
if issymmetric(A) && cond(A)>kmax 
    lambda_max = max(real(eig(A)));
    lambda_min = min(real(eig(A)));
    if  (lambda_max-kmax*lambda_min)>0
        tau2 = (lambda_max-kmax*lambda_min)/(kmax - 1);
    else
        tau2 = 0;
    end
    B = A + tau2*eye(size(A,1));
else
    B= A;
end
if ~isreal(B)
    disp('bug')
end
%
% try chol(K);
% catch ME
%     disp('Matrix is not symmetric positive definite, added jitter of 1e-08 to the diagonal')
%     K= K + 1E-8*eye(size(K));
% end
