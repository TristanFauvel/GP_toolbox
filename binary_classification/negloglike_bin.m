function [negL, dnegL] = negloglike_bin(theta, xtrain, ctrain,kernelfun, varargin)
% c = 0 or 1
opts = namevaluepairtostruct(struct( ...
    'tol', 1e-7, ...
    'MaxIt', 1e3,...
    'modeltype', 'exp_prop', ...
    'regularization', 'nugget' ...
    ), varargin);

UNPACK_STRUCT(opts, false)

ctrain = ctrain(:);
% covariance and derivative
if nargout>1
    [K, dK] = kernelfun(theta, xtrain, xtrain, true, regularization);
else
    K = kernelfun(theta, xtrain, xtrain, true, regularization);
end
% toc

if isnan(sum(K(:)))
    disp('stop')
end

n = size(xtrain, 2);

if strcmp(modeltype, 'laplace')
    
    % find y= maximum of p(y_tr|c_tr, xtrain), (19.5.19 in Barber book)
    ystar = siteparams(ctrain, K, 'tol', tol, 'MaxIt', MaxIt, ...
        modeltype, post, regularization);
    
    % compute 'noise' matrix (Eq 19.5.17 in Barber book)
    D = diag(logistic(ystar).*(1-logistic(ystar)));
    
    % precompute C\y as you'll use it a lot in the future
    yoverC = pinv(K)*ystar;
    
    %%% negative log-likelihood (19.5.33 in Barber book)
    negL = -ctrain'*ystar ...
        + sum(log1pexp(ystar)) ...
        + 0.5*ystar'*yoverC...
        + 0.5*logdet(eye(n) + K*D);
    
elseif strcmp(modeltype, 'exp_prop')
    [nu_tilde, tau_tilde, nu_cav, tau_cav] ...
        = siteparams(ctrain, K, 'MaxIt', MaxIt, 'tol', tol);

    Stilde_half = diag(sqrt(tau_tilde));
    B = eye(n)+Stilde_half*K*Stilde_half;

    Lchol = chol(B);
    
    mu_cav = nu_cav./tau_cav;
    sigma2_cav = 1./tau_cav;
    
    T = diag(tau_cav);
    Stilde = diag(tau_tilde);
    
    z = (2*ctrain-1).*mu_cav./sqrt(1+sigma2_cav);
    
    negL = -(0.5*sum(log(1+tau_tilde./tau_cav)) - sum(log(diag(Lchol))) ...
        +  sum(logcdfnormal(z)) ...
        + 0.5*nu_tilde'*(K-K*(Stilde_half/B)*Stilde_half*K-diag(1./(tau_tilde+tau_cav)))*nu_tilde ...
        + 0.5*mu_cav'*T*diag(1./(tau_tilde+tau_cav))*(Stilde*mu_cav-2*nu_tilde));    

else
    error('modeltype must be ''laplace'' or ''exp_prop''')
end

% disp("derivative of negative log-likelihood")

if nargout>1    
    if strcmp(modeltype, 'laplace')        
        % dy/dtheta (Eq. 6.94 or Bishop book, Eq 19.5.38 of Barber book)
        dy = @(dK) (eye(n)+K*D)\(dK*(ctrain-logistic(ystar)));
        
        % d logdet(I+K*D)/dy (Eq 6.92 of Bishop book, Eq 19.55.36 of Barber book)
        term2 = diag((eye(n)+K*D)\K).*logistic(ystar).*(1-logistic(ystar)).*(1-2*logistic(ystar));
        
        % put together to get derivative of dL/dtheta (Eq 6.91 of Bishop book, Eq 19.5.35 of Bishop book)
        dnegL_fun = @(dK) - 0.5*yoverC'*dK*yoverC ...
            + 0.5*trace((eye(n)+K*D)\(dK*D))  ...
            + 0.5*term2'*dy(dK);

        % compute dL/dtheta for each hyperparameter theta.
        dnegL = zeros(size(theta));
        for i = 1:numel(theta)
            dnegL(i) = dnegL_fun(dK(:, :, i));
        end
        
    elseif strcmp(modeltype, 'exp_prop')        
        if nargout>1
            b=(eye(n)-(Stilde_half/B)*Stilde_half*K)*nu_tilde;
            
            dnegL= -0.5*permute( sum(sum((b*b' - (Stilde_half/B)*Stilde_half).*dK, 1), 2), [3 1 2]);
        end        
    else
        error('modeltype must be ''laplace'' or ''exp_prop''')
    end
end
% toc

return
