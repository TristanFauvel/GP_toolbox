function [negL, dnegL] = negloglike_bin(theta, xtrain, ctrain, model, varargin)
% c = 0 or 1
opts = namevaluepairtostruct(struct( ...
    'tol', 1e-7, ...
    'MaxIt', 1e3 ...
    ), varargin);

UNPACK_STRUCT(opts, false)

kernelfun = model.kernelfun;
modeltype = model.modeltype;
regularization = model.regularization;
link = model.link;
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
    ystar = siteparams(ctrain, K, link, 'tol', tol, 'MaxIt', MaxIt);
    
    % compute 'noise' matrix (Eq 19.5.17 in Barber book)
    
    if strcmp(func2str(link),'logistic')
        D = diag(logistic(ystar).*(1-logistic(ystar)));
        dloglike = (ctrain - logistic(ystar));
    elseif strcmp(func2str(link),'normcdf')
        y = 2*c-1;
        a = normpdf(ystar);
        b = normcdf(y.*ystar);
        D = diag((a./b).^2 + y.*ystar.*(a./b));
        dloglike = y.*a./b;
    end
    
    % precompute C\y as you'll use it a lot in the future
    yoverC = pinv(K)*ystar;
    
    %%% negative log-likelihood (19.5.33 in Barber book)
    
    if strcmp(func2str(link),'logistic')
        negL = -ctrain'*ystar ...
            + sum(log1pexp(ystar)) ...
            + 0.5*ystar'*yoverC...
            + 0.5*logdet(eye(n) + K*D);
    elseif strcmp(func2str(link),'normcdf')
        negL =  sum(ctrain.*log(normcdf(ystar)) + (1-ctrain).*log(1-normcdf(ystar))) ...
            + 0.5*ystar'*yoverC...
            + 0.5*logdet(eye(n) + K*D);
    end
    
    
    
    
elseif strcmp(modeltype, 'exp_prop')
    [nu_tilde, tau_tilde, nu_cav, tau_cav] ...
        = siteparams(ctrain, K, link, 'MaxIt', MaxIt, 'tol', tol);
    
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
        dy = @(dK) (eye(n)+K*D)\(dK*dloglike);
        
        % d logdet(I+K*D)/dy (Eq 6.92 of Bishop book, Eq 19.55.36 of Barber book)
        if strcmp(func2str(link),'logistic')
            term2 = diag((eye(n)+K*D)\K).*logistic(ystar).*(1-logistic(ystar)).*(1-2*logistic(ystar));
            
        elseif strcmp(func2str(link),'normcdf')
            d3logphi = @(f) (4.*2.^(1./2).*exp(-(3.*f.^2.*(2.*c - 1).^2)./2).*(2.*c - 1).^3)./(pi.^(3./2).*(erf((2.^(1./2).*f.*(2.*c - 1))./2) + 1).^3) - ...
                (2.^(1./2).*exp(-(f.^2.*(2.*c - 1).^2)./2).*(2.*c - 1).^3)./(pi.^(1./2).*(erf((2.^(1./2).*f.*(2.*c - 1))./2) + 1)) + ...
                (6.*f.*exp(-f.^2.*(2.*c - 1).^2).*(2.*c - 1).^4)./(pi.*(erf((2.^(1./2).*f.*(2.*c - 1))./2) + 1).^2) + ...
                (2.^(1./2).*f.^2.*exp(-(f.^2.*(2.*c - 1).^2)./2).*(2.*c - 1).^5)./(pi.^(1./2).*(erf((2.^(1./2).*f.*(2.*c - 1))./2) + 1));

            term2 = diag((eye(n)+K*D)\K).*d3logphi(ystar);
        end
        % put together to get derivative of dL/dtheta (Eq 6.91 of Bishop book, Eq 19.5.35 of Barber book)
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
