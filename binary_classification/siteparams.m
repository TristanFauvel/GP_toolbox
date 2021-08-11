function [siteparam1, siteparam2, siteparam3, siteparam4, mu, Sigma] = siteparams(c, K, link, varargin)

opts = namevaluepairtostruct(struct( ...
    'tol', 1e-7, ...
    'MaxIt', 1e3...
    ), varargin);

UNPACK_STRUCT(opts, false)
n = numel(c);

ystar = zeros(n, 1);
tau_tilde =  zeros(n, 1);
nu_tilde =  zeros(n, 1);
    

c = c(:);

if strcmp(func2str(link),'logistic')
    dloglike = @(ystar) c-logistic(ystar);
elseif strcmp(func2str(link),'normcdf')
    y = 2*c-1;
    dloglike =@(ystar)  y.*normpdf(ystar)./ normcdf(y.*ystar);
end

if nargout == 1 %% Laplace approximation
    
    dy = 100; count = 0;
    while dy>tol && count<MaxIt
        count = count+1;
        
        if strcmp(func2str(link),'logistic')
            D = diag(logistic(ystar).*(1-logistic(ystar)));
        elseif strcmp(func2str(link),'normcdf')
            a = normpdf(ystar);
            b = normcdf(y.*ystar);
            D = diag((a./b).^2 + y.*ystar.*(a./b));
        end
        ystar =  (eye(n)+K*D)\(K*(dloglike(ystar)+D*ystar));
        dy = mean((ystar(:)-K*dloglike(ystar)).^2);
    end
    
    siteparam1 = ystar;
    
    
elseif nargout > 1 %% Expectation propagation
    if strcmp(func2str(link),'logistic')
       error('Expectation propagation not implement for a logistic link function') 
    end
    
    Sigma = K;
    mu = zeros(n, 1);
    
    tau_tilde0 = tau_tilde;  nu_tilde0 = nu_tilde;
    
    tau = zeros(n, 1);    nu = zeros(n, 1);
    
    dL = 100; count =0; dparams = 100;
    
    while dparams>tol && count<MaxIt
        count = count+1;
        
        for i = 1:n
            
            tau(i) = 1/Sigma(i,i)     - tau_tilde(i);
            nu(i)  = 1/Sigma(i, i)*mu(i)-nu_tilde(i);
            
            [mu_hat, sigma2_hat] = compute_moments(nu(i)./tau(i), 1./tau(i), c(i));
            
            dtau_tilde   = 1./sigma2_hat - tau(i) - tau_tilde(i);
            tau_tilde(i) = tau_tilde(i) + dtau_tilde;
            
            nu_tilde(i)   = 1./sigma2_hat.*mu_hat - nu(i);
            
            Sigma = Sigma - Sigma(:, i)*Sigma(:, i)'./(1./dtau_tilde + Sigma(i, i));
            
            mu = Sigma*nu_tilde;
            
        end
        
        Stilde_half = diag(sqrt(tau_tilde));
        B = eye(n)+Stilde_half*K*Stilde_half;
        
        L = chol(B);
        V = L'\Stilde_half*K;
        
        Sigma = K-V'*V;
        
        mu = Sigma*nu_tilde;
        
        dparams = mean( (tau_tilde-tau_tilde0).^2 ) ...
            + mean( (nu_tilde-nu_tilde0).^2 );
        
        nu_tilde0 = nu_tilde; tau_tilde0 = tau_tilde;
        
    end
    
    siteparam1 = nu_tilde;
    siteparam2 = tau_tilde;
    siteparam3 = nu;
    siteparam4 = tau;
    
end
end

function [mu_hat, sigma2_hat] = compute_moments(mu, sigma2, c)


z = (2*c-1).*mu./sqrt(1+sigma2);

Zhat    = normcdf(z);
normz   = normpdf(z);

mu_hat = mu + (2*c-1).*sigma2.*normz./(Zhat.*sqrt(1+sigma2));

sigma2_hat = sigma2 ...
    - (sigma2.^2.*normz./(Zhat.*(1+sigma2))).*(z + normz./Zhat);

end





