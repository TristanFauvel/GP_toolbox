classdef gp_classification_model < gpmodel
    properties
        link
        modeltype
        type
    end
    methods
        function model = gp_classification_model(D, meanfun, kernelfun, regularization, hyps, lb, ub, type, link, modeltype)
            model = model@gpmodel(D, meanfun, kernelfun, regularization, hyps, lb, ub);
            model.link = link;
            model.modeltype = modeltype;
            model.type = type;
        end
        function [output1,  mu_y, sigma2_y, Sigma2_y, dmuc_dx, dmuy_dx, dsigma2y_dx, dSigma2y_dx, var_muc, dvar_muc_dx,post] =  prediction(model, hyps, xtrain, ctrain, xtest, post)
            theta = hyps.cov;
            % mu_c corresponds to P(x1 > x2)
            kernelfun = model.kernelfun;
            modeltype = model.modeltype;
            regularization = model.regularization;
            link = model.link;
            
            if ~isempty(xtest) && size(xtrain,1) ~= size(xtest,1)
                error("Dimensions of test and training sets are not consistent")
            end
            dmuc_dx = [];
            dmuy_dx=[];
            dsigma2y_dx=[];
            dSigma2y_dx= [];
            MaxIt = 1e3;
            tol = 1e-7;
            
            if isempty(post) || nargout == 11 || isempty(xtest)
                clear('post')
                comp_post = true;
            else
                comp_post = false;
            end
            
            ctrain = ctrain(:);
            if any(isnan([xtest(:);xtrain(:)]))
                error('x is NaN')
            end
            [D,n] = size(xtrain);
            
            if isempty(ctrain)
                mu_y = zeros(size(xtest));
                Sigma2_y = kernelfun(theta, xtest, xtest, false, 'false');
                sigma2_y= diag(Sigma2_y);
                mu_c = 0.5*ones(size(xtest));
                output1 = mu_c;
                return
            end
            
            if any(ctrain ~= 0 & ctrain ~= 1)
                error('c should be either 0 or 1')
            end
            % kernel functions
            if comp_post
                K = kernelfun(theta, xtrain, xtrain, true, regularization);
                post.K = K;
            else
                K= post.K;
            end
            if ~isempty(xtest)
                ntst = size(xtest, 2);
                
                if nargout>4  && nargout ~= 9
                    [k, ~, dk_dx] = kernelfun(theta, xtrain, xtest, false, 'false');
                    post.dk_dx = dk_dx;
                else
                    k = kernelfun(theta, xtrain, xtest, false, 'false');
                end
                post.k = k;
            end
            
            if ~isempty(xtest)
                if nargout>4  && nargout ~= 9
                    [ks, ~, dks_dx] = kernelfun(theta, xtest, xtest, false, 'false');
                    diag_ks= diag(ks);
                    n=size(xtest,2);
                    L= reshape(dks_dx, n*n, [], D);
                    ddiagks_dx = L(1:n+1:end,:,:);
                else
                    
                    try ks = kernelfun(theta, xtest, xtest, false, 'false');
                        diag_ks= diag(ks);
                    catch
                        diag_ks = zeros(ntst,1);
                        for l = 1:ntst
                            diag_ks(l) = kernelfun(theta, xtest(:,l), xtest(:,l), false, 'false');
                        end
                    end
                end
            end
            
            if strcmp(modeltype, 'laplace')
                if comp_post
                    % find y= maximum of p(y_tr|c_tr, xtrain), (19.5.19 in Barber book)
                    ystar = siteparams(ctrain, K, link, 'tol', tol, 'MaxIt', MaxIt);
                    % compute 'noise' matrix (Eq 19.5.17 in Barber book)
                    
                    if strcmp(func2str(link),'logistic')
                        D = diag(logistic(ystar).*(1-logistic(ystar)));
                        invS = diag(1./(logistic(ystar).*(1-logistic(ystar))+1e-12));
                        dloglike = (ctrain - logistic(ystar));
                    elseif strcmp(func2str(link),'normcdf')
                        y = 2*ctrain-1;
                        a = normpdf(ystar);
                        b = normcdf(y.*ystar);
                        D = diag((a./b).^2 + y.*ystar.*(a./b));
                        dloglike = y.*a./b;
                        
                        invS = diag(1./((a./b).^2 + y.*ystar.*(a./b)+1e-12));
                    end
                    
                    invKS = inv(K + invS);
                    
                    post.ystar = ystar;
                    post.invKS = invKS;
                    post.invS = invS;
                    post.dloglike =  dloglike;
                    post.D = D;
                else
                    ystar = post.ystar;
                    invKS = post.invKS;
                    invS = post.invS;
                    dloglike = post.dloglike;
                end
                
            elseif strcmp(modeltype, 'exp_prop')
                if comp_post
                    [nu_tilde, tau_tilde] = siteparams(ctrain, K, link, 'modeltype', 'exp_prop', 'tol', tol);
                    
                    invS = diag(1./tau_tilde);
                    invKS = inv(K+invS);
                    
                    dloglike = invKS*(invS*nu_tilde); % mu_y = (k'*invKS)*invS*nu_tilde;
                    
                    post.nu_tilde = nu_tilde;
                    post.invS = invS;
                    post.invKS = invKS;
                    post.tau_tilde = tau_tilde;
                    post.dloglike = dloglike;
                else
                    invKS = post.invKS;
                    nu_tilde = post.nu_tilde;
                    invS = post.invS;
                    dloglike = post.dloglike;
                end
            else
                error('modeltype must be ''laplace'' or ''exp_prop''')
            end
            
            
            if ~isempty(xtest)
                
                % mean of p(y|x, c_tr, xtrain) (Eq 19.5.24 in Barber book)
                mu_y = k'*dloglike;
                
                % standard deviation of p(y|x, c_tr, xtrain) (Eq 19.5.26 in Barber book)
                %sigma2_y = diag_ks-sum((k'/(K + invS)).*k', 2);
                sigma2_y = diag_ks-sum((k'*invKS).*k', 2);
                sigma2_y(sigma2_y<0) = 0;
                if nargout>=4
                    %Sigma2_y = ks-k'*inv(K + invS)*k;
                    Sigma2_y = ks-k'*invKS*k;
                    Sigma2_y = (Sigma2_y + Sigma2_y')/2;
                    if strcmp(model.regularization, 'nugget')
                        Sigma2_y = nugget_regularization(Sigma2_y);
                    end
                end
                % mean of p(c|x, c_tr, xtrain) (Eq 19.5.27 in Barber book)
                if strcmp(func2str(link),'logistic')
                    nu = sqrt(pi)/4;
                    mu_c = 0.5+0.5*erf(nu.*mu_y./sqrt(1+2*nu*sigma2_y));
                elseif strcmp(func2str(link),'normcdf')
                    mu_c =  normcdf(mu_y./sqrt(1+sigma2_y));
                end
                
                if nargout>4 && nargout ~= 9
                    if strcmp(func2str(link),'logistic')
                        sigma_y = sqrt(sigma2_y);
                        % derivative of mu_c with respect to mu_y and sigma_y
                        arg = nu.*mu_y.*(1+2*nu*sigma_y.^2).^(-0.5);
                        
                        darg_muy = nu./sqrt(1+2*nu*sigma_y.^2);
                        
                        darg_sigma2y = -nu^2.*mu_y.*(1 + 2*nu*sigma_y.^2).^(-1.5);
                        
                        dmuc_dmuy = exp(-arg.^2)./sqrt(pi).*darg_muy;
                        
                        dmuc_dsigma2y = exp(-arg.^2)./sqrt(pi).*darg_sigma2y;
                        
                        
                    elseif strcmp(func2str(link),'normcdf')
                        
                        % derivative of mu_c with respect to mu_y and sigma_y
                        arg = normpdf(mu_y./sqrt(1+sigma2_y));
                        
                        dmuc_dmuy = arg./sqrt(1+sigma2_y);
                        
                        dmuc_dsigma2y = -0.5*arg.*mu_y./((1+sigma2_y).^(1.5));
                    else
                        error('modeltype must be ''laplace'' or ''exp_prop''')
                    end
                    % derivative of mu_y with respect to x
                    dmuy_dx = squeeze(mtimesx(dk_dx, 'T', post.dloglike));
                    
                    dsigma2y_dx = 2*squeeze(ddiagks_dx) -2*squeeze(sum(mtimesx(dk_dx, 'T', invKS).*k', 2));
                    dSigma2y_dx = squeeze(dks_dx) - 2*mtimesx(mtimesx(dk_dx, 'T', invKS),k);
                    
                    % derivative of mu_c with respect to x
                    dmuc_dx = dmuc_dsigma2y.*dsigma2y_dx + dmuc_dmuy.*dmuy_dx;
                    
                end
                
                
                if nargout >= 9
                    if strcmp(func2str(link),'normcdf')
                        h = mu_y./sqrt(1+sigma2_y);
                        a = 1./sqrt(1+2*sigma2_y);
                        
                        [tfn_output, dTdh, dTda] = tfn(h, a);
                        
                        if nargout > 9
                            indexat = @(expr, varargin) expr(varargin{:});
                            dmcdx = indexat(reshape(dmuc_dx, size(xtest,2)^2, D),1:size(xtest,2)+1:(size(xtest,2)^2), 1:D);
                            dmydx =  indexat(reshape(dmuy_dx, size(xtest,2)^2, D),1:size(xtest,2)+1:(size(xtest,2)^2), 1:D);
                            ds2dx = indexat(reshape(dsigma2y_dx, size(xtest,2)^2, D),1:size(xtest,2)+1:(size(xtest,2)^2), 1:D);
                            dhdx = dmydx./sqrt(1+sigma2_y) - 0.5*mu_y.*ds2dx.*((1+sigma2_y).^(-3/2));
                            dadx =  -ds2dx.*((1+2*sigma2_y).^(-1.5));
                            dTdx = dTdh.*dhdx +dTda.*dadx;
                            dvar_muc_dx = dmcdx-2*mu_c.*dmcdx-2*dTdx;
                        end
                        var_muc = (mu_c - 2*tfn_output) - mu_c.^2;
                        
                        var_muc(sigma2_y==0) = 0;
                    else
                        var_muc = [];
                        warning('The computation of $V(\pi(x))$ is not implemented for this link function')
                    end
                end
            end
            
            if isempty(xtest)
                output1 = post;
            else
                output1 = mu_c;
            end
        end
        
        function [negL, dnegL] = negloglike(model, hyps, xtrain, ctrain)
            theta = hyps.cov;
            tol = 1e-7;
            MaxIt =  1e3;
            
            
            kernelfun = model.kernelfun; 
            regularization = model.regularization; 
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
            
            if strcmp(model.modeltype, 'laplace')
                
                % find y= maximum of p(y_tr|c_tr, xtrain), (19.5.19 in Barber book)
                ystar = siteparams(ctrain, K, model.link, 'tol', tol, 'MaxIt', MaxIt);
                
                % compute 'noise' matrix (Eq 19.5.17 in Barber book)
                
                if strcmp(func2str(model.link),'logistic')
                    D = diag(logistic(ystar).*(1-logistic(ystar)));
                    dloglike = (ctrain - logistic(ystar));
                elseif strcmp(func2str(model.link),'normcdf')
                    y = 2*ctrain-1;
                    a = normpdf(ystar);
                    b = normcdf(y.*ystar);
                    D = diag((a./b).^2 + y.*ystar.*(a./b));
                    dloglike = y.*a./b;
                end
                
                % precompute C\y as you'll use it a lot in the future
                yoverC = pinv(K)*ystar;
                
                %%% negative log-likelihood (19.5.33 in Barber book)
                
                if strcmp(func2str(model.link),'logistic')
                    negL = -ctrain'*ystar ...
                        + sum(log1pexp(ystar)) ...
                        + 0.5*ystar'*yoverC...
                        + 0.5*logdet(eye(n) + K*D);
                elseif strcmp(func2str(model.link),'normcdf')
                    negL =  sum(ctrain.*log(normcdf(ystar)) + (1-ctrain).*log(1-normcdf(ystar))) ...
                        + 0.5*ystar'*yoverC...
                        + 0.5*logdet(eye(n) + K*D);
                end
                
                
                
                
            elseif strcmp(model.modeltype, 'exp_prop')
                [nu_tilde, tau_tilde, nu_cav, tau_cav] ...
                    = siteparams(ctrain, K, model.link, 'MaxIt', MaxIt, 'tol', tol);
                
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
                if strcmp(model.modeltype, 'laplace')
                    % dy/dtheta (Eq. 6.94 or Bishop book, Eq 19.5.38 of Barber book)
                    dy = @(dK) (eye(n)+K*D)\(dK*dloglike);
                    
                    % d logdet(I+K*D)/dy (Eq 6.92 of Bishop book, Eq 19.55.36 of Barber book)
                    if strcmp(func2str(model.link),'logistic')
                        term2 = diag((eye(n)+K*D)\K).*logistic(ystar).*(1-logistic(ystar)).*(1-2*logistic(ystar));
                        
                    elseif strcmp(func2str(model.link),'normcdf')
                        d3logphi = @(f) (4.*2.^(1./2).*exp(-(3.*f.^2.*(2.*ctrain - 1).^2)./2).*(2.*ctrain - 1).^3)./(pi.^(3./2).*(erf((2.^(1./2).*f.*(2.*ctrain - 1))./2) + 1).^3) - ...
                            (2.^(1./2).*exp(-(f.^2.*(2.*ctrain - 1).^2)./2).*(2.*ctrain - 1).^3)./(pi.^(1./2).*(erf((2.^(1./2).*f.*(2.*ctrain - 1))./2) + 1)) + ...
                            (6.*f.*exp(-f.^2.*(2.*ctrain - 1).^2).*(2.*ctrain - 1).^4)./(pi.*(erf((2.^(1./2).*f.*(2.*ctrain - 1))./2) + 1).^2) + ...
                            (2.^(1./2).*f.^2.*exp(-(f.^2.*(2.*ctrain - 1).^2)./2).*(2.*ctrain - 1).^5)./(pi.^(1./2).*(erf((2.^(1./2).*f.*(2.*ctrain - 1))./2) + 1));
                        
                        term2 = diag((eye(n)+K*D)\K).*d3logphi(ystar);
                    end
                    % put together to get derivative of dL/dtheta (Eq 6.91 of Bishop book, Eq 19.5.35 of Barber book)
                    dnegL_fun = @(dK) - 0.5*yoverC'*dK*yoverC ...
                        + 0.5*trace((eye(n)+K*D)\(dK*D))  ...
                        + 0.5*term2'*dy(dK);
                    
                    % compute dL/dtheta for each hyperparameter theta.
                    dnegL = zeros(model.ncov_hyp + model.nmean_hyp,1);
                    for i = 1:model.ncov_hyp
                        dnegL(i) = dnegL_fun(dK(:, :, i));
                    end
                    
                elseif strcmp(model.modeltype, 'exp_prop')
                    if nargout>1
                        dnegL = zeros(model.ncov_hyp + model.nmean_hyp,1);
                        b=(eye(n)-(Stilde_half/B)*Stilde_half*K)*nu_tilde;
                        
                        dnegL(1:model.ncov_hyp) = -0.5*permute( sum(sum((b*b' - (Stilde_half/B)*Stilde_half).*dK, 1), 2), [3 1 2]);
                    end
                else
                    error('modeltype must be ''laplace'' or ''exp_prop''')
                end
            end
        end
    end
end