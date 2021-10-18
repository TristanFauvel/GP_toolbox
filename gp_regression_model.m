classdef gp_regression_model < gpmodel   
    properties
       type = 'regression'; 
    end
    methods
        %%
        function [output1, sigma2_y,dmu_dx, dsigma2_dx, Sigma2_y, dSigma2_dx, post] =  prediction(model, theta, xtrain, ytrain, xtest, post)      
            regularization = model.regularization;
            kernelfun = model.kernelfun;
            meanfun = model.meanfun;
            
            if  isempty(post) || nargout == 7
                clear('post')
                comp_post=  true;
            else
                comp_post= false;
            end
            
            if isempty(xtrain)
                mu_y = meanfun(xtest, theta.mean)';
                Sigma2_y = kernelfun(theta.cov, xtest, xtest, true, regularization);
                sigma2_y = diag(Sigma2_y);
                dmu_dx = [];
                dsigma2_dx = [];
                dSigma2_dx = [];
                post = [];
                output1 = mu_y;
            else
                
                ytrain = ytrain(:);
                ntest = size(xtest,2);
                nd=size(xtrain,1);
                if comp_post
                    K = kernelfun(theta.cov, xtrain, xtrain, true, regularization);
                    L = chol(K)'; %in Rasmussen and Williams L = cholesky(K) is a lower triangular matrix, whereas in matlab it is an upper one
                    invK = inv(K);
                    post.invK = invK;
                    post.L = L;
                    post.K = K;
                else
                    K = post.K;
                    L = post.L;
                    invK = post.invK;
                end
                
                if ~isempty(xtest)
                    if nargout>2
                        [k, ~, dk_dx] = kernelfun(theta.cov, xtrain, xtest, false, 'no');
                    else
                        k = kernelfun(theta.cov, xtrain, xtest, false, 'no');
                        
                    end
                    
                    prior_mean_tr=meanfun(xtrain, theta.mean)';
                    
                    if nargout > 2
                        [prior_mean_test, ~, dprior_mean_test_dx]=meanfun(xtest, theta.mean);
                        [ks, ~,dks_dx] = kernelfun(theta.cov, xtest, xtest, false, 'no');
                        V= reshape(dks_dx, ntest*ntest, [],nd);
                        ddiagks_dx = V(1:ntest+1:end,:,:);
                        
                    else
                        prior_mean_test=meanfun(xtest, theta.mean);
                        ks = kernelfun(theta.cov, xtest, xtest, false, 'no');
                        
                    end
                    
                    alpha= L'\(L\(ytrain-prior_mean_tr(:))); %see Algorithm 2.1 (p19), eq 2.25 in Williams and Rasmussen, alpha = inv(K)*(ytrain-mean)
                    mu_y = prior_mean_test(:) + k'*alpha;
                    
                    diag_ks= diag(ks);
                    
                    if nargout > 1
                        v= L\k;
                        sigma2_y=diag_ks-diag(v'*v);
                    end
                    if nargout > 4
                        Sigma2_y =  ks- v'*v;  % full covariance matrix (see 2.24 Rasmussen and Williams), = ks - (k'/K)*k;
                        if isequal(xtrain,xtest)
                            Sigma2_y = (Sigma2_y + Sigma2_y')/2;
                        end
                    end
                    
                    if nargout > 2
                        %koverK = (k'*invK); % (k'/K);
                        dmu_dx = squeeze(mtimesx(dk_dx, 'T', K\(ytrain-prior_mean_tr(:)))); %dprior_mean_test_dx +
                        dsigma2_dx = 2*squeeze(ddiagks_dx) -2*squeeze(sum(mtimesx(dk_dx, 'T', invK).*k', 2));
                        
                        dSigma2_dx =  dks_dx - 2*mtimesx(mtimesx(dk_dx, 'T', invK),k);
                    end
                    output1 = mu_y;
                else
                    output1 = post;
                end
            end
            
            sigma2_y(sigma2_y<0)=0;
        end
        
        %%
        function [negL, dnegL] = negloglike(model, theta, x_tr, y_tr)
            
            N = numel(y_tr);
            % covariance and derivative
            if nargout>1
                [K, dK] = model.kernelfun(theta.cov, x_tr,x_tr, true, model.regularization);
            else
                K = model.kernelfun(theta.cov, x_tr,x_tr, true, model.regularization);
            end
            
            if nargout>1
                [prior_mean, dprior_mean_dtheta]=model.meanfun(x_tr, theta.mean);
            else
                prior_mean= model.meanfun(x_tr, theta.mean);
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
        end
        
                       
        function [mu_y,  dmuy_dx] = to_maximize_mean_GP(theta, xtrain, ytrain, x, kernelfun, post)
            
            if any(isnan(x(:)))
                error('x is NaN')
            end            
            [mu_y, sigma2_y, dmu_dx] =  model.prediction(theta, xtrain, ytrain, x, post);
                       
            mu_y = -mu_y;
            dmuy_dx= -squeeze(dmu_dx);            
        end        
    end
end
