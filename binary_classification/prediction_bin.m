function [output1,  mu_y, sigma2_y, Sigma2_y, dmuc_dx, dmuy_dx, dsigma2y_dx, dSigma2y_dx, var_muc, dvar_muc_dx,post] =  prediction_bin(theta, xtrain, ctrain, xtest, kernelfun, modeltype, post, regularization)
% mu_c corresponds to P(x1 > x2)

if ~isempty(xtest) && size(xtrain,1) ~= size(xtest,1)
    error("Dimensions of test and training sets are not consistent")
end

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

    if nargout>4
        [k, ~, dk_dx] = kernelfun(theta, xtrain, xtest, false, 'false');
        post.dk_dx = dk_dx;
    else
        k = kernelfun(theta, xtrain, xtest, false, 'false');
    end
    post.k = k;
end

if ~isempty(xtest)
    if nargout>4
        [ks, ~, dks_dx] = kernelfun(theta, xtest, xtest, false, 'false');  %%%% Careful : normally here it's 'false', but I put 'true'in order to get a correct computation of the derivative, only ok if noise = 0
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
        ystar = siteparams(ctrain, K, 'tol', tol, 'MaxIt', MaxIt);
        % compute 'noise' matrix (Eq 19.5.17 in Barber book)
        invS = diag(1./(logistic(ystar).*(1-logistic(ystar))+1e-12));
        invKS = inv(K + invS);

        post.ystar = ystar;
        post.invKS = invKS;
        post.invS = invS;
    else
        ystar = post.ystar;
        invKS = post.invKS;
        invS = post.invS;
    end
    if ~isempty(xtest)

        % mean of p(y|x, c_tr, xtrain) (Eq 19.5.24 in Barber book)
        mu_y = k'*(ctrain - logistic(ystar));

        % standard deviation of p(y|x, c_tr, xtrain) (Eq 19.5.26 in Barber book)
        %sigma2_y = diag_ks-sum((k'/(K + invS)).*k', 2);
        sigma2_y = diag_ks-sum((k'*invKS).*k', 2);
        if nargout>=4
            %Sigma2_y = ks-k'*inv(K + invS)*k;
            Sigma2_y = ks-k'*invKS*k;
            if isequal(xtrain,xtest)
                Sigma2_y = (Sigma2_y + Sigma2_y')/2;
            end
        end
        % mean of p(c|x, c_tr, xtrain) (Eq 19.5.27 in Barber book)
        nu = sqrt(pi)/4;
        mu_c = 0.5+0.5*erf(nu.*mu_y./sqrt(1+2*nu*sigma2_y));
    end
elseif strcmp(modeltype, 'exp_prop')
    if comp_post
        [nu_tilde, tau_tilde] = siteparams(ctrain, K, 'modeltype', 'exp_prop', 'tol', tol);

        invS = diag(1./tau_tilde);
        invKS = inv(K+invS);

        post.nu_tilde = nu_tilde;
        post.invS = invS;
        post.invKS = invKS;
        post.tau_tilde = tau_tilde;
    else
        invKS = post.invKS;
        nu_tilde = post.nu_tilde;
        invS = post.invS;
    end
    if ~isempty(xtest)
        mu_y = (k'*invKS)*invS*nu_tilde;
        sigma2_y = diag_ks- sum((k'*invKS).*k', 2);

        if nargout>=4
            Sigma2_y = ks - k'*invKS*k;
            if isequal(xtrain,xtest)
                Sigma2_y = (Sigma2_y + Sigma2_y')/2;
            end
        end
        mu_c = normcdf(mu_y./sqrt(1+sigma2_y));
    end
else
    error('modeltype must be ''laplace'' or ''exp_prop''')
end

if ~isempty(xtest)
    if nargout>4
        if strcmp(modeltype, 'laplace')
            sigma_y = sqrt(sigma2_y);
            % derivative of mu_c with respect to mu_y and sigma_y
            arg = nu.*mu_y.*(1+2*nu*sigma_y.^2).^(-0.5);

            darg_muy = nu./sqrt(1+2*nu*sigma_y.^2);

            %darg_sigmay = -2*nu^2.*mu_y.*sigma_y.*(1 + 2*nu*sigma_y.^2).^(-1.5);
            darg_sigma2y = -nu^2.*mu_y.*(1 + 2*nu*sigma_y.^2).^(-1.5);

            dmuc_dmuy = exp(-arg.^2)./sqrt(pi).*darg_muy;

            %dmuc_dsigmay = exp(-arg.^2)./sqrt(pi).*darg_sigmay;
            dmuc_dsigma2y = exp(-arg.^2)./sqrt(pi).*darg_sigma2y;

            % derivative of mu_y with respect to x
            dmuy_dx = squeeze(mtimesx(dk_dx, 'T',ctrain-logistic(ystar)));

            % derivative of sigma_y with respect to x
            %dsigma2y_dx = -2*ddiagks_dx -sum(mtimesx(dk_dx, 'T', invKS).*k', 2)./sigma_y;

        elseif strcmp(modeltype, 'exp_prop')

            % derivative of mu_c with respect to mu_y and sigma_y
            arg = normpdf(mu_y./sqrt(1+sigma2_y));

            dmuc_dmuy = arg./sqrt(1+sigma2_y);

            dmuc_dsigma2y = -0.5*arg.*mu_y./((1+sigma2_y).^(1.5));
            % derivative of mu_y with respect to x
            %         dmuy_dx = squeeze(mtimesx(dk_dx, 'T', (K+invS)\(invS*nu_tilde)));
            dmuy_dx = squeeze(mtimesx(dk_dx, 'T', post.invKS*(invS*nu_tilde)));
        else
            error('modeltype must be ''laplace'' or ''exp_prop''')
        end
           dsigma2y_dx = 2*squeeze(ddiagks_dx) -2*squeeze(sum(mtimesx(dk_dx, 'T', invKS).*k', 2));
            dSigma2y_dx = squeeze(dks_dx) - 2*mtimesx(mtimesx(dk_dx, 'T', invKS),k);

                    % derivative of mu_c with respect to x
            dmuc_dx = dmuc_dsigma2y.*dsigma2y_dx + dmuc_dmuy.*dmuy_dx;

    end


    if nargout >= 9 
        if strcmp(modeltype, 'exp_prop')
        h = mu_y./sqrt(1+sigma2_y);
        a = 1./sqrt(1+2*sigma2_y);

        [tfn_output, dTdh, dTda] = tfn(h, a);

        indexat = @(expr, varargin) expr(varargin{:});
        dmcdx = indexat(reshape(dmuc_dx, size(xtest,2)^2, D),1:size(xtest,2)+1:(size(xtest,2)^2), 1:D);
        dmydx =  indexat(reshape(dmuy_dx, size(xtest,2)^2, D),1:size(xtest,2)+1:(size(xtest,2)^2), 1:D);
        ds2dx = indexat(reshape(dsigma2y_dx, size(xtest,2)^2, D),1:size(xtest,2)+1:(size(xtest,2)^2), 1:D);
        dhdx = dmydx./sqrt(1+sigma2_y) - 0.5*mu_y.*ds2dx.*((1+sigma2_y).^(-3/2));
        dadx =  -ds2dx.*((1+2*sigma2_y).^(-1.5));
        dTdx = dTdh.*dhdx +dTda.*dadx;
        dvar_muc_dx = dmcdx-2*mu_c.*dmcdx-2*dTdx;

        var_muc = (mu_c - 2*tfn_output) - mu_c.^2;

        var_muc(sigma2_y==0) = 0;
        else
            var_muc = [];
        end
    end
end

if isempty(xtest)
    output1 = post;
else
    output1 = mu_c;
end
return
