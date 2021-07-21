function [mu_c,  mu_y, sigma2_y, Sigma2_y, dmuc_dx, dmuy_dx, dsigma2y_dx, dSigma2y_dx, var_muc, dvar_muc_dx, post] =  prediction_bin(theta, xtrain, ctrain, xtest, kernelfun, varargin)
opts = namevaluepairtostruct(struct( ...
    'tol', 1e-7, ...
    'MaxIt', 1e3,...
    'post', [], ...
    'modeltype', 'exp_prop' ...
    ), varargin);

UNPACK_STRUCT(opts, false)

if isempty(post)
    clear('post')
end
if size(xtrain,1) ~= size(xtest,1)
    error("Dimensions of test and training sets are not consistent")
end

ctrain = ctrain(:);
if any(isnan([xtest(:);xtrain(:)]))
    error('x is NaN')
end
[d, n] = size(xtrain);

if any(ctrain ~= 1 & ctrain~= 0)
    error('c should be +1 or 0')
end
% kernel functions
if ~exist('post','var')
    K = kernelfun(theta, xtrain, xtrain,true);
    K =(K+K')/2;
else
    K = post.K;
end

if nargout>4
    [k,~,dk_dx] = kernelfun(theta, xtrain, xtest, false);
else
    k = kernelfun(theta, xtrain, xtest, false);
end
ks = kernelfun(theta, xtest, xtest,false);

if strcmp(modeltype, 'laplace')
    if ~exist('post','var')
        % find y= maximum of p(y_tr|ctrain, xtrain), (19.5.19 in Barber book)
        ystar = siteparams(ctrain, K, 'tol', tol, 'MaxIt', MaxIt);
        % compute 'noise' matrix (Eq 19.5.17 in Barber book)
        invD = diag(1./(logistic(ystar).*(1-logistic(ystar))+1e-12));

        post.ystar = ystar;
        post.invD = invD;
        post.invKD = inv(K + invD);
    else
        ystar = post.ystar;
        invD = post.invD;
        invKD =post.invKD;
    end
    % mean of p(y|x, ctrain, xtrain) (Eq 19.5.24 in Barber book)
    mu_y = k'*(ctrain - logistic(ystar));


    % standard deviation of p(y|x, ctrain, xtrain) (Eq 19.5.26 in Barber book)
    %     sigma2_y = diag(ks)-sum((k'/(K + invD)).*k', 2);
    sigma2_y = diag(ks)-sum((k'*invKD).*k', 2);

    if nargout>=4
        Sigma2_y = ks-k'*invKD*k;
        if isequal(xtrain,xtest)
            Sigma2_y = (Sigma2_y + Sigma2_y')/2;
        end
    end
    % mean of p(c|x, ctrain, xtrain) (Eq 19.5.27 in Barber book)
    nu = sqrt(pi)/4;
    mu_c = 0.5+0.5*erf(nu.*mu_y./sqrt(1+2*nu*sigma2_y));

    post.latent_y_mean= ystar;
    post.latent_y_cov = inv(inv(K) + inv(invD));
elseif strcmp(modeltype, 'exp_prop')

    if ~exist('post','var')
        [nu_tilde, tau_tilde, ~, ~, mu, Sigma] = siteparams(ctrain, K, 'modeltype', 'exp_prop', 'tol', tol);

        invS = diag(1./tau_tilde);
        invKS = inv(K+invS);

        post.invS = invS;
        post.invKS = invKS;
        post.nu_tilde = nu_tilde;
        post.tau_tilde = tau_tilde;
        post.Sigma =Sigma;
        post.mu = mu;
    else
        invS =  post.invS;
        invKS =  post.invKS;
        nu_tilde =  post.nu_tilde;
        tau_tilde =  post.tau_tilde;
        Sigma = post.Sigma;
        mu =  post.mu;
    end

    %     mu_y = (k'/(K+invS))*invS*nu_tilde;
    mu_y = (k'*invKS)*invS*nu_tilde;
    %sigma2_y = diag(ks)- sum((k'/(K+invS)).*k', 2);
    sigma2_y = diag(ks)- sum((k'*invKS).*k', 2);

    if nargout>=4
        Sigma2_y = ks - k'*invKS*k;
        if isequal(xtrain,xtest)
            Sigma2_y = (Sigma2_y + Sigma2_y')/2;
        end
    end

    mu_c = normcdf(mu_y./sqrt(1+sigma2_y));

    if nargout  == 11
        post.latent_y_cov = Sigma; %inv(inv(K) + invS);
        post.latent_y_mean= mu; %  post.latent_y_cov*invS*nu_tilde;
    end
else
    error('modeltype must be ''laplace'' or ''exp_prop''')
end


if nargout>4 && nargout~=8
    if strcmp(modeltype, 'laplace')

        % derivative of mu_c with respect to mu_y and sigma_y
        arg = nu.*mu_y.*(1+2*nu*sigma_y.^2).^(-0.5);

        darg_muy = nu./sqrt(1+2*nu*sigma_y.^2);

        darg_sigmay = -2*nu^2.*mu_y.*sigma_y.*(1 + 2*nu*sigma_y.^2).^(-1.5);

        dmuc_dmuy = exp(-arg.^2)./sqrt(pi).*darg_muy;

        dmuc_dsigmay = exp(-arg.^2)./sqrt(pi).*darg_sigmay;

        % derivative of mu_y with respect to x
        dmuy_dx = mtimesx(dk_dx, 'T',ctrain-logistic(ystar));

        % derivative of sigma_y with respect to x
        dsigma2y_dx = -sum(mtimesx(dk_dx, 'T', invKD).*k', 2)./sigma_y;

        % derivative of mu_c with respect to x
        dmuc_dx = dmuc_dsigmay.*dsigma2y_dx + dmuc_dmuy.*dmuy_dx;

    elseif strcmp(modeltype, 'exp_prop')

        % derivative of mu_c with respect to mu_y and sigma_y
        arg = normpdf(mu_y./sqrt(1+sigma2_y));

        dmuc_dmuy = arg./sqrt(1+sigma2_y);

        dmuc_dsigma2y = -0.5*arg.*mu_y./(1+sigma2_y).^(1.5);

        % derivative of mu_y with respect to x
        %dmuy_dx = mtimesx(dk_dx, 'T', (K+invS)\(invS*nu_tilde));
        dmuy_dx = mtimesx(dk_dx, 'T', invKS*(invS*nu_tilde));

        dsigma2y_dx = -2*sum(mtimesx(dk_dx, 'T', invKS).*k', 2);
        %                 dsigma2y_dx = -2*sum(mtimesx(dk_dx, 'T', pinv(K + invS)).*k', 2);

        % derivative of mu_c with respect to x
        dmuc_dx = dmuc_dsigma2y.*dsigma2y_dx + dmuc_dmuy.*dmuy_dx;

    else
        error('modeltype must be ''laplace'' or ''exp_prop''')
    end
end

if nargout >= 8
    %     tfn_output = NaN(numel(mu_c),1);
    dTdx = NaN(numel(mu_c),d);
    dvar_muc_dx =  NaN(numel(mu_c),d);
    dhdx = NaN(numel(mu_c),d);
    dadx = NaN(numel(mu_c),d);
    h = mu_y./sqrt(1+sigma2_y);
    a = 1./sqrt(1+2*sigma2_y);

    [tfn_output, dTdh, dTda] = tfn(h, a);

    for i= 1:numel(mu_c)
        dmcdx = squeeze(dmuc_dx(i,1,i,:))';
        dmydx = squeeze(dmuy_dx(i,1,i,:))';
        ds2dx = squeeze(dsigma2y_dx(i,1,i,:))';
        dhdx(i,:) = dmydx./sqrt(1+sigma2_y(i)) - 0.5*mu_y(i)*ds2dx.*((1+sigma2_y(i)).^(-3/2));
        dadx(i,:) = -ds2dx*((1+2*sigma2_y(i)).^(-1.5));

        dTdx(i,:) = dTdh(i)*dhdx(i,:) +dTda(i)*dadx(i,:);
        dvar_muc_dx(i,:) = dmcdx-2*mu_c(i)*dmcdx-2*dTdx(i,:);
    end
    var_muc = (mu_c - 2*tfn_output) - mu_c.^2;

    var_muc(sigma2_y==0) = 0;
end
dSigma2y_dx = [];

post.K = K;

return
