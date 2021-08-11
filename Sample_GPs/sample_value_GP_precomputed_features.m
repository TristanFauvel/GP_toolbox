function [sample_g, dsample_g_dx, decomposition] = sample_value_GP_precomputed_features(phi, dphi_dx, phi_pref, dphi_pref_dx, xtrain, ctrain, theta, model, approximation, post)
%% SSGP : Method based on the Sparse-Spectrum GP, Lazaro-Gredilla 2010
%% RRGP: Method based on the Reduced-Rank GP, Solin 2019
D = (size(xtrain,1))/2; %dimension
% xtrain : (2*dimension)*n
%y_data : n x 1
regularization = model.regularization;
condition = model.condition;
decoupled_bases = approximation.decoupled_bases;
base_kernelfun = model.base_kernelfun;
kernelfun = model.kernelfun;

if isempty(post)
    [~,  mu_y, ~, Sigma2_y,~,~,~,~,~,~,post] =prediction_bin(theta, xtrain, ctrain, xtrain, model, post);
else
    [~,  mu_y, ~, Sigma2_y] =prediction_bin(theta, xtrain, ctrain, xtrain, model, post);
end

Sigma2_y = nugget_regularization(Sigma2_y);
y_data =  mvnrnd(mu_y(:), Sigma2_y)';

Phi = phi_pref(xtrain); % n x m
nfeatures = size(Phi,2);

x0 = condition.x0;
if isfield(condition, 'y0')
    if decoupled_bases == 0
        error('You should use the decoupled bases method')
    end
    %base_kernelfun is the base kernel, without conditioning.
    w =randn(nfeatures,1);
    
    sample_prior = @(x) (phi(x)*w)';
    
    K = base_kernelfun(theta,condition.x0, condition.x0, 'true', regularization);
    
    v1 =  (K\(condition.y0 - sample_prior(condition.x0)'))';
    update_1 =  @(x) v1*base_kernelfun(theta,condition.x0,x, 0, regularization);
    
    cond_sample_prior = @(x) sample_prior(x) + update_1(x);
    cond_sample_prior_pref = @(x) cond_sample_prior(x(1:D,:)) - cond_sample_prior(x(D+1:end,:));
    
    K = post.K;
    v2 =  (K\(y_data - cond_sample_prior_pref(xtrain)'))';
    update_2 =  @(x) v2*kernelfun(theta,xtrain,x, 0, regularization);
    sample_g = @(x) cond_sample_prior(x) + update_2([x;condition.x0.*ones(D,size(x,2))]);
    dsample_g_dx = @(x) dcond_prior_dx(x, D, w, v1, dphi_dx, base_kernelfun, theta, condition.x0, regularization) + dupdate_dx(x, x0, D, v2, theta, xtrain, kernelfun, regularization); % Dxntest
    
    
    decomposition.update_1 = update_1;
    decomposition.update_2 = update_2;
    decomposition.sample_prior = sample_prior;
    decomposition.cond_sample_prior = cond_sample_prior;
    
else
    %     warning('There is no conditioning on the value function offset')
    if decoupled_bases
        w =randn(nfeatures,1);
        sample_prior = @(x) (phi_pref(x)*w)';
        
        v =  (post.K\(y_data - sample_prior(xtrain)'))';
        update =  @(x) v*kernelfun(theta,xtrain,x, 0, regularization);
        sample_g = @(x) sample_prior([x;x0.*ones(D,size(x,2))]) + update([x;x0.*ones(D,size(x,2))]);
        dsample_g_dx = @(x) dprior_dx(x, x0, D, w, dphi_pref_dx) + dupdate_dx(x, x0, D, v, theta, xtrain, kernelfun, regularization);
        decomposition.update = update;
        decomposition.sample_prior = sample_prior;
    else
        decomposition = [];
        sig = 1e-12;
        %A=Phi'*Phi +sig*eye(size(Phi,2)); % A = nugget_regularization(A);
        theta = sample_theta(nfeatures, Phi', sig, xtrain',y_data);
        sample_g= @(x) phi_pref([x;x0.*ones(D,size(x,2))])*theta;
        dsample_g_dx = @(x) dphi_pref_dx([x;x0.*ones(D,size(x,2))])*theta;
    end
end
end
function dudx = dupdate_dx(x, x0, D, v, theta, xtrain, kernelfun, regularization)
if size(x,2) >1
    error('Derivative only implemented for size(x,2) == 1')
end
[~, ~, dkdx]= kernelfun(theta,xtrain,[x;x0.*ones(D,size(x,2))], 0, regularization);
dkdx = squeeze(dkdx); % ntr*d
dkdx = dkdx(:,1:D);
dudx =mtimesx(v,dkdx)';
end

function dpdx = dprior_dx(x, x0, D,w, dphi_pref_dx)
dpdx =dphi_pref_dx([x;x0.*ones(D,size(x,2))])*w; %D x n
end

function dpdx = dcond_prior_dx(x, D, w, v1, dphi_dx, base_kernelfun, theta, x0, regularization)
if size(x,2) >1
    error('Derivative only implemented for size(x,2) == 1')
end
[~, ~, dkdx]= base_kernelfun(theta,x0,x, 0, regularization); % 1 x 1 x 1 x D
if ~ (isequal(size(dkdx), [1, 1, 1, D])  || isequal(size(dkdx), [1, 1]))
    error('Error in the kernel derivative dimensions')
end
dkdx = squeeze(dkdx)';
dpdx = dphi_dx(x)'*w + mtimesx(v1,dkdx)';
end

function theta= sample_theta(nFeatures, Z, sigma0, xx, yy) %From Hernandez-Lobato
% Draw the coefficient theta.
noise = randn(nFeatures, 1);
if (size(xx, 1) < nFeatures)
    % We adopt the formula $theta \sim \N(Z(Z'Z + \sigma^2 I)^{-1} y,
    % I-Z(Z'Z + \sigma^2 I)Z')$.
    Sigma = Z' * Z + sigma0 * eye(size(xx, 1));
    mu = Z*chol2invchol(Sigma)*yy;
    [U, D] = eig(Sigma);
    D = diag(D);
    R = (sqrt(D) .* (sqrt(D) + sqrt(sigma0))).^-1;
    theta = noise - (Z * (U * (R .* (U' * (Z' * noise))))) + mu;
else
    % $theta \sim \N((ZZ'/\sigma^2 + I)^{-1} Z y / \sigma^2,
    % (ZZ'/\sigma^2 + I)^{-1})$.
    Sigma = chol2invchol(Z*Z' / sigma0 + eye(nFeatures));
    mu = Sigma * Z * yy / sigma0;
    theta = mu + noise * chol(Sigma);
end
end
