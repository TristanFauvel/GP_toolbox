function [sample_g, dsample_g_dx] = sample_GP_precomputed_features(theta,phi, dphi_dx, xtrain, ytrain, kernelname, decoupled_bases, kernelfun)
%% SSGP : Method based on the Sparse-Spectrum GP, Lazaro-Gredilla 2010
%% RRGP: Method based on the Reduced-Rank GP, Solin 2019
%D = size(xtrain,1); %dimension

% xtrain : (2*dimension)*n
%ytrain : n x 1

Phi = phi(xtrain); % n x m
nfeatures = size(Phi,2);

if decoupled_bases
    noise = 0;
    w =randn(nfeatures,1);
    sample_prior = @(x) (phi(x)*w)';
    if contains(kernelname, 'wnoise')
        sigma = sqrt(exp(theta(end)));
        noise =  sigma*randn(N,1);
    end
     if isempty(ytrain)
        sample_g = @(x) sample_prior(x);
        dsample_g_dx = @(x) dprior_dx(x,w, dphi_dx)';        
    else
        K = kernelfun(theta,xtrain,xtrain, 1);
        v =  (K\(ytrain(:) - sample_prior(xtrain)'+noise))';
        update =  @(x) v*kernelfun(theta,xtrain,x, 0);
        sample_g = @(x) sample_prior(x) + update(x);
        dsample_g_dx = @(x) dprior_dx(x,w, dphi_dx)' + dupdate_dx(x, v, theta, xtrain, kernelfun)';
     end
else
    sig = 1e-9;
    A=Phi'*Phi+sig*eye(nfeatures); % A = nugget_regularization(A);
    
    theta = mvnrnd(A\(Phi'*ytrain), sig*inv(A))';
    

    %mu = Phi'*Phi+sig*eye(nfeatures)\(Phi'*ytrain);
    
%    theta = mvnrnd(mu, Sigma);

%Hernandez lobato: method: 
%       theta = sample_theta(nfeatures, Phi', sig, xtrain',ytrain);

sample_g= @(x) phi(x)*theta; % + meanfun(theta, x');
    dsample_g_dx = @(x) dphi_dx(x)*theta; % + dmean_dx(meanfun, theta, x);
    
end
end
function dudx = dupdate_dx(x, v, theta, xtrain, kernelfun)
[~, ~, dkdx]= kernelfun(theta,xtrain,x, 0);
dkdx = squeeze(dkdx);
dudx =mtimesx(v,dkdx);
end

function dpdx = dprior_dx(x,w, dphi_dx)
dpdx = w'*dphi_dx(x);
end

function dmdx = dmean_dx(meanfun,theta, x)
[~,~, dmdx] = meanfun(theta, x');
end

function theta= sample_theta(nFeatures, Z, sigma0, xx, yy)
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
