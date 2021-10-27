function [sample_g, dsample_g_dx, sample_prior, update] = sample_value_GP_precomputed_features(approximation, theta, xtrain_norm, ctrain, post);(phi_pref, dphi_pref_dx, x_data, c_data, x0, hyp,kernelfun, varargin)
opts = namevaluepairtostruct(struct( ...
    'decoupled_bases', 0, ...
    'modeltype', 'exp_prop', ...
    'post', [] ...
    ), varargin);
UNPACK_STRUCT(opts, false)

%% SSGP : Method based on the Sparse-Spectrum GP, Lazaro-Gredilla 2010
%% RRGP: Method based on the Reduced-Rank GP, Solin 2019
D = size(x_data,1)/2; %dimension
% x_data : (2*dimension)*n
%y_data : n x 1

Phi = phi_pref(x_data); % n x m
    nfeatures = size(Phi,2);

[mu_c,  mu_y, sigma2_y, Sigma2_y, ~, ~, ~, ~, ~, ~, post] =model.prediction(hyp, x_data, c_data, x_data, post);
Sigma2_y = (Sigma2_y + Sigma2_y')/2;

y_data =  mvnrnd(mu_y(:), Sigma2_y)';

sample_prior = [];
update = [];
if decoupled_bases
    nfeatures = size(Phi,2);
    w =randn(nfeatures,1);
    sample_prior = @(x) (phi_pref(x)*w)';
        
    v =  (post.K\(y_data - sample_prior(x_data)'))';
    update =  @(x) v*kernelfun(hyp.cov,x_data,x, 0);
    sample_g = @(x) sample_prior([x;x0.*ones(D,size(x,2))]) + update([x;x0.*ones(D,size(x,2))]);
    dsample_g_dx = @(x) dprior_dx(x, x0, D, w, dphi_pref_dx) + dupdate_dx(x, x0, D, v, hyp, x_data, kernelfun);
else
    sig = 1e-12;
    A=Phi'*Phi +sig*eye(size(Phi,2)); % A = nugget_regularization(A);
    theta = sample_theta(nfeatures, Phi', sig, x_data',y_data);
    sample_g= @(x) phi_pref([x;x0.*ones(D,size(x,2))])*theta;
    dsample_g_dx = @(x) dphi_pref_dx([x;x0.*ones(D,size(x,2))])*theta;
end
return

end
function dudx = dupdate_dx(x, x0, D, v, hyp, x_data, kernelfun)
if size(x,2) >1
    error('Derivative only implemented for size(x,2) == 1')
end
[~, ~, dkdx]= kernelfun(hyp,x_data,[x;x0.*ones(D,size(x,2))], 0);
dkdx = squeeze(dkdx); % ntr*d
dkdx = dkdx(:,1:D);
dudx =mtimesx(v,dkdx)';
end

function dpdx = dprior_dx(x, x0, D,w, dphi_pref_dx)
dpdx =dphi_pref_dx([x;x0.*ones(D,size(x,2))])*w; %D x n
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
