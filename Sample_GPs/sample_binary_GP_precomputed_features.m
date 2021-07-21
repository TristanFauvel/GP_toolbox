function [sample_g, dsample_g_dx, decomposition] = sample_binary_GP_precomputed_features(phi, dphi_dx, xtrain, ctrain, theta,kernelfun, decoupled_bases, modeltype, post)
%% SSGP : Method based on the Sparse-Spectrum GP, Lazaro-Gredilla 2010
%% RRGP: Method based on the Reduced-Rank GP, Solin 2019

% xtrain : (2*dimension)*n
%ytrain : n x 1

if isempty(post)
[~,  mu_y, ~, Sigma2_y,~,~,~,~,~,~,post] =prediction_bin(theta, xtrain, ctrain, xtrain, kernelfun, 'modeltype', modeltype, 'post', post);
else
    [~,  mu_y, ~, Sigma2_y] =prediction_bin(theta, xtrain, ctrain, xtrain, kernelfun, 'modeltype', modeltype, 'post', post);
end

Sigma2_y = nugget_regularization(Sigma2_y);


Phi = phi(xtrain); % n x m

ytrain =  mvnrnd(mu_y(:), Sigma2_y)';

if decoupled_bases
    nfeatures = size(Phi,2);
    w =randn(nfeatures,1);
    sample_prior = @(x) (phi(x)*w)';

    K = kernelfun(theta,xtrain,xtrain, 1); %ntr x ntr

    v =  (K\(ytrain - sample_prior(xtrain)'))'; %1 x ntr
    update =  @(x) v*kernelfun(theta,xtrain,x, 0); % 1 x ntest
    sample_g = @(x) sample_prior(x) + update(x); % 1 x ntest
    dsample_g_dx = @(x) dprior_dx(x, w, dphi_dx)' + dupdate_dx(x, v, theta, xtrain, kernelfun)'; % D x ntest x ntest
    decomposition.sample_prior = sample_prior;
    decomposition.update = update;
else
    A=Phi'*Phi +1e-12*eye(size(Phi,2)); % A = nugget_regularization(A);
    theta = A\(Phi'*(ytrain));
    sample_g= @(x) phi(x)*theta;
    dsample_g_dx = @(x) dphi_dx(x)*theta;
    decomposition = [];
end
return

end
function dudx = dupdate_dx(x, v, theta, xtrain, kernelfun)
if size(x,2) >1
    error('Derivative only implemented for size(x,2) == 1')
end
[k, ~, dkdx]= kernelfun(theta,xtrain,x, 0);
dkdx = squeeze(dkdx);
dudx =mtimesx(v,dkdx);
end

function dpdx = dprior_dx(x,w, dphi_dx)
dpdx = w'*dphi_dx(x); %1 x D  does not work with SSGP and ARD, but works
% with Matern
% dpdx = w'*dphi_dx(x)'; %1 x D  works with SSGP and ARD, but not with Matern and RRGP

end
