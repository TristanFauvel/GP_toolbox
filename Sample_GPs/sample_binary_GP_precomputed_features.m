function [sample_g, dsample_g_dx, decomposition] = sample_binary_GP_precomputed_features(xtrain, ctrain, theta, model, approximation, post)
%% SSGP : Method based on the Sparse-Spectrum GP, Lazaro-Gredilla 2010
%% RRGP: Method based on the Reduced-Rank GP, Solin 2019

% xtrain : (2*dimension)*n
%ytrain : n x 1
decoupled_bases = approximation.decoupled_bases;
kernelfun = model.kernelfun;
regularization = model.regularization;
phi=approximation.phi;
dphi_dx = approximation.dphi_dx;

if isempty(post)
post =prediction_bin(theta, xtrain, ctrain, [], model, post);
end
[~,  mu_y, ~, Sigma2_y] =prediction_bin(theta, xtrain, ctrain, xtrain, model, post);


Sigma2_y = nugget_regularization(Sigma2_y);


Phi = phi(xtrain); % n x m

ytrain =  mvnrnd(mu_y(:), Sigma2_y)';

if decoupled_bases
    nfeatures = size(Phi,2);
    w =randn(nfeatures,1);
    sample_prior = @(x) (phi(x)*w)';

    K = kernelfun(theta,xtrain,xtrain, 1, regularization); %ntr x ntr

    v =  (K\(ytrain - sample_prior(xtrain)'))'; %1 x ntr
    update =  @(x) v*kernelfun(theta,xtrain,x, 0, regularization); % 1 x ntest
    sample_g = @(x) sample_prior(x) + update(x); % 1 x ntest
    dsample_g_dx = @(x) dprior_dx(x, w, dphi_dx) + dupdate_dx(x, v, theta, xtrain, kernelfun,regularization); % D x ntest x ntest
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

function dudx = dupdate_dx(x, v, theta, xtrain, kernelfun, regularization)
% if size(x,2) >1
%     error('Derivative only implemented for size(x,2) == 1')
% end
% [k, ~, dkdx]= kernelfun(theta,xtrain,x, 0, regularization);
% dkdx = squeeze(dkdx);
% dudx =mtimesx(v,dkdx);

dudx = NaN(size(x,1), size(x,2), size(x,2));
for i = 1:size(x,2)
    [k, ~, dkdx]= kernelfun(theta,xtrain,x(:,i), 0, regularization);
    dkdx = squeeze(dkdx);
    dudx(:,i,i) =mtimesx(v,dkdx);
end
if i==1
    dudx = dudx';
end
end

function dpdx = dprior_dx(x,w, dphi_dx)
dpdx = NaN(size(x,1),size(x,2),size(x,2));
for i = 1:size(x,2)
    dpdx(:,i,i) = w'*dphi_dx(x(:,i));
end

if i==1
    dpdx = dpdx';
end
end
