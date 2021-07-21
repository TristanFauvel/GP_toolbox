% This is an example of GP regression in one dimensional space

clear all;
close all;

addpath(genpath('../../GP_toolbox'))
addpath(genpath('../../BO_toolbox'))

%% Data definition
% Discretize space 
n=1000; 
x = linspace(-10,10,n);

% Function we want to sample from : 
y = 3*x;
y(y<0) = 0;

y = y+ randn(1,n);
% 
% Number of points in the training set
ntr =9;

% Seed (for reproducibility)
rng(6)

% Sample training data
i_tr= randsample(n,ntr);
x_tr = x(:,i_tr);
y_tr = y(:, i_tr);

% Choose test data
x_test = x;
y_test = y;

%% Regression

% Prior mean of the gaussian process

meanfunc = {@meanConst}; hyp.mean =0;
% covfunc = {@covPoly ,2};
covfunc = {@covNNone};
hyp.cov = [0;0];% log([1 1 1 1 1 1 1 1 1 1 1])';

feval(covfunc{:})

% covPoly(mode,par,d,hyp,x,z)
 
likfunc = @likGauss;
sn = 1; hyp.lik = log(sn);
% m = hyp.cov; V = eye(2);
% pG = {@priorGaussMulti ,m,V};
% prior.cov = pG;
inf_func = {@infExact};


K = feval(covfunc{:}, hyp.cov, x');
mu = feval(meanfunc{:}, hyp.mean, x');
sample  = mvnrnd(mu,K);

figure()
plot(x, sample)

% inf_func = {@infPrior,@infExact, prior};
[a b c d]  = gp(hyp, inf_func, meanfunc, covfunc, likfunc, x_tr', y_tr', x_test');

figure()
plot(x_test, c)

sample_x = x;
sample_y = mvnrnd(meanfun(sample_x,theta.mean), kernelfun(theta.cov,sample_x, sample_x, false));
figure()
scatter(sample_x, sample_y)

% Compute the posterior distribution

x_test_norm = x_test;
x_norm = x_tr;
mean_y = mean(y_tr);
y_norm = y_tr ; %- mean_y;

[mu_y, sigma2_y]= prediction(theta, x_norm, y_norm, x_test_norm, kernelfun, meanfun);
%mu_y = mu_y + mean_y;

graphics_style_paper;
% Plot the result
fig = figure();
fig.Color =  [1 1 1];
% plot(x,y); hold on;
plot(x,mu_y,'r'); hold on;
scatter(x_tr, y_tr,markersize, 'k', 'filled') ; hold on;
errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
box off

% Uncomment this if you want to use multistart optimization to estimate the hyperparameters using type-2 maximum likelihood.

%% Global optimization of hyperparameters
% gs = GlobalSearch;
% gs.Display='iter';
% lb= -100*ones(1,numel(hyp));
% ub= 100*ones(1,numel(hyp));
% opts = optimoptions(@fmincon,'Algorithm','sqp');
% myobjfun = @(hyp) minimize_negloglike(hyp, x_tr, y_tr, kernelfun, meanfun, ncov_hyp, nmean_hyp);
% problem = createOptimProblem('fmincon','objective', myobjfun,'x0',hyp, 'lb',lb,'ub',ub, 'options', opts);
% ms = MultiStart;
% hyp = run(ms,problem,200);
% theta.cov = hyp(1:ncov_hyp);
% theta.mean = hyp(ncov_hyp+1:ncov_hyp+nmean_hyp);


%% Local optimization of hyperparameters
update = 'all';
options=[];
hyp=[theta.cov, theta.mean];
hyp = minFunc(@(hyp)minimize_negloglike(hyp, x_norm, y_norm, kernelfun, meanfun, ncov_hyp, nmean_hyp, update), hyp, options);
theta.cov = hyp(1:ncov_hyp);
theta.mean = hyp(ncov_hyp+1:ncov_hyp+nmean_hyp);

%% Prediction with the new hyperparameters
[mu_y, sigma2_y]= prediction(theta, x_norm, y_norm, x_test_norm, kernelfun, meanfun);
mu_y = mu_y + mean_y;

fig = figure();
fig.Color =  [1 1 1];
plot(x,y); hold on;
scatter(x_tr, y_tr) ; hold on;
errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
box off

