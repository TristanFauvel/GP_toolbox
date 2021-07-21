% This is an example of GP regression in one dimensional space

clear all;
close all;

addpath(genpath('../../GP_toolbox'))
addpath(genpath('../../BO_toolbox'))

objective ='grlee12';
if strcmp(objective, 'forretal08')
    xbounds = [0,1];
    d=1; %dimension of the input space
elseif strcmp(objective, 'grlee12')
    xbounds = [0.5,2.5];
    d=1; %dimension of the input space
end
objective = str2func(objective);
%% Data definition
% Discretize space
n=1000;
x = linspace(xbounds(1), xbounds(2),n);

% Function we want to sample from :
y = objective(x);

% Number of points in the training set
ntr =250;

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
meanfun= @constant_mean;

%Choice of the kernel : ARD kernel
kernelfun = @ARD_kernelfun;
ncov_hyp=2; % number of hyperparameters for the covariance function

% kernelfun = @Rational_Quadratic_kernelfun;
% ncov_hyp=3; % number of hyperparameters for the covariance function

nmean_hyp=1; % number of hyperparameters for the mean function

% Initialize hyperparameters
% theta.cov = [1,2,1];
theta.cov = [1,2];
theta.mean = zeros(nmean_hyp,1);

% Compute the posterior distribution

x_test_norm =(x_test- xbounds(1))./(xbounds(2)-xbounds(1));
x_norm = (x_tr- xbounds(1))./(xbounds(2)-xbounds(1));
mean_y = mean(y_tr);
y_norm = y_tr - mean_y;

[mu_y, sigma2_y]= prediction(theta, x_norm, y_norm, x_test_norm, kernelfun, meanfun);
mu_y = mu_y + mean_y;

graphics_style_paper;
% Plot the result
fig = figure();
fig.Color =  [1 1 1];
plot(x,y); hold on;
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
scatter(x_tr, y_tr, markersize, 'k', 'filled') ; hold on;
errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
title('Hyperparameters inferred using ML');
box off


