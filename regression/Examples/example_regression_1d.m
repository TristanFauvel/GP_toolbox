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
objective = objective();
objective = @(x) objective.do_eval(x);

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
xtrain = x(:,i_tr);
ytrain = y(:, i_tr);

% Choose test data
x_test = x;
y_test = y;

%% Regression

% Prior mean of the gaussian process
meanfun= @constant_mean;

%Choice of the kernel : ARD kernel
kernelfun = @ARD_kernelfun;
kernelname = 'ARD';

% Initialize hyperparameters
% theta.cov = [1,2,1];
D = 1;
regularization = 'nugget';
type = 'regression';
hyps.ncov_hyp =2; % number of hyperparameters for the covariance function
hyps.nmean_hyp =1; % number of hyperparameters for the mean function
hyps.hyp_lb = -10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
hyps.hyp_ub = 10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
lb = xbounds(1); ub = xbounds(2); %bounds of the space
model = gp_regression_model(D, meanfun, kernelfun, regularization, hyps, lb, ub, kernelname);

% Compute the posterior distribution
theta.cov = [1,2];
theta.mean = zeros(hyps.nmean_hyp,1);

x_test_norm =(x_test- xbounds(1))./(xbounds(2)-xbounds(1));
x_norm = (xtrain- xbounds(1))./(xbounds(2)-xbounds(1));
mean_y = mean(ytrain);
y_norm = ytrain - mean_y;

[mu_y, sigma2_y]= model.prediction(theta, x_norm, y_norm, x_test_norm, []);
mu_y = mu_y + mean_y;

graphics_style_paper;
% Plot the result
fig = figure();
fig.Color =  background_color;
plot(x,y); hold on;
plot(x,mu_y,'r'); hold on;
scatter(xtrain, ytrain,markersize, 'k', 'filled') ; hold on;
errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
box off


%% Optimization of hyperparameters
update = 'all';
options=[];

hyp = model.model_selection(x_norm, y_norm, theta,  update);

%% Prediction with the new hyperparameters
[mu_y, sigma2_y]= model.prediction(theta, x_norm, y_norm, x_test_norm, []);
mu_y = mu_y + mean_y;

fig = figure();
fig.Color =  background_color;
plot(x,y); hold on;
scatter(xtrain, ytrain, markersize, 'k', 'filled') ; hold on;
errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
title('Hyperparameters inferred using ML');
box off


