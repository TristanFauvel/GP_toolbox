clear all;
close all;
addpath(genpath('/home/tristan/Desktop/GP_toolbox'))
n=100;

% objective = 'forretal08';
% objective_fun  = @forretal08;

objective = 'grlee12';
objective_fun  = @grlee12;

meanfun= @constant_mean;

%Gaussian kernel
ncov_hyp=3;
nmean_hyp=1;
kernelfun= @Gaussian_kernelfun_wnoise;

theta.mean = zeros(nmean_hyp,1);

theta.cov = rand(ncov_hyp,1);

d=2;
if strcmp(objective, 'forretal08')
    x = linspace(0,1, n);
    d=1; %dimension of the input space
    theta.cov  = [2.6933, 8.6930,-29.2431]';
elseif strcmp(objective, 'grlee12')
    x = linspace(0,1, n);
    d=1; %dimension of the input space
    theta.cov  = [4.0826, 6.3708, -29.2431]';
elseif strcmp(objective, 'levy')
    x1range= linspace(-10,10, n);
    x2range=x1range;
    [x1, x2] = ndgrid(x1range);
    x=  [x2(:), x1(:)]';
    theta.cov  = [-0.5387, 5.6751,-88.9788]';
elseif strcmp(objective, 'goldpr')
    x1range= linspace(-2,2, n);
    x2range=x1range;
    [x1, x2] = ndgrid(x1range);
    x=  [x2(:), x1(:)]';
    theta.cov  = [-1 , 0, -10]'; %[-1 ; 0; -10]
elseif strcmp(objective, 'camel6')
    x1range = linspace(-3,3, n);
    x2range = linspace(-2,2, n);
    [x1, x2] = ndgrid(x1range,x2range);
    x=  [x2(:), x1(:)]';
    theta = [ -0.5642, 5.7548,-98.9687]';
end


y = objective_fun(x);


rng(10)

ntr = 5;
i_tr= randsample(n,ntr);
% i_tr= 1:10:n^d;
ntr=numel(i_tr);
xtrain = x(:,i_tr);
y_tr = y(:, i_tr)+ randn(1,ntr);

x_test = x;
y_test = y;


[mu_y, sigma2_y]= prediction(theta, xtrain, y_tr, x_test, kernelfun, meanfun);

if d==1
    fig = figure();
    fig.Color =  background_color;
    plot(x,y); hold on;
    scatter(xtrain, y_tr) ; hold on;
    errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
    pbaspect([1,1,1])
elseif d==2
    fig = figure();
    fig.Color =  background_color;
    subplot(1,2,1)
    imagesc(x1range, x2range, reshape(y, n,n));
    pbaspect([1,1,1])
    title('True function')
    subplot(1,2,2)
    imagesc(x1range, x2range, reshape(mu_y', n,n));
    pbaspect([1,1,1])
    title('Regression function')
end

update = 'cov';

%% Local optimization of hyperparameters
options=[];
hyp = minFunc(@(hyp)negloglike_zero_mean(hyp, xtrain, y_tr, kernelfun, meanfun, ncov_hyp, nmean_hyp, update), theta.cov, options);
theta.cov = hyp;

%% Prediction with the new hyperparameters
[mu_y, sigma2_y]= prediction(theta, xtrain, y_tr, x_test, kernelfun, meanfun);

if d==1
    fig = figure();
    fig.Color =  background_color;
    plot(x,y); hold on;
    scatter(xtrain, y_tr) ; hold on;
    errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
elseif d==2
    fig = figure();
    fig.Color =  background_color;
    subplot(1,2,1)
    imagesc(x1range, x2range, reshape(y, n,n));
    pbaspect([1,1,1])
    title('True function')
    subplot(1,2,2)
    imagesc(x1range, x2range, reshape(mu_y', n,n));
    pbaspect([1,1,1])
    title('Regression function')
end
% update = 'all';
% neg = @(hyp) minimize_negloglike(hyp, xtrain, y_tr, kernelfun, meanfun, ncov_hyp, nmean_hyp, update);
% result= squeeze(test_deriv(neg, hyp, 1e-12));
% [negL, dnegL] =minimize_negloglike(hyp, xtrain, y_tr, kernelfun, meanfun,ncov_hyp, nmean_hyp, update);
% R=sqrt((dnegL-squeeze(result)).^2);
% figure()
% plot(result); hold on;
% plot(dnegL); hold off;
%
% [prior_mean, dprior_mean_dtheta]=meanfun(xtrain, theta.mean);
% pm= @(theta) meanfun(xtrain, theta);
% result= squeeze(test_deriv(pm, theta.mean, 1e-12));
% figure()
% plot(result); hold on;
% plot(dprior_mean_dtheta); hold off;
%
%
%
% k=@(theta) kernelfun(theta, x_tr);
% result= squeeze(test_deriv(k, theta.cov, 1e-8));
% [K, dK] = kernelfun(theta.cov, x_tr);
% R=sqrt((dK(:)-result(:)).^2);
% figure()
% plot(result(:)); hold on;
% plot(dK(:)); hold off;
%
% k=@(x_test) kernelfun(theta.cov, xtrain, x_test);
% result= squeeze(test_deriv(k, x_test(1), 1e-8));
% [K, dK, dK_dx] = kernelfun(theta.cov, xtrain, x_test(1));
% R=sqrt((dK_dx(:)-result(:)).^2);
% max(R)
% figure()
% plot(result(:)); hold on;
% plot(dK_dx(:)); hold off;
%
