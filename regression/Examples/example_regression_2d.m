clear all;
close all;

  objective = camel6();
lb = objective.xbounds(:,1);
ub = objective.xbounds(:,2);
D = 2;
%%
rng(1)
addpath(genpath('/home/tfauvel/Desktop/optim_retina/GP_toolbox'))
n=100;
x1range = linspace(lb(1),ub(1),n);
x2range = linspace(lb(2),ub(2),n);
[p,q]= meshgrid(x1range,x2range);
x=[p(:),q(:)]';
y = objective.do_eval(x); % sin(x(1,:)).*sin(x(2,:));

ntr =400;
i_tr= randsample(n*n,ntr);
x_tr = x(:,i_tr);
y_tr = y(:, i_tr);

x_test = x;
y_test = y;

x_norm = (x_tr - lb)./(ub - lb);
x_test_norm = (x_test - lb)./(ub - lb);
mean_y = mean(y_tr);
y_norm = y_tr - mean_y;

meanfun= @constant_mean;
kernel_name = 'ARD';

if strcmp(kernel_name, 'ARD')
    %ARD kernel
    ncov_hyp=1+D;
    nmean_hyp=1;
    kernelfun = @ARD_kernelfun;
elseif strcmp(kernel_name, 'Gaussian_wnoise')
    %Gaussian kernel
    ncov_hyp=3;
    nmean_hyp=1;
    kernelfun= @Gaussian_kernelfun_wnoise;
end

theta.cov = 10*ones(ncov_hyp,1);%rand(ncov_hyp,1);
theta.mean = zeros(nmean_hyp,1);

regularization = 'nugget';
hyps.ncov_hyp =ncov_hyp; % number of hyperparameters for the covariance function
hyps.nmean_hyp =nmean_hyp; % number of hyperparameters for the mean function
hyps.hyp_lb = -10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
hyps.hyp_ub = 10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);

model = gp_regression_model(D, meanfun, kernelfun, regularization, hyps, lb, ub);


[mu_y, sigma2_y]= model.prediction(theta, x_norm, y_norm, x_test_norm,[]);
mu_y = mu_y + mean_y;

sigma2_y(sigma2_y<0) = 0;
figure()
subplot(3,1,1)
imagesc(x1range, x2range, reshape(y, n, n)); hold on;
scatter(x_tr(1,:), x_tr(2,:)) ; hold off;
pbaspect([1,1,1]);
subplot(3,1,2)
imagesc(x1range, x2range, reshape(mu_y, n, n)); hold on;
pbaspect([1,1,1]);
subplot(3,1,3)
imagesc(x1range, x2range, reshape(sqrt(sigma2_y), n, n)); hold on;
pbaspect([1,1,1]);



hyp=[theta.cov; theta.mean];
update = 'all';

theta = model.model_selection(x_norm, y_norm, update);

%% Prediction with the new hyperparameters
[mu_y, sigma2_y]= model.prediction(theta, x_norm, y_norm, x_test_norm, []);
mu_y = mu_y + mean_y;

figure()
subplot(3,1,1)
imagesc(x1range, x2range, reshape(y, n, n)); hold on;
scatter(x_tr(1,:), x_tr(2,:)) ; hold off;
subplot(3,1,2)
imagesc(x1range, x2range, reshape(mu_y, n, n)); hold on;
subplot(3,1,3)
imagesc(x1range, x2range, reshape(sqrt(sigma2_y), n, n)); hold on;



