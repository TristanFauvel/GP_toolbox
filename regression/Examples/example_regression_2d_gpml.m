clear all;
close all;

rng(1)
addpath(genpath('/home/tfauvel/Desktop/optim_retina/GP_toolbox'))
n=50;
x1range = linspace(-10,10,n);
x2range = linspace(-10,10,n);
[p,q]= meshgrid(x1range,x2range);
x=[p(:),q(:)]';
y = sin(x(1,:)).*sin(x(2,:));

ntr =50;
i_tr= randsample(n*n,ntr);
x_tr = x(:,i_tr);
y_tr = y(:, i_tr);

x_test = x;
y_test = y;

[x_norm, min_x, max_x] = minmax_normalize(x_tr);
x_test_norm = (x_test - min_x)./(max_x-min_x);
mean_y = mean(y_tr);
y_norm = y_tr - mean_y;

nd=2;

% kernel_name = 'Gaussian_wnoise';
% % kernel_name = 'ARD';
% 
% if strcmp(kernel_name, 'ARD')
%     %ARD kernel
%     ncov_hyp=2+nd;
%     nmean_hyp=1;
%     kernelfun = @ARD_kernelfun_wnoise;
% elseif strcmp(kernel_name, 'Gaussian_wnoise')
%     %Gaussian kernel
%     ncov_hyp=3;
%     nmean_hyp=1;
%     kernelfun= @Gaussian_kernelfun_wnoise;
% end
% 

%% Regression

% Prior mean of the gaussian process

meanfunc = {@meanConst}; hyp.mean =0;
covfunc = {@covPoly ,2};
feval(covfunc{:})
hyp.cov = [0;4];% log([1 1 1 1 1 1 1 1 1 1 1])';


% covfunc = {@covPoly, 'ard', [1,0;0,2], 2};

% covPoly(mode,par,d,hyp,x,z)
 
likfunc = @likGauss;
sn = 0.1; hyp.lik = log(sn);
% m = hyp.cov; V = eye(2);
% pG = {@priorGaussMulti ,m,V};
% prior.cov = pG;
inf_func = {@infExact};


K = feval(covfunc{:}, hyp.cov, x');
mu = feval(meanfunc{:}, hyp.mean, x');
sample  = mvnrnd(mu,K);

figure()
surf(x1range, x2range, reshape(sample,n,n))


[a b mu_y sigma2_y]  = gp(hyp, inf_func, meanfunc, covfunc, likfunc, x_tr', y_tr', x_test');

figure()
subplot(3,1,1)
imagesc(x1range, x2range, reshape(y, n, n)); hold on;
scatter(x_tr(1,:), x_tr(2,:)) ; hold off;
subplot(3,1,2)
imagesc(x1range, x2range, reshape(mu_y, n, n)); hold on;
subplot(3,1,3)
imagesc(x1range, x2range, reshape(sqrt(sigma2_y), n, n)); hold on;





%% Prediction with the new hyperparameters
hyp= minimize(hyp, @gp, -400, inf_func, meanfunc, covfunc, likfunc, x_tr', y_tr');
[a b mu_y sigma2_y]  = gp(hyp, inf_func, meanfunc, covfunc, likfunc, x_tr', y_tr', x_test');

figure()
subplot(3,1,1)
imagesc(x1range, x2range, reshape(y, n, n)); hold on;
scatter(x_tr(1,:), x_tr(2,:)) ; hold off;
subplot(3,1,2)
imagesc(x1range, x2range, reshape(mu_y, n, n)); hold on;
subplot(3,1,3)
imagesc(x1range, x2range, reshape(sqrt(sigma2_y), n, n)); hold on;



