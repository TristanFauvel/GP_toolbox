clear all;
close all;


%%
j = 3;
objectives = {'levy', 'goldpr', 'camel6'};
objective = objectives{j};
if strcmp(objective, 'levy')  
%     4.9244
%     5.4415
%     6.8752
    d=2; %dimension of the input space
    xbounds = [-10, 10; -10, 10];
elseif strcmp(objective, 'goldpr')
    d=2; %dimension of the input space
    xbounds = [-2, 2; -2, 2];
%      4.3574
%     6.5804
%    10.0000

elseif strcmp(objective, 'camel6')
    d=2; %dimension of the input space
    xbounds = [-3, 3; -2, 2];
%     4.0870
%     2.3590
%     8.6282
end

objective = str2func(objective);
lb = xbounds(:,1);
ub =xbounds(:,2);
%%
rng(1)
addpath(genpath('/home/tfauvel/Desktop/optim_retina/GP_toolbox'))
n=100;
x1range = linspace(lb(1),ub(1),n);
x2range = linspace(lb(2),ub(2),n);
[p,q]= meshgrid(x1range,x2range);
x=[p(:),q(:)]';
y = objective(x); % sin(x(1,:)).*sin(x(2,:));

ntr =400;
i_tr= randsample(n*n,ntr);
x_tr = x(:,i_tr);
y_tr = y(:, i_tr);

x_test = x;
y_test = y;

x_norm = (x_tr model.lb)./(model.ub - model.lb);
x_test_norm = (x_test model.lb)./(model.ub - model.lb);
mean_y = mean(y_tr);
y_norm = y_tr - mean_y;

meanfun= @constant_mean;


kernel_name = 'ARD';

if strcmp(kernel_name, 'ARD')
    %ARD kernel
    ncov_hyp=1+d;
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

[mu_y, sigma2_y]= prediction(theta, x_norm, y_norm, x_test_norm, kernelfun, meanfun);
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

%% Global optimization of hyperparameters
gs = GlobalSearch;
gs.Display='iter';
lb= -10*ones(1,numel(hyp));
ub= 10*ones(1,numel(hyp));
opts = optimoptions(@fmincon,'Algorithm','sqp');
myobjfun = @(hyp) minimize_negloglike(hyp, x_norm, y_norm, kernelfun, meanfun, ncov_hyp, nmean_hyp, update);
problem = createOptimProblem('fmincon','objective', myobjfun,'x0',hyp, 'lb',lb,'ub',ub, 'options', opts);
ms = MultiStart;
hyp = run(ms,problem,200);
theta.cov = hyp(1:ncov_hyp);
theta.mean = hyp(ncov_hyp+1:ncov_hyp+nmean_hyp);

% %% Local optimization of hyperparameters
% options=[];
% hyp = minFunc(@(hyp)minimize_negloglike(hyp, x_norm, y_norm, kernelfun, meanfun, ncov_hyp, nmean_hyp, update), hyp, options);
% theta.cov = hyp(1:ncov_hyp);
% theta.mean = hyp(ncov_hyp+1:ncov_hyp+nmean_hyp);

%% Prediction with the new hyperparameters
[mu_y, sigma2_y]= prediction(theta, x_norm, y_norm, x_test_norm, kernelfun, meanfun);
mu_y = mu_y + mean_y;

figure()
subplot(3,1,1)
imagesc(x1range, x2range, reshape(y, n, n)); hold on;
scatter(x_tr(1,:), x_tr(2,:)) ; hold off;
subplot(3,1,2)
imagesc(x1range, x2range, reshape(mu_y, n, n)); hold on;
subplot(3,1,3)
imagesc(x1range, x2range, reshape(sqrt(sigma2_y), n, n)); hold on;



