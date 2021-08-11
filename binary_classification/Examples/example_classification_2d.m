%% GP classification with 2D inputs

clear all;
close all;

add_gp_module;
rng(11)
graphics_style_paper;


n=30;
modeltype = 'exp_prop';
post = [];
regularization = 'nugget';

kernelfun = @ARD_kernelfun;
d= 2; % Dimension of the input space

theta_true = [rand(d,1);2];

% Sample a function corresponding to the true kernel 
x_range= linspace(0,10,n);
[p,q] = meshgrid(x_range);
x= [p(:),q(:)]';
y = mvnrnd(constant_mean(x,0), kernelfun(theta_truex,x, true, 'nugget')); %generate a function
y=y-mean(y);

link = @normcdf;
p= link(y);

% Sample data points
ntr =100;
itrain= randsample(n^d,ntr);
xtrain = x(:,itrain);
ytrain = y(:, itrain);
ctrain = rand(1,ntr)<link(ytrain);


% Define test points
xtest = x;
ytest = y;
ptest = p;

%% With true hyperparameters
theta = theta_true;

[mu_c,  mu_y, sigma2_y, Sigma2_y, dmuc_dx, dmuy_dx, dsigma2y_dx,var_muc, dvar_muc_dx]= prediction_bin(theta, xtrain, ctrain, xtest, model, post);
Xlim= [0,10];
Ylim = [-5,5];

mr = 4;
mc = 2;
i = 1;
h=figure(1);
h.Color =  [1 1 1];
subplot(mr, mc,i)
imagesc(x_range, x_range, reshape(y, n, n)); hold on;
colorbar
set(gca,'YDir','normal')
pbaspect([1,1,1])
title('True function')

i = i +1;
subplot(mr, mc,i) 
imagesc(x_range, x_range, reshape(p, n, n)); hold on;
colorbar
set(gca,'YDir','normal')
pbaspect([1,1,1])
title('True P(c=1)')


i = i +1;
subplot(mr, mc,i)
imagesc(x_range, x_range, reshape(mu_y, n, n)); hold off;
colorbar
set(gca,'YDir','normal')
pbaspect([1,1,1])
title('$\mu_y$, true hyperparameters')

i = i +1;
subplot(mr, mc,i)
imagesc(x_range, x_range, reshape(mu_c, n, n)); hold on;
colorbar()
set(gca,'YDir','normal')
scatter(xtrain(1,:), xtrain(2,:), 'MarkerFaceColor', colo(end,:), 'MarkerEdgeColor', 'none') ; hold on;
pbaspect([1,1,1]); 
title('$\mu_c$, true hyperparameters')


%% With wrong hyperparameters
theta = 0*rand(size(theta_true));

[mu_c,  mu_y, sigma2_y]= prediction_bin(theta, xtrain, ctrain, xtest, model, post);

i = i +1;
subplot(mr, mc, i)
imagesc(x_range, x_range, reshape(mu_y, n, n)); hold off;
set(gca,'YDir','normal')
pbaspect([1,1,1])
title('$\mu_y$, wrong hyperparameters')
colorbar

i = i +1;
subplot(mr, mc,i)
imagesc(x_range, x_range, reshape(mu_c, n, n)); hold on;
colorbar()
set(gca,'YDir','normal')
title('Inferred probability')
scatter(xtrain(1,:), xtrain(2,:), 'MarkerFaceColor', colo(end,:), 'MarkerEdgeColor', 'none') ; hold on;
pbaspect([1,1,1]); 
title('$\mu_c$, wrong hyperparameters')

%% Local optimization of hyperparameters
options=[];
theta = minFunc(@(hyp)negloglike_bin(hyp, xtrain, ctrain, kernelfun), theta, options);

% Prediction with the new hyperparameters
[mu_c,  mu_y, sigma2_y]= prediction_bin(theta, xtrain, ctrain, xtest, model, post);

i = i +1;
subplot(mr, mc,i)
imagesc(x_range, x_range, reshape(mu_y, n, n)); hold on;
set(gca,'YDir','normal')
title('Inferred probability')
pbaspect([1,1,1]); 
title('$\mu_y$, hyperparameters inferred with maximum likelihood')
colorbar

i = i +1;
subplot(mr, mc,i)
imagesc(x_range, x_range, reshape(mu_c, n, n)); hold on;
set(gca,'YDir','normal')
title('Inferred probability')
scatter(xtrain(1,:), xtrain(2,:), 'MarkerFaceColor', colo(end,:), 'MarkerEdgeColor', 'none') ; hold on;
pbaspect([1,1,1]); 
title('$\mu_c$, hyperparameters inferred with maximum likelihood')
colorbar



