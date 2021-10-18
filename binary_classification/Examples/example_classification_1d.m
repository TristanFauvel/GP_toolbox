%% Gaussian process classification

close all;
add_gp_module;
rng(11)
graphics_style_paper;


% Define the likelihood approximation method
modeltype = 'laplace'; % 'exp_prop' or 'laplace'
post = [];
regularization = 'none';

% Choose a link function
link = @logistic;
link = @normcdf;

% Choose a kernel
kernelfun = @ARD_kernelfun; 

% Kernel hyperparameters
theta_true= [3;3];

% Generate a function by sampling from the GP prior
n=100;
x = linspace(0,1,n);
y = mvnrnd(constant_mean(x,0), kernelfun(theta_true, x,x, false, 'none')); % Latent function
p= link(y); % P(c = 1)

% Create data points
ntr =120; % Number of data points
i_tr= randsample(n,ntr,'true');
x_tr = x(:,i_tr);
y_tr = y(:, i_tr);
c_tr = p(i_tr)>rand(1,ntr);

% Define test points
x_test = x;
y_test = y;

%% GP classification with the correct hyperparameters
theta.cov =theta_true ; 
theta.mean = 0;
% Compute the predictive distribution
D = 1;
meanfun = 0;
type = 'classification'; 
 hyps.ncov_hyp =2; % number of hyperparameters for the covariance function
hyps.nmean_hyp =1; % number of hyperparameters for the mean function
   hyps.hyp_lb = -10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
hyps.hyp_ub = 10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
lb = 0;
ub = 1;

model = gp_classification_model(D, meanfun, kernelfun, regularization, hyps, lb,ub, type, link, modeltype);

[mu_c,  mu_y, sigma2_y, Sigma2_y, dmuc_dx, dmuy_dx, dsigma2y_dx,var_muc, dvar_muc_dx]= model.prediction(theta, x_tr, c_tr, x_test, post);


%% Plotting !
mr = 3;
mc = 2;
Xlim= [min(x),max(x)];
Ylim = [-10,10];

h=figure(1);
h.Color =  [1 1 1];
subplot(mr,mc,1)
b = plot_gp(x, mu_y, sigma2_y, C(1,:), linewidth); hold on
a = plot(x,y,'LineWidth',1.5,'Color', C(2,:)); hold off;
legend([a,b], 'True activation function', 'Inferred latent function')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',Ylim)
grid off
box off
title('True hyperparameters');


subplot(mr,mc,2)
scatter(x_tr, c_tr, 'MarkerFaceColor', C(1,:), 'MarkerEdgeColor', 'none') ; hold on;
plot(x_test, mu_c, 'LineWidth',1.5,'Color', C(2,:)) ; hold on;
plot(x_test, p, 'LineWidth',1.5,'Color',  C(1,:)) ; hold off;
legend('Data points', 'Inferred probability', 'True probability')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',[0,1])
grid off
box off
title('True hyperparameters');

% GP classification with the wrong hyperparameters
theta.cov = rand(size(theta_true));

[mu_c,  mu_y, sigma2_y]= model.prediction(theta, x_tr, c_tr, x_test, post);

subplot(mr,mc,3)
b = plot_gp(x, mu_y, sigma2_y, C(1,:), linewidth); hold on
a = plot(x,y,'LineWidth',1.5,'Color', C(2,:)); hold off;
legend([a,b], 'True activation function', 'Inferred activation function')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',Ylim)
grid off
box off
title('Wrong hyperparameters');

subplot(mr,mc,4)
scatter(x_tr, c_tr, 'MarkerFaceColor', C(1,:), 'MarkerEdgeColor', 'none') ; hold on;
plot(x_test, mu_c, 'LineWidth',1.5,'Color', C(2,:)) ; hold on;
plot(x_test, p, 'LineWidth',1.5,'Color',  C(1,:)) ; hold off;
legend('Data points', 'Inferred probability', 'True probability')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',[0,1])
grid off
box off
title('Wrong hyperparameters');


%% Local optimization of hyperparameters
theta = model.model_selection(x_tr, c_tr, theta, 'cov');
 
%% Prediction with the new hyperparameters
[mu_c,  mu_y, sigma2_y]= model.prediction(theta, x_tr, c_tr, x_test, post);

subplot(mr,mc,5)
b = plot_gp(x, mu_y, sigma2_y, C(1,:), linewidth); hold on
a = plot(x,y,'LineWidth',1.5,'Color', C(2,:)); hold off;
legend([a,b], 'True activation function', 'Inferred activation function')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',Ylim)
grid off
box off
title('Hyperparameters inferred using maximum likelihood');

subplot(mr,mc,6)

scatter(x_tr, c_tr, 'MarkerFaceColor', C(1,:), 'MarkerEdgeColor', 'none') ; hold on;
plot(x_test, mu_c, 'LineWidth',1.5,'Color', C(2,:)) ; hold on;
plot(x_test, p, 'LineWidth',1.5,'Color',  C(1,:)) ; hold off;
legend('Data points', 'Inferred probability', 'True probability')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',[0,1])
grid off
box off
title('Hyperparameters inferred using maximum likelihood');


