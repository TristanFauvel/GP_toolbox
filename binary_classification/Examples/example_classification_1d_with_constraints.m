%% Gaussian process classification with constraints (g(x0) = y0)

clear all;
close all;
add_modules;
n=1000;
graphics_style_paper
rng(11)
modeltype = 'exp_prop';

link = @normcdf;

kernelfun = @ARD_kernelfun;
theta_true= [1;1];

x = linspace(0,1,n);
y = mvnrnd(constant_mean(x,0), kernelfun(theta_true, x,x)); %generate a function
y=y-y(1);

p= link(y);
c = rand(1,n)<p;
ntr =30;

i_tr= randsample(n,ntr);
x_tr = x(:,i_tr);
y_tr = y(:, i_tr);
c_tr = c(:, i_tr);

x_test = x;
y_test = y;
c_test = c;

theta =theta_true ; % rand(size(theta_true));

x0 = x(1); %constraint
y0 = y(1);
[mu_c,  mu_y, sigma2_y, Sigma2_y, ~, dmuy_dx, dsigma2y_dx,var_muc]= prediction_bin_with_constraints(theta, x_tr, c_tr, x_test, kernelfun, x0, y0, 'modeltype', 'laplace');


h = figure();
h.Color =  [1 1 1];
h.Name = 'GP classification with true hyperparameters';
subplot(1,3,1)
scatter(x_tr, c_tr, 'MarkerFaceColor', colo(end,:), 'MarkerEdgeColor', 'none') ; hold on;
plot(x_test, mu_c, 'LineWidth',1.5,'Color', colo(22,:)) ; hold on;
plot(x_test, p, 'LineWidth',1.5,'Color',  colo(end,:)) ; hold off;
legend('Data points', 'Inferred probability', 'True probability')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',[0,1])
grid off
box off
subplot(1,3,2)
plot(x_test, var_muc, 'LineWidth',1.5,'Color',  colo(end,:)) ; hold off;
legend('Data points', 'Inferred probability', 'True probability')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',[0,1])
grid off
box off
subplot(1,3,3)
Xlim= [min(x),max(x)];
Ylim = [-5,5];
plot(x,y,'LineWidth',1.5,'Color', colo(end,:)); hold on;
errorshaded(x, mu_y,sigma2_y, 'LineWidth',1.5,'Color', colo(22,:)); hold off
legend('True activation function', 'Inferred activation function')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',Ylim)
grid off
box off




%%
theta = rand(size(theta_true));

[mu_c,  mu_y, sigma2_y]= prediction_bin(theta, x_tr, c_tr, x_test, kernelfun, 'modeltype', modeltype);
Fontsize = 14;
set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex');  
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',Fontsize)
set(0,'defaulttextFontSize',Fontsize)

Ylim = [-5,5];
h=figure(3);
h.Color =  [1 1 1];
h.Name = 'Activation function';
plot(x,y,'LineWidth',1.5,'Color', colo(end,:)); hold on;
plot(x, mu_y,'LineWidth',1.5,'Color', colo(22,:)); hold off
legend('True activation function', 'Inferred activation function')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',Ylim)
grid off
box off

h=figure(4);
h.Color =  [1 1 1];
h.Name = 'GP classification before hyperparameters optimization';
scatter(x_tr, c_tr, 'MarkerFaceColor', colo(end,:), 'MarkerEdgeColor', 'none') ; hold on;
plot(x_test, mu_c, 'LineWidth',1.5,'Color', colo(22,:)) ; hold on;
plot(x_test, p, 'LineWidth',1.5,'Color',  colo(end,:)) ; hold off;
legend('Data points', 'Inferred probability', 'True probability')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',[0,1])
grid off
box off

%% Local optimization of hyperparameters
options=[];
theta = minFunc(@(hyp)negloglike_bin(hyp, x_tr, c_tr, kernelfun), theta, options);

%% Prediction with the new hyperparameters
[mu_c,  mu_y, sigma2_y]= prediction_bin(theta, x_tr, c_tr, x_test, kernelfun);

h=figure(5);
h.Color =  [1 1 1];
h.Name = 'GP classification after hyperparameters optimization';
scatter(x_tr, c_tr, 'MarkerFaceColor', colo(end,:), 'MarkerEdgeColor', 'none') ; hold on;
plot(x_test, mu_c, 'LineWidth',1.5,'Color', colo(22,:)) ; hold on;
plot(x_test, p, 'LineWidth',1.5,'Color',  colo(end,:)) ; hold off;
legend('Data points', 'Inferred probability', 'True probability')
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',[0,1])
grid off
box off


