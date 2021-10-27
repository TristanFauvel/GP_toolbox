clear all;
close all;
addpath(genpath('/home/tristan/Desktop/GP_toolbox'))

n=1000;
rng(10)

meanfun= @constant_mean;

ncov_hyp=3;
nmean_hyp=1;
kernelfun = @Rational_Quadratic_kernelfun;
kernelname = 'rational_quadratic';
theta_true = [1;2;0];

ub = 5;
lb = 0;
x = linspace(lb, ub,n);

y = mvnrnd(meanfun(x,0), kernelfun(theta_true, x,x, true, 'nugget'));


D = 1;
regularization = 'nugget';
type = 'regression';
hyps.ncov_hyp =3; % number of hyperparameters for the covariance function
hyps.nmean_hyp =1; % number of hyperparameters for the mean function
hyps.hyp_lb = -10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
hyps.hyp_ub = 10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
model = gp_regression_model(D, meanfun, kernelfun, regularization, hyps, lb, ub, kernelname);



h=figure(1);
h.Name = 'True function';
h.Color =  [1 1 1];
plot(x, y,'LineWidth',1.5)
xlabel('x','Fontsize',14)
ylabel('f(x)','Fontsize',14)
grid off
box off

maxiter = 20;
idx= randsample(n,maxiter);

x_test = x;
y_test = y;
theta.cov = theta_true;
theta.mean = 0;

c= othercolor('GnBu7');

    
%% Plot the prior
mu_y= meanfun(x_test, theta.mean);
sigma2_y = diag(kernelfun(theta.cov, x_test,x_test));

Fontsize =14;
Xlim =[0,5];

h=figure(21);
h.Name = 'Bayesian optimisation before hyperparameters optimization';
h.Color =  [1 1 1];
plot(x,y,'LineWidth',1.5, 'Color', c(end,:)); hold on;
errorshaded(x, mu_y, sqrt(sigma2_y), 'Color',  c(22,:),'LineWidth', 1.5, 'DisplayName','Prediction', 'Fontsize', 14, ...   
    'Xlim', [0,5]); hold off
xlabel('x','Fontsize',Fontsize)
ylabel('f(x)','Fontsize',Fontsize)
xlim([0,5])
grid off
box off

    
for ntr= 5:maxiter
    i_tr = idx(1:ntr);
    xtrain = x(:,i_tr);
    ytrain = y(:, i_tr); %+ 0.1*randn(1,ntr);
    
    if ntr ==5        
        h=figure(22);
        h.Name = 'Bayesian optimisation before hyperparameters optimization';
        h.Color =  [1 1 1];
        scatter(xtrain, ytrain, 'MarkerFaceColor', c(end,:))
        xlabel('x','Fontsize',Fontsize)
        ylabel('f(x)','Fontsize',Fontsize)
        xlim(Xlim)
        ylim([-4,3])
        set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim)    
        grid off
        box off
    end
    
    %     theta.cov = rand(ncov_hyp,1);
    %     theta.mean = zeros(nmean_hyp,1);
   
    [mu_y, sigma2_y]= model.prediction(theta, xtrain, ytrain, x_test, []);
    
    h=figure(2);
    h.Name = 'Bayesian optimisation before hyperparameters optimization';
    h.Color =  [1 1 1];
    plot(x,y,'LineWidth',1.5); hold on;
    errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', c(22,:),'DisplayName','Prediction'); hold on
    scatter(xtrain, ytrain, 'MarkerFaceColor', c(end,:), 'MarkerEdgeColor', 'none'); hold off;
    xlabel('x','Fontsize',Fontsize)
    ylabel('f(x)','Fontsize',Fontsize)
    set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim)    
    grid off
    box off
    
   theta = model.model_selection(xtrain, ytrain, theta, 'cov');
    
    %% Prediction with the new hyperparameters
    [mu_y, sigma2_y]= model.prediction(theta, xtrain, ytrain, x_test, []);
    
    K= kernelfun(theta_true, x(1:10:end),x(1:10:end));
    
    xplot =  x(1:10:end);
    [~,~,~, ~, K] =  model.prediction(theta, xtrain, ytrain, xplot, []);
    
    h=figure(2);
    h.Name = 'Bayesian optimisation after hyperparameters optimization';
    h.Color =  [1 1 1];
    plot(x,y,'LineWidth',1.5); hold on;
    scatter(xtrain, ytrain) ; hold on;
    errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
    xlabel('x','Fontsize',14)
    ylabel('f(x)','Fontsize',14)
    ylim([-5 5])
    box off     
       
    h=figure(3);
    set(h, 'Position', get(0, 'Screensize'));
    h.Name = 'Bayesian optimisation after hyperparameters optimization';
    h.Color =  [1 1 1];
    
    f2 =subplot(1,2,2);
    imagesc(K)
    colormap(c)
    grid off
    set(gca,'visible','off')
    box off
    pbaspect([1 1 1])
    
    f1=subplot(1,2,1);
    h.Name = 'Bayesian optimisation after hyperparameters optimization';
    h.Color =  [1 1 1];
    plot(x,y,'LineWidth',1.5); hold on;
    scatter(xtrain, ytrain) ; hold on;
    errorshaded(x, mu_y, sqrt(sigma2_y), 'Color', 'red','DisplayName','Prediction'); hold off
    ylim([-5 5])
    box off 
    set(gca,'XColor', 'none','YColor','none')
    

end
