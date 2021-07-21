%% I use this function to find the hyperparameters for the benchmarks
close all;
clear all
add_modules;
n=100;



rng(11)


objective = 'Ursem_waves';
kernelname = 'ARD';
kernelfun = @ARD_kernelfun;

[g, theta, lb, ub, lb_norm, ub_norm, theta_lb, theta_ub] = define_benchmark_PBO(objective, kernelname);
d = numel(lb);

if d==1
    x = linspace(lb(1), ub(2), n);
elseif d==2
    x_range= linspace(lb(1), ub(1),n);
    y_range= linspace(lb(2), ub(2),n);
    
    [p,q] = meshgrid(x_range, y_range);
    x= [p(:),q(:)]';
end

ntr =1000;

g = str2func(objective);

i_tr= randsample(n^d,ntr);
xtrain = x(:,i_tr);
ytrain = g(xtrain);
ytrain_norm = ytrain - mean(ytrain);

min_x = lb;
max_x = ub;
xtrain_norm = (xtrain - min_x)./(max_x-min_x);

x_norm = (x - min_x)./(max_x-min_x);

ncov_hyp = d+1;

hyp.cov = 10*ones(ncov_hyp,1);%rand(ncov_hyp,1);
hyp.mean = 0 ;


init_guess = [];
theta_lb = -15*ones(size(theta'));
theta_ub = 15*ones(size(theta'));
options_theta = [];
update = 'cov';
meanfun= @constant_mean;
nmean_hyp = 1;
theta = multistart_minConf(@(hyp)minimize_negloglike(hyp, xtrain_norm, ytrain_norm, kernelfun, meanfun, ncov_hyp, nmean_hyp, update), [theta_lb; 0], [theta_ub; 0],10, init_guess, options_theta);
hyp.cov = theta(1:end-1);



nsub = 100;
idsub = randsample(ntr, nsub);
% idsub = 1:ntr;
[mu_y, sigma2_y]= prediction(hyp, xtrainnorm(:,idsub), ytrain_norm(idsub), x_norm, kernelfun, meanfun);
y = g(x);
y = y - mean(y);

graphics_style_paper;
if d==1
    figure()
    subplot(1,3,1)
    plot(x, y)
    pbaspect([1 1 1])
    
    subplot(1,3,2)
    plot(x, mu_y)
    pbaspect([1 1 1])
    
    subplot(1,3,3)
    plot(x, sigma2_y)
    pbaspect([1 1 1])
    
elseif d ==2
    
    Fontsize = 14;
    set(0,'defaulttextInterpreter','latex')
    set(0,'defaultAxesTickLabelInterpreter','latex');
    set(0,'defaultLegendInterpreter','latex');
    set(0,'defaultAxesFontSize',Fontsize)
    set(0,'defaulttextFontSize',Fontsize)
    figure()
    subplot(1,3,1)
    imagesc(x_range, y_range, reshape(y, n, n)); hold on;
    scatter(xtrain(1,idsub), xtrain(2,idsub), 25, 'k','filled') ; hold off;
    pbaspect([1 1 1])
        colorbar()
    subplot(1,3,2)
    imagesc(x_range, y_range, reshape(mu_y, n, n)); hold on;
    pbaspect([1 1 1])
        colorbar()
    subplot(1,3,3)
    imagesc(x_range, y_range, reshape(sigma2_y, n, n)); hold on;
    pbaspect([1 1 1])
    colorbar()
end
