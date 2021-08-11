clear all
close all;
rng(3)
n = 100;
x = linspace(0,1,n);
graphics_style_paper;
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

% Prior mean of the gaussian process
regularization = 'nugget';
meanfun= @constant_mean;

kernelfun =@ARD_kernelfun;
theta.cov = [3,2];
theta.mean = 0;
Kard = kernelfun(theta.cov,x,x, 'true', regularization);
mu_y_ard = meanfun(x,theta.mean);
y = mvnrnd(mu_y_ard, Kard);

ntr = 2; 
i_tr= randsample(n,ntr);
% i_tr(3)=100 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ntr = 3; %%%%%%%%%%%
xtrain = x(:,i_tr);
ytrain = y(:, i_tr);


[mu_y, sigma2_y,dmu_dx, sigma2_dx, Sigma2_y, dSigma2_dx, post] = prediction(theta, xtrain, ytrain, x, kernelfun, meanfun, [], regularization); 
sample_post = mvnrnd(mu_y, Sigma2_y);

theta_w = theta;
theta_w.cov = [4,3];
[mu_y_w, sigma2_y_w,~,~, Sigma2_y_w] = prediction(theta_w, xtrain, ytrain, x, kernelfun, meanfun, [], regularization);
sample_post_w = mvnrnd(mu_y_w, Sigma2_y_w);


hyp.mean = 0;
hyp.cov = theta.cov;
update = 'cov';
ncov_hyp = 2;
nmean_hyp = 1;
update = 'cov';
init_guess = hyp;
theta_lb = -10*ones(1,2);
theta_ub = 10*ones(1,2);
options_theta.method = 'lbfgs';
hyp.cov = multistart_minConf(@(hyp)minimize_negloglike(hyp, xtrain, ytrain, kernelfun, meanfun, ncov_hyp, nmean_hyp, update), [theta_lb,0], [theta_ub,0],10, [], options_theta);
hyp.cov = hyp.cov(1:2);
[mu_y_ml, sigma2_y_ml,~,~, Sigma2_y_ml] = prediction(hyp, xtrain, ytrain, x, kernelfun, meanfun, [], regularization);
sample_post_ml = mvnrnd(mu_y_ml, Sigma2_y_ml);



mr = 1;
mc = 3;
legend_pos = [-0.18,1];

fig=figure('units','centimeters','outerposition',1+[0 0 fwidth fheight(1)]);
fig.Color =  [1 1 1];
layout = tiledlayout(mr,mc, 'TileSpacing', 'tight', 'padding','tight');
i = 0;


nexttile();
i=i+1;
plot_gp(x,mu_y_w, sqrt(sigma2_y_w), C(1,:),linewidth);
plot(x, y, 'Color',  C(2,:),'LineWidth', linewidth); hold on;
% errorshaded(x,mu_y_w, sqrt(sigma2_y_w), 'Color',  C(1,:),'LineWidth', linewidth, 'Fontsize', Fontsize); hold on
%errorshaded(x,mu_y_w, sqrt(sigma2_y_w), 'Color',  C(1,:),'LineWidth', linewidth, 'Fontsize', Fontsize); hold on
plot(xtrain, ytrain, 'ro', 'MarkerSize', 10, 'color', C(2,:)); hold on;
scatter(xtrain, ytrain, 2*markersize, C(2,:), 'filled'); hold on;
%title('Posterior distribution with wrong hyperparameters','Fontsize',Fontsize, 'interpreter', 'latex')
box off
xlabel('$x$')
ylabel('$f(x)$')
yl = get(gca,'Ylim');
set(gca, 'Xlim', [0,1], 'Xtick', [0,0.5,1], 'Ytick', floor([yl(1), 0, yl(2)]), 'Fontsize', Fontsize);
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)

nexttile();
i=i+1;
plot_gp(x,mu_y_ml, sqrt(sigma2_y_ml), C(1,:),linewidth);
plot(x, y, 'Color',  C(2,:),'LineWidth', linewidth); hold on;
% errorshaded(x,mu_y, sqrt(sigma2_y), 'Color',  C(1,:),'LineWidth', linewidth, 'Fontsize', Fontsize); hold on
plot(xtrain, ytrain, 'ro', 'MarkerSize', 10, 'color', C(2,:)); hold on;
scatter(xtrain, ytrain, 2*markersize, C(2,:), 'filled'); hold on;
ylabel('$f(x)$')
box off
%title('Posterior distribution','Fontsize',Fontsize, 'interpreter', 'latex')
xlabel('$x$')
set(gca, 'Xlim', [0,1], 'Xtick', [0,0.5,1], 'Ylim', yl, 'Ytick', floor([yl(1), 0, yl(2)]), 'Fontsize', Fontsize');
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)


nexttile();
i=i+1;
plot_gp(x,mu_y, sqrt(sigma2_y), C(1,:),linewidth);
plot(x, y, 'Color',  C(2,:),'LineWidth', linewidth); hold on;
% errorshaded(x,mu_y, sqrt(sigma2_y), 'Color',  C(1,:),'LineWidth', linewidth, 'Fontsize', Fontsize); hold on
plot(xtrain, ytrain, 'ro', 'MarkerSize', 10, 'color', C(2,:)); hold on;
scatter(xtrain, ytrain, 2*markersize, C(2,:), 'filled'); hold on;
ylabel('$f(x)$')
box off
%title('Posterior distribution','Fontsize',Fontsize, 'interpreter', 'latex')
xlabel('$x$')
set(gca, 'Xlim', [0,1], 'Xtick', [0,0.5,1], 'Ylim', yl, 'Ytick', floor([yl(1), 0, yl(2)]), 'Fontsize', Fontsize');
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)


colormap(cmap)
box off

%% Now, I plot the prior distribution corresponding to each kernel, along with samples from this distribution


figname  = 'GP_regression_hyps';
folder = ['/home/tfauvel/Documents/PhD/Figures/Thesis_figures/Chapter_1/',figname];
savefig(fig, [folder,'/', figname, '.fig'])
exportgraphics(fig, [folder,'/' , figname, '.pdf']);
exportgraphics(fig, [folder,'/' , figname, '.png'], 'Resolution', 300);


