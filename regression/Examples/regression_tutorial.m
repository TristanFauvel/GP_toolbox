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
x_tr = x(:,i_tr);
y_tr = y(:, i_tr);

sample_prior = mvnrnd(mu_y_ard, Kard);

model.kernelfun = kernelfun;
model.meanfun = meanfun;
model.regularization = regularization;
[mu_y, sigma2_y,dmu_dx, sigma2_dx, Sigma2_y, dSigma2_dx, post] = prediction(theta, x_tr, y_tr, x, model, []); 
sample_post = mvnrnd(mu_y, Sigma2_y);

theta_w = theta;
theta_w.cov = [4,3];
[mu_y_w, sigma2_y_w,~,~, Sigma2_y_w] = prediction(theta_w, x_tr, y_tr, x, model, []);
sample_post_w = mvnrnd(mu_y, Sigma2_y);

mr = 1;
mc = 3;
legend_pos = [-0.18,1];

fig=figure('units','centimeters','outerposition',1+[0 0 fwidth fheight(1)]);
fig.Color =  [1 1 1];
layout = tiledlayout(mr,mc, 'TileSpacing', 'tight', 'padding','tight');
i = 0;

nexttile();
i=i+1;
plot_gp(x,mu_y_ard, sqrt(diag(Kard)), C(1,:),linewidth);
plot(x, y, 'Color',  C(2,:),'LineWidth', linewidth); hold on;
plot(x, sample_prior, 'Color',  'k','LineWidth', linewidth/2); hold off;

box off
xlabel('$x$')
ylabel('$f(x)$')
yl = get(gca,'Ylim');
set(gca, 'Xlim', [0,1], 'Xtick', [0,0.5,1], 'Ytick', floor([yl(1), 0, yl(2)]), 'Fontsize', Fontsize);
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)


nexttile();
i=i+1;
plot_gp(x,mu_y, sqrt(sigma2_y), C(1,:),linewidth);
plot(x, y, 'Color',  C(2,:),'LineWidth', linewidth); hold on;
% errorshaded(x,mu_y, sqrt(sigma2_y), 'Color',  C(1,:),'LineWidth', linewidth, 'Fontsize', Fontsize); hold on
plot(x_tr, y_tr, 'ro', 'MarkerSize', 10, 'color', C(2,:)); hold on;
scatter(x_tr, y_tr, 2*markersize, C(2,:), 'filled'); hold on;
plot(x, sample_post, 'Color',  'k','LineWidth', linewidth/2); hold off;
ylabel('$f(x)$')
box off
%title('Posterior distribution','Fontsize',Fontsize, 'interpreter', 'latex')
xlabel('$x$')
set(gca, 'Xlim', [0,1], 'Xtick', [0,0.5,1], 'Ylim', yl, 'Ytick', floor([yl(1), 0, yl(2)]), 'Fontsize', Fontsize');
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)

nexttile();
i=i+1;
ax1 = imagesc(x, x, Sigma2_y); hold on;
xlabel('$x$','Fontsize',Fontsize)
ylabel('$x''$','Fontsize',Fontsize)
title('Cov$(f(x),f(x'')|\mathcal{D})$','Fontsize',Fontsize, 'interpreter', 'latex')
set(gca,'YDir','normal')
pbaspect([1 1 1])
set(gca, 'Xlim', [0,1], 'Xtick', [0,0.5,1],'Ylim', [0,1], 'Ytick', [0,0.5,1], 'Fontsize', Fontsize');
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
cb = colorbar;
cb_lim = get(cb, 'Ylim');
cb_tick = get(cb, 'Ytick');
set(cb,'Ylim', cb_lim, 'Ytick', [cb_tick(1),cb_tick(end)]);

colormap(cmap)
box off

%% Now, I plot the prior distribution corresponding to each kernel, along with samples from this distribution


figname  = 'GP_regression';
folder = ['/home/tfauvel/Documents/PhD/Figures/Thesis_figures/Chapter_1/',figname];
savefig(fig, [folder,'/', figname, '.fig'])
exportgraphics(fig, [folder,'/' , figname, '.pdf']);
exportgraphics(fig, [folder,'/' , figname, '.png'], 'Resolution', 300);



