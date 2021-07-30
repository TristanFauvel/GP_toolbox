clear all
add_gp_module
figure_path = '/home/tfauvel/Documents/PhD/Figures/Thesis_figures/Chapter_1/';
close all
rng(1)
graphics_style_paper;
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

rng(4)
link = @normcdf; %inverse link function
regularization = 'nugget';
post = [];
%% Define the range of parameters
n = 50;
x = linspace(0,1, n);
d =1;
ntr =10;

[p,q]= meshgrid(x);
x2d = [p(:), q(:)]';


modeltype = 'exp_prop'; % Approximation method
x0 = 0;

base_kernelname = 'Matern52';
original_kernelfun =  @Matern52_kernelfun;%kernel used within the preference learning kernel, for subject = computer
theta= [-1;1];

base_kernelname = 'ARD_kernelfun';
original_kernelfun =  @ARD_kernelfun;%kernel used within the preference learning kernel, for subject = computer
theta= [3,0];

kernel_approximation = 'RRGP';

base_kernelfun = @(theta, xi, xj, training, reg) conditioned_kernelfun(theta, original_kernelfun, xi, xj, training, x0, reg);
kernelfun_cond= @(theta, xi, xj, training, reg) conditional_preference_kernelfun(theta, original_kernelfun, xi, xj, training, reg,x0);
kernelfun= @(theta, xi, xj, training, reg) preference_kernelfun(theta, original_kernelfun, xi, xj, training, reg);

link = @normcdf; %inverse link function for the classification model

% gfunc = @(x) forretal08(x)/10;
% gfunc = @(x) normpdf(x, 0.5, 0.2);
% g = gfunc(x)-gfunc(x0);

g = mvnrnd(zeros(1,n),base_kernelfun(theta, x, x, 'false', 'no'));
% g = g-g(1);

f = g-g';
f= f(:);

rd_idx = randsample(size(x2d,2), ntr, 'true');
xtrain= x2d(:,rd_idx);
ytrain= f(rd_idx);
ctrain = link(ytrain)>rand(ntr,1);


[mu_c,  mu_f, sigma2_f] = prediction_bin(theta, xtrain(:,1:ntr), ctrain(1:ntr), x2d, kernelfun, modeltype, post, regularization);
[mu_c_x0,  mu_g, sigma2_g, Sigma2_g] = prediction_bin(theta, xtrain(:,1:ntr), ctrain(1:ntr), [x; x0*ones(1,n^d)], kernelfun, modeltype, post, regularization);
    
[mu_c_cond,  mu_f_cond, sigma2_f_cond] = prediction_bin(theta, xtrain(:,1:ntr), ctrain(1:ntr), x2d, kernelfun_cond, modeltype, post, regularization);
[mu_c_cond_x0,  mu_g_cond, sigma2_g_cond, Sigma2_g_cond, dmuc_dx, dmuy_dx, dsigma2y_dx, dSigma2y_dx, var_muc] = prediction_bin(theta, xtrain(:,1:ntr), ctrain(1:ntr), [x; x0*ones(d,n^d)], kernelfun_cond, modeltype, post, regularization);


%% Find the true global optimum of g
[gmax, id_xmax] = max(g);
xmax = x(id_xmax);
% legend_pos = [-0.2,1];
legend_pos = [0.05,1];

mr = 2;
mc = 2;
i = 0;
fig=figure('units','centimeters','outerposition',1+[0 0 fwidth fheight(mr)]);
fig.Color =  [1 1 1];


tiledlayout(mr,mc, 'TileSpacing' , 'tight', 'Padding', 'tight')
nexttile
i=i+1;

plot_gp(x, mu_g, sigma2_g, C(2,:), linewidth); hold on;
plot(x, g,'color', 'k', 'linewidth', linewidth); hold off;

xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$g(x)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)

ylim = get(gca, 'Ylim');
%title('Inferred value function $g(x)$','Fontsize', Fontsize)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off
% pbaspect([1 1 1])

nexttile
i=i+1;

plot_gp(x, mu_g_cond, sigma2_g_cond, C(1,:), linewidth); hold on;
plot(x, g,'color', 'k', 'linewidth', linewidth); hold on;

xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$g(x)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
set(gca, 'Ylim', ylim, 'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
%title('Inferred value function $g(x)$','Fontsize', Fontsize)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off
% pbaspect([1 1 1])

nexttile
i=i+1;

plot(x, mu_g,'color', C(2,:), 'linewidth', linewidth); hold on;
plot(x, mu_g_cond, 'color', C(1,:), 'linewidth', linewidth); hold off;
xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$\mu_g(x)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
%title('Inferred value function $g(x)$','Fontsize', Fontsize)
legend('Without conditionning', 'Conditionning  on $g(x_0) = 0$')
legend box off
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off
% pbaspect([1 1 1])


nexttile
i=i+1;

plot(x, sigma2_g,'color', C(2,:), 'linewidth', linewidth); hold on;
plot(x, sigma2_g_cond, 'color', C(1,:), 'linewidth', linewidth); hold off;

xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$\sigma^2(x)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
%title('Inferred value function $g(x)$','Fontsize', Fontsize)
legend('Without conditionning', 'Conditionning  on $g(x_0) = 0$')
legend box off
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off
% pbaspect([1 1 1])



figname  = 'preference_learning_GP_comparison';
folder = [figure_path,figname];
savefig(fig, [folder,'/', figname, '_', num2str(ntr), '.fig'])
exportgraphics(fig, [folder,'/' , figname,  '_',num2str(ntr),'.pdf']);
exportgraphics(fig, [folder,'/' , figname,  '_',num2str(ntr),'.png'], 'Resolution', 300);




mr = 1;
mc = 2;
fig=figure('units','centimeters','outerposition',1+[0 0 fwidth fheight(mr)]);
fig.Color =  [1 1 1];

tiledlayout(mr,mc, 'TileSpacing' , 'tight', 'Padding', 'tight')
plot(x, mu_c_x0,'color', C(2,:), 'linewidth', linewidth); hold on;
plot(x, mu_c_cond_x0,'color', C(1,:), 'linewidth', linewidth); hold on;
plot(x, normcdf(g),'color', 'k', 'linewidth', linewidth); hold off;
xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$g(x)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)

ylim = get(gca, 'Ylim');
box off



mr = 2;
mc = 2;
fig=figure('units','centimeters','outerposition',1+[0 0 fwidth fheight(mr)]);
fig.Color =  [1 1 1];
tiledlayout(mr,mc, 'TileSpacing' , 'tight', 'Padding', 'tight')
nexttile()
ax1 = imagesc(x, x, reshape(mu_c, n,n)); hold on;
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

nexttile()
ax1 = imagesc(x, x, reshape(mu_c_cond, n,n)); hold on;
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

nexttile()
ax1 = imagesc(x, x, reshape(sigma2_f, n,n)); hold on;
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

nexttile()
ax1 = imagesc(x, x, reshape(sigma2_f_cond, n,n)); hold on;
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


