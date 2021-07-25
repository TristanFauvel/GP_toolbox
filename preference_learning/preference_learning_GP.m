clear all
add_gp_module
figure_path = '/home/tfauvel/Documents/PhD/Figures/Thesis_figures/Chapter_1/';
close all
rng(1)
graphics_style_paper;
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

link = @normcdf; %inverse link function

%% Define the range of parameters
n = 50;
x = linspace(0,1, n);
d =1;
ntr = 5;

[p,q]= meshgrid(x);
x2d = [p(:), q(:)]';


modeltype = 'exp_prop'; % Approximation method
original_kernelfun =  @Matern52_kernelfun;%kernel used within the preference learning kernel, for subject = computer
x0 = 0;

base_kernelname = 'Matern52';


kernel_approximation = 'RRGP';

base_kernelfun = @(theta, xi, xj, training, reg) conditioned_kernelfun(theta, original_kernelfun, xi, xj, training, x0, reg);
% kernelfun = @(theta, xi, xj, training) preference_kernelfun(theta, base_kernelfun, xi, xj, training);
kernelfun_cond= @(theta, xi, xj, training) conditional_preference_kernelfun(theta, original_kernelfun, xi, xj, training, x0);
kernelfun= @(theta, xi, xj, training) preference_kernelfun(theta, original_kernelfun, xi, xj, training);

link = @normcdf; %inverse link function for the classification model

% gfunc = @(x) forretal08(x)/10;
% gfunc = @(x) normpdf(x, 0.5, 0.2);
% g = gfunc(x)-gfunc(x0);

theta= [-1;1];
g = mvnrnd(zeros(1,n),base_kernelfun(theta, x, x, 'false', 'no'));
%g = g-g(1);

f = g'-g;
f= f(:);

nsamp= 500;
rd_idx = randsample(size(x2d,2), nsamp, 'true');
xtrain= x2d(:,rd_idx);
ytrain= f(rd_idx);
ctrain = link(ytrain)>rand(nsamp,1);


[mu_c,  mu_f, sigma2_f] = prediction_bin(theta, xtrain(:,1:ntr), ctrain(1:ntr), x2d, kernelfun, modeltype, post, regularization);
[~,  mu_g, sigma2_g, Sigma2_g] = prediction_bin(theta, xtrain(:,1:ntr), ctrain(1:ntr), [x; x0*ones(1,n^d)], kernelfun, modeltype, post, regularization);
mu_g = -mu_g; %(because prediction_bin considers P(x1 > x2);
    
[mu_c_cond,  mu_f_cond, sigma2_f_cond] = prediction_bin(theta, xtrain(:,1:ntr), ctrain(1:ntr), x2d, kernelfun_cond, modeltype, post, regularization);
[mu_c_cond_x0,  mu_g_cond, sigma2_g_cond, Sigma2_g_cond, dmuc_dx, dmuy_dx, dsigma2y_dx, dSigma2y_dx, var_muc] = prediction_bin(theta, xtrain(:,1:ntr), ctrain(1:ntr), [x; x0*ones(d,n^d)], kernelfun_cond, modeltype, post, regularization);
mu_g_cond = -mu_g_cond; %(because prediction_bin considers P(x1 > x2);


%% Find the true global optimum of g
[gmax, id_xmax] = max(g);
xmax = x(id_xmax);
legend_pos = [-0.2,1];

mr = 1;
mc = 3;
i = 0;
fig=figure('units','centimeters','outerposition',1+[0 0 width height(mr)]);
fig.Color =  [1 1 1];


tiledlayout(mr,mc, 'TileSpacing' , 'tight', 'Padding', 'tight')

% nexttile([1,2])
cl = [0,1];
nexttile
i=i+1;
imagesc(x, x, reshape(link(f),n,n), cl); hold on;
xlabel('$x$','Fontsize', Fontsize)
ylabel('$x''$','Fontsize', Fontsize)
title('$P(x''>x)$','Fontsize', Fontsize)
set(gca,'YDir','normal')
set(gca,'XTick',[0 0.5 1],'YTick',[0 0.5 1],'Fontsize', Fontsize)
pbaspect([1 1 1])
% c = colorbar;
% c.Limits = [0,1];
% set(c, 'XTick', [0,1]);
colormap(cmap)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)

% nexttile([1,2])
nexttile
i=i+1;
imagesc(x, x, reshape(mu_c_cond, n,n),cl); hold on;
scatter(xtrain(1, ctrain(1:ntr)==1),xtrain(2, ctrain(1:ntr)==1), markersize, 'o', 'k','filled'); hold on;
scatter(xtrain(1, ctrain(1:ntr)==0),xtrain(2, ctrain(1:ntr)==0), markersize, 'o','k'); hold off;
xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$x''$', 'Fontsize', Fontsize)
title('$\mu_c(x,x'')$','Fontsize', Fontsize)
% title('$P(x''>x | \mathcal{D})$')
set(gca,'YDir','normal')
set(gca,'XTick',[0 0.5 1],'YTick',[0 0.5 1], 'Fontsize', Fontsize)
pbaspect([1 1 1])
% c = colorbar;
% c.Limits = [0,1];
% set(c, 'XTick', [0,1]);
colormap(cmap)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)

% nexttile([1,2])
% nexttile
% i=i+1;
% imagesc(x, x, reshape(mu_f, n,n)); hold on;
% scatter(xtrain(1, ctrain(1:ntr)==1),xtrain(2, ctrain(1:ntr)==1), markersize, 'o', 'k','filled'); hold on;
% scatter(xtrain(1, ctrain(1:ntr)==0),xtrain(2, ctrain(1:ntr)==0), markersize, 'o','k'); hold off;
% xlabel('$x$')
% ylabel('$x''$')
% title('$\mu_f(x,x'')$')
% set(gca,'YDir','normal')
% set(gca,'XTick',[0 0.5 1],'YTick',[0 0.5 1])
% pbaspect([1 1 1])
% c = colorbar;
%text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
% colormap(cmap)

% nexttile([1,2])
nexttile
i=i+1;

plot_gp(x, mu_g_cond, sigma2_g_cond, C(1,:), linewidth); 

xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$g(x)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
%title('Inferred value function $g(x)$','Fontsize', Fontsize)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off
pbaspect([1 1 1])



figname  = 'preference_learning_GP';
folder = [figure_path,figname];
savefig(fig, [folder,'\', figname, '.fig'])
exportgraphics(fig, [folder,'\' , figname, '.pdf']);
exportgraphics(fig, [folder,'\' , figname, '.png'], 'Resolution', 300);

