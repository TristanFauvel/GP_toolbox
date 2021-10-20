clear all
add_gp_module
figure_path = '/home/tfauvel/Documents/PhD/Figures/Thesis_figures/Chapter_1/';

graphics_style_paper
close all
rng(1)
regularization = 'nugget';
post = [];
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

link = @normcdf; %inverse link function

%% Define the range of parameters
n = 50;
lb = 0;
ub = 1;

x = linspace(lb, ub, n);
d =1;
ntr = 5;

[p,q]= meshgrid(x);
x2d = [p(:), q(:)]';

x0 = x(:,1);

modeltype = 'exp_prop'; % Approximation method
base_kernelfun =  @Matern52_kernelfun;%kernel used within the preference learning kernel, for subject = computer
base_kernelname = 'Matern52';
approximationimation = 'RRGP';
kernelfun = @(theta, xi, xj, training, reg) preference_kernelfun(theta, base_kernelfun, xi, xj, training, reg);
link = @normcdf; %inverse link function for the classification model

% gfunc = @(x) forretal08(x)/10;
% gfunc = @(x) normpdf(x, 0.5, 0.2);
% g = gfunc(x)-gfunc(x0);

theta.cov= [-1;1];
g = mvnrnd(zeros(1,n),base_kernelfun(theta.cov, x, x, 'false', regularization));
g = g-g(1);

f = g'-g;
f= f(:);

nsamp= 500;
rd_idx = randsample(size(x2d,2), nsamp, 'true');
xtrain= x2d(:,rd_idx);
ytrain= f(rd_idx);
ctrain = link(ytrain)>rand(nsamp,1);



%%
regularization = 'nugget';
model.kernelfun = kernelfun;
model.link = link;
model.modeltype = modeltype;
meanfun = 0;
type = 'preference';
hyps.ncov_hyp =2; % number of hyperparameters for the covariance function
hyps.nmean_hyp =0; % number of hyperparameters for the mean function
hyps.hyp_lb = -10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
hyps.hyp_ub = 10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
D = 1;
 
model = gp_classification_model(D, meanfun, kernelfun, regularization, hyps, lb, ub, type, link, modeltype);


%%
[mu_c,  mu_f, sigma2_f] = model.prediction(theta, xtrain(:,1:ntr), ctrain(1:ntr), x2d, post);

[~,  mu_g, sigma2_g, Sigma2_g] = model.prediction(theta, xtrain(:,1:ntr), ctrain(1:ntr), [x; x0*ones(1,n^d)], post);
mu_g = -mu_g; %(because prediction_bin considers P(x1 > x2);
    


%% Find the true global optimum of g
[gmax, id_xmax] = max(g);
xmax = x(id_xmax);
legend_pos = [-0.18,1.15];

fig=figure('units','centimeters','outerposition',1+[0 0 16 1.2/2*16]);
fig.Color =  background_color;

mr = 1;
mc = 3;
i = 0;
tiledlayout(mr,mc, 'TileSpacing' , 'tight', 'Padding', 'tight')

% nexttile([1,2])
cl = [0,1];
nexttile
i=i+1;
imagesc(x, x, reshape(link(f),n,n), cl); hold on;
xlabel('$x$')
ylabel('$x''$')
title('$P(x''>x)$')
set(gca,'YDir','normal')
set(gca,'XTick',[0 0.5 1],'YTick',[0 0.5 1])
pbaspect([1 1 1])
% c = colorbar;
% c.Limits = [0,1];
% set(c, 'XTick', [0,1]);
colormap(cmap)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)

% nexttile([1,2])
nexttile
i=i+1;
imagesc(x, x, reshape(mu_c, n,n),cl); hold on;
scatter(xtrain(1, ctrain(1:ntr)==1),xtrain(2, ctrain(1:ntr)==1), markersize, 'o', 'k','filled'); hold on;
scatter(xtrain(1, ctrain(1:ntr)==0),xtrain(2, ctrain(1:ntr)==0), markersize, 'o','k'); hold off;
xlabel('$x$')
ylabel('$x''$')
title('$\mu_c(x,x'')$')
% title('$P(x''>x | \mathcal{D})$')
set(gca,'YDir','normal')
set(gca,'XTick',[0 0.5 1],'YTick',[0 0.5 1])
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
options.handle = fig;
options.alpha = 0.2;
options.line_width = linewidth;
options.color_area = cmap(1,:);%[23, 20, 196]./255;    % Blue theme
options.color_line = cmap(1,:);%[23, 20, 196]./255;
options.x_axis = x;
options.line_width = linewidth;
h1 = plot_area(mu_g', sqrt(sigma2_g'), options); hold on;
plot(x,g, 'Color',  cmap(end,:),'linewidth',linewidth); hold on
xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$g(x)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1])
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3))
%title('Inferred value function $g(x)$','Fontsize', Fontsize)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off
pbaspect([1 1 1])


