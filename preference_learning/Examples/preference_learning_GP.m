clear all
add_gp_module
close all
rng(1)
graphics_style_paper;

post = [];
%% Define the range of parameters
n = 50;
lb = 0;
ub = 1;
x = linspace(lb,ub, n);
d =1;
ntr = 5;

[p,q]= meshgrid(x);
x2d = [p(:), q(:)]';


modeltype = 'exp_prop'; % Approximation method
kernelfun =  @Matern52_kernelfun;%kernel used within the preference learning kernel, for subject = computer
x0 = 0;

condition.x0 = x0;
condition.y0 = 0;

kernelname = 'Matern52';


approximationimation = 'RRGP';
 
link = @normcdf; %inverse link function for the classification model

theta.mean =0;
theta.cov= [-1;1];
g = mvnrnd(zeros(1,n),kernelfun(theta.cov, x, x, 'false', 'no'));
g = g-g(1);

f = g-g';
f= f(:);

nsamp= 500;
rd_idx = randsample(size(x2d,2), nsamp, 'true');
xtrain= x2d(:,rd_idx);
ytrain= f(rd_idx);
ctrain = link(ytrain)>rand(nsamp,1);

regularization = 'nugget';
meanfun = 0;
type = 'preference';
hyps.ncov_hyp =2; % number of hyperparameters for the covariance function
hyps.nmean_hyp =0; % number of hyperparameters for the mean function
hyps.hyp_lb = -10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
hyps.hyp_ub = 10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
D = 1;
 
model = gp_preference_model(D, meanfun, kernelfun, regularization, hyps, lb, ub, type, link, modeltype, kernelname,  []);
model_cond = gp_preference_model(D, meanfun, kernelfun, regularization, hyps, lb, ub, type, link, modeltype, kernelname,  condition);

[mu_c,  mu_f, sigma2_f] = model.prediction(theta, xtrain(:,1:ntr), ctrain(1:ntr), x2d, post);
[~,  mu_g, sigma2_g, Sigma2_g] = model.prediction(theta, xtrain(:,1:ntr), ctrain(1:ntr), [x; x0*ones(1,n^d)], post);
     
[mu_c_cond,  mu_f_cond, sigma2_f_cond] = model_cond.prediction(theta, xtrain(:,1:ntr), ctrain(1:ntr), x2d, post);
[mu_c_cond_x0,  mu_g_cond, sigma2_g_cond, Sigma2_g_cond, dmuc_dx, dmuy_dx, dsigma2y_dx, dSigma2y_dx, var_muc] = model_cond.prediction(theta, xtrain(:,1:ntr), ctrain(1:ntr), [x; x0*ones(d,n^d)], post);
 

%% Find the true global optimum of g
[gmax, id_xmax] = max(g);
xmax = x(id_xmax);
legend_pos = [-0.2,1];

mr = 1;
mc = 3;
i = 0;
fig=figure('units','centimeters','outerposition',1+[0 0 fwidth fheight(mr)]);
fig.Color =  background_color;


tiledlayout(mr,mc, 'TileSpacing' , 'tight', 'Padding', 'tight')
cl = [0,1];
nexttile
i=i+1;
imagesc(x, x, reshape(link(f),n,n), cl); hold on;
xlabel('$x$','Fontsize', Fontsize)
ylabel('$x''$','Fontsize', Fontsize)
title('$P(x>x'')$','Fontsize', Fontsize)
set(gca,'YDir','normal')
set(gca,'XTick',[0 0.5 1],'YTick',[0 0.5 1],'Fontsize', Fontsize)
pbaspect([1 1 1])
c = colorbar;
c.FontName = 'CMU Serif';
c.FontSize = Fontsize;
c.Limits = [0,1];
set(c, 'XTick', [0,1]);
colormap(cmap)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)

nexttile
i=i+1;
imagesc(x, x, reshape(mu_c_cond, n,n),cl); hold on;
p1 = scatter(xtrain(1, ctrain(1:ntr)==1),xtrain(2, ctrain(1:ntr)==1), markersize, 'o', 'k','filled'); hold on;
p2 = scatter(xtrain(1, ctrain(1:ntr)==0),xtrain(2, ctrain(1:ntr)==0), markersize, 'o','k'); hold off;
xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$x''$', 'Fontsize', Fontsize)
title('$\mu_c(x,x'')$','Fontsize', Fontsize)

set(gca,'YDir','normal')
set(gca,'XTick',[0 0.5 1],'YTick',[0 0.5 1], 'Fontsize', Fontsize)
pbaspect([1 1 1])
c = colorbar;
c.FontName = 'CMU Serif';
c.FontSize = Fontsize;
c.Limits = [0,1];
set(c, 'XTick', [0,1]);
colormap(cmap)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
legend([p1, p2], '$c=1$', '$c=0$','Fontsize',Fontsize, 'Location', 'northeast')
legend boxoff


nexttile
i=i+1;

h1 = plot_gp(x, mu_g_cond, sigma2_g_cond, C(1,:), linewidth); hold on;
h2 = plot(x, g, 'color', C(2,:), 'linewidth', linewidth);  hold off;

xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$g(x)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off
pbaspect([1 1 1])
legend([h1 h2], 'Posterior GP','True value function')
legend box off

