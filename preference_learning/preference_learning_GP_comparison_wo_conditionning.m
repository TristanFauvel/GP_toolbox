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
%% Define the range of parameters
n = 50;
x = linspace(0,1, n);
d =1;
ntr =5;

[p,q]= meshgrid(x);
x2d = [p(:), q(:)]';


modeltype = 'exp_prop'; % Approximation method
x0 = 0.5;

base_kernelname = 'Matern52';
original_kernelfun =  @Matern52_kernelfun;%kernel used within the preference learning kernel, for subject = computer
theta= [-1;1];

base_kernelname = 'ARD_kernelfun';
original_kernelfun =  @ARD_kernelfun;%kernel used within the preference learning kernel, for subject = computer
theta= [3,0];

approximationimation = 'RRGP';

base_kernelfun = @(theta, xi, xj, training, reg) conditioned_kernelfun(theta, original_kernelfun, xi, xj, training, x0, reg);
kernelfun_cond= @(theta, xi, xj, training, reg) conditional_preference_kernelfun(theta, original_kernelfun, xi, xj, training, reg,x0);
kernelfun= @(theta, xi, xj, training, reg) preference_kernelfun(theta, original_kernelfun, xi, xj, training, reg);

 

model.link = link;
model.modeltype = modeltype;
model.kernelfun = kernelfun;
model.regularization = regularization;

model_cond = model;
model_cond.kernelfun = kernelfun_cond;

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


[mu_c,  mu_f, sigma2_f,~, ~, ~, ~, ~, var_muc_f] = prediction_bin(theta, xtrain(:,1:ntr), ctrain(1:ntr), x2d, model, []);
[mu_c_x0,  mu_g, sigma2_g, Sigma2_g, ~, ~, ~, ~, var_muc] = prediction_bin(theta, xtrain(:,1:ntr), ctrain(1:ntr), [x; x0*ones(1,n^d)], model, []);
    
[mu_c_cond,  mu_f_cond, sigma2_f_cond,~, ~, ~, ~, ~, var_muc_fcond] = prediction_bin(theta, xtrain(:,1:ntr), ctrain(1:ntr), x2d, model_cond, []);
[mu_c_cond_x0,  mu_g_cond, sigma2_g_cond, Sigma2_g_cond, ~, ~, ~, ~, var_muc_cond] = prediction_bin(theta, ...
    xtrain(:,1:ntr), ctrain(1:ntr), [x; x0*ones(d,n^d)], model_cond, []);


%% Find the true global optimum of g
[gmax, id_xmax] = max(g);
xmax = x(id_xmax);
% legend_pos = [-0.2,1];
legend_pos = [0.05,1];

mr = 3;
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
legend('Without conditionning', ['Conditionning  on $g(x_0) = ', num2str(x0), '$'])
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
legend('Without conditionning', ['Conditionning  on $g(x_0) = ', num2str(x0), '$'])
legend box off
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off
% pbaspect([1 1 1])

nexttile
i=i+1;

plot(x, var_muc,'color', C(2,:), 'linewidth', linewidth); hold on;
plot(x, var_muc_cond, 'color', C(1,:), 'linewidth', linewidth); hold off;

xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$V(\Phi(f(x)))$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
%title('Inferred value function $g(x)$','Fontsize', Fontsize)
legend('Without conditionning', ['Conditionning  on $g(x_0) = ', num2str(x0), '$'])
legend box off
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off


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



mr = 3;
mc = 2;
fig=figure('units','centimeters','outerposition',1+[0 0 fwidth fheight(mr)]);
fig.Color =  [1 1 1];
tiledlayout(mr,mc, 'TileSpacing' , 'tight', 'Padding', 'tight')
nexttile()
ax1 = imagesc(x, x, reshape(mu_c, n,n)); hold on;
xlabel('$x$','Fontsize',Fontsize)
ylabel('$x''$','Fontsize',Fontsize)
title('$\mu_c(x)$','Fontsize',Fontsize, 'interpreter', 'latex')
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
title('$\mu_c(x)$','Fontsize',Fontsize, 'interpreter', 'latex')
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



nexttile()
ax1 = imagesc(x, x, reshape(var_muc_f, n,n)); hold on;
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
ax1 = imagesc(x, x, reshape(var_muc_fcond, n,n)); hold on;
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

%%

 

acquisition_funs = {'active_sampling_binary'};
acquisition_fun = @active_sampling_binary;
acquisition_name = 'BALD';
maxiter = 80;%total number of iterations

nreplicates = 40; %20;

nacq = numel(acquisition_funs);


% wbar = waitbar(0,'Computing...');
rescaling = 0;
if rescaling ==0
    load('benchmarks_table.mat')
else
    load('benchmarks_table_rescaled.mat')
end
objectives = benchmarks_table.fName; %; 'Ursem_waves';'forretal08'; 'camel6';'goldpr'; 'grlee12';'forretal08'};
nobj =numel(objectives);
seeds = 1:nreplicates;
update_period = maxiter+2;
conditions = [0,1];

cum_regret_maxvar= NaN(nreps,maxiter+1);
cum_regret_rand= NaN(nreps,maxiter+1);
link = @normcdf;
modeltype = 'exp_prop';
theta = [1,-1]';

D = 1;
kernelname = 'ARD';
model.base_kernelfun = @ARD_kernelfun;
model.kernelname = kernelname;
model.modeltype = modeltype;
model.link = link;
model.regularization = 'nugget';
model.lb = zeros(D,1);
model.ub = ones(D,1);
model.lb_norm = zeros(D,1);
model.ub_norm = ones(D,1);

model.D = D;
for c = 1:2
    clear('xtrain', 'xtrain_norm', 'ctrain', 'score');
    condition = conditions(c);
    for r=1:nreplicates  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        seed  = seeds(r);
        objective = GPnd(D, theta, kernelname,seed);
        g = @(x) objective.do_eval(x);
        %             waitbar(((a-1)*nreplicates+r)/(nreplicates*nacq),wbar,'Computing...');
        [xtrain{r}, xtrain_norm{r}, ctrain{r}, score{r}] =  AL_preference_loop(acquisition_fun, seed, maxiter, theta, g, update_period, model,condition);
    end
    
    if c ==1
        score_wo = score;
    else
        score_w = score;
    end
end

scores{1} = cell2mat(score_wo');
scores{2} = cell2mat(score_w');
legends = {'Without conditionning', 'With conditionning'};

fig=figure('units','centimeters','outerposition',1+[0 0 fwidth fheight(mr)]);
fig.Color =  [1 1 1];
tiledlayout(mr, mc, 'TileSpacing', 'tight', 'padding','compact');
options.handle = fig;
options.alpha = 0.2;
options.error= 'sem';
options.line_width = linewidth/2;
options.semilogy = false;
options.cmap = C;
options.colors = colororder;
plots =  plot_areaerrorbar_grouped(scores, options);
legend(plots, legends, 'Fontsize', Fontsize);

