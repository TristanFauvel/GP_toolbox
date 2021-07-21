%% Gaussian process classification

clear all;
close all;
add_gp_module;

figure_path = '/home/tfauvel/Documents/PhD/Figures/Thesis_figures/Chapter_1/';
n=100;

rng(2)
graphics_style_paper;
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

colo= othercolor('GnBu7');
modeltype = 'exp_prop'; % or 'laplace'

link = @normcdf;

kernelfun = @ARD_kernelfun; 
theta_true= [3;3];

% Generate a function 
x = linspace(0,1,n);
y = mvnrnd(constant_mean(x,0), kernelfun(theta_true, x,x, 'false', 'false')); 
y=y-mean(y);

p= link(y);


figure();
subplot(1,2,1);
plot(x,y);
subplot(1,2,2);
plot(x,p);

ntr =10; % Number of data points

i_tr= randsample(n,ntr,'true');
xtrain = x(:,i_tr);
y_tr = y(:, i_tr);
ctrain = p(i_tr)>rand(1,ntr);

x_test = x;
y_test = y;

% GP classification with the correct hyperparameters
theta =theta_true ; % rand(size(theta_true));

[mu_c,  mu_y, sigma2_y, Sigma2_y, dmuc_dx, dmuy_dx, dsigma2y_dx, dSigma2y_dx, var_muc, dvar_muc_dx]= prediction_bin(theta, xtrain, ctrain, x_test, kernelfun, 'modeltype', modeltype);
fun = @(x_test) prediction_bin(theta, xtrain, ctrain, x_test, kernelfun, 'modeltype', modeltype);
dF = test_matrix_deriv(fun, x_test, 1e-8);

%Xlim= [min(x),max(x)];
%Ylim = [-5,5];

legend_pos = [-0.18,1];

mr = 1;
mc = 3;
fig=figure('units','centimeters','outerposition',1+[0 0 width height(mr)]);
fig.Color =  [1 1 1];
layout = tiledlayout(mr,mc, 'TileSpacing', 'tight', 'padding','compact');

i = 0;

nexttile();
i=i+1;
%errorshaded(x,mu_c, sqrt(var_muc), 'Color',  C(1,:),'LineWidth', linewidth, 'Fontsize', Fontsize); hold off
plot(x_test, p, 'LineWidth',linewidth,'Color',  C(2,:)) ; hold on;
plot(x_test, mu_c, 'LineWidth',linewidth,'Color', C(1,:)) ; hold on;
scatter(xtrain(ctrain == 1), ctrain(ctrain == 1), markersize, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k') ; hold on;
scatter(xtrain(ctrain == 0), ctrain(ctrain == 0), markersize, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k') ; hold off;
legend('$P(c=1)$', '$\mu_c(x)$','$c=1$', '$c=0$','Fontsize',Fontsize, 'Location', 'northeast')
xlabel('$x$','Fontsize',Fontsize)
legend boxoff
grid off
box off
pbaspect([1 1 1])
set(gca, 'Fontsize', Fontsize);
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)


nexttile();
i=i+1;
p1 = plot_gp(x,mu_y, sqrt(sigma2_y), C(1,:),linewidth);
p2 = plot(x,y,'LineWidth',linewidth,'Color', C(2,:)); hold off;
% errorshaded(x,mu_y, sqrt(sigma2_y), 'Color',  C(1,:),'LineWidth', linewidth, 'Fontsize', Fontsize); hold off
% legend('True function', 'Inferred function','Fontsize',Fontsize)
legend([p1, p2], {'$\mu_y(x)$','$g(x)$'},'Fontsize',Fontsize, 'Location', 'northeast')
xlabel('$x$','Fontsize',Fontsize)
legend boxoff
grid off
box off
pbaspect([1 1 1])
set(gca, 'Fontsize', Fontsize);
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)

nexttile();
i=i+1;
ax1 = imagesc(x, x, Sigma2_y); hold on;
xlabel('$x$','Fontsize',Fontsize)
ylabel('$x''$','Fontsize',Fontsize)
%title('Posterior covariance','Fontsize',Fontsize, 'interpreter', 'latex')
set(gca,'YDir','normal')
pbaspect([1 1 1])
cb = colorbar;
cb_lim = get(cb, 'Ylim');
cb_tick = get(cb, 'Ytick');
set(cb,'Ylim', cb_lim, 'Ytick', [cb_tick(1),cb_tick(end)]);
colormap(cmap)
set(gca, 'Fontsize', Fontsize);
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
title('Cov$(f(x),f(x'')|\mathcal{D})$')

figname  = 'GP_classification';
folder = [figure_path,figname];
savefig(fig, [folder,'/', figname, '.fig'])
exportgraphics(fig, [folder,'/' , figname, '.pdf']);
exportgraphics(fig, [folder,'/' , figname, '.png'], 'Resolution', 300);


nsamps =10000;
samples= mvnrnd(mu_y, Sigma2_y, nsamps);
mu_c_samples = normcdf(samples);
figure()
plot(x, mu_c_samples, 'color', C(1,:), 'linewidth', 1)




mr = 3;
mc = 3;
i = 0;
fig=figure('units','centimeters','outerposition',1+[0 0 width height(mr)]);
fig.Color =  [1 1 1];
tiledlayout(mr,mc, 'TileSpacing' , 'tight', 'Padding', 'tight')
nexttile()
i=i+1;
plot_gp(x,mu_y, sigma2_y, C(1,:), linewidth); 
xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$g(x)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
%title('Inferred value function $g(x)$','Fontsize', Fontsize)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off
% pbaspect([1 1 1])

nexttile()
i=i+1;
plot(x,mu_c, 'color', C(1,:)); 
xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$\mu_c(x)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
%title('Inferred value function $g(x)$','Fontsize', Fontsize)
text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
box off
% pbaspect([1 1 1])

nexttile()
i=i+1;
Hg = 0.5*log(2*pi*sigma2_y)+0.5;
plot(x, Hg, 'color', C(1,:))
xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$H(g(x)|\mathcal{D})$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
box off
title('Epistemic uncertainty')

h = @(p) -p.*log(p+eps) - (1-p).*log(1-p+eps);
nexttile()
i=i+1;
H = h(mu_c);
plot(x, H, 'color', C(1,:))
xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$H(c|\mu_c)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
box off
title('Aleatoric uncertainty')

nexttile()
i=i+1;
I= adaptive_sampling(theta, xtrain, ctrain, x, kernelfun, modeltype, []);
I = -I;
plot(x, I, 'color', C(1,:))
xlabel('$x$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
box off
ylabel('$I(c,g(x)|\mathcal{D})$')
title('Information gain')

nexttile()
i=i+1;
plot(x, var_muc, 'color', C(1,:))
xlabel('$x$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
box off
title('$V(\mu_c)$')


nexttile()
i=i+1;

Hc = mean(h(mu_c_samples),1);
plot(x,Hc, 'color', C(1,:))
title('$E(h(c|\mu_c))$') % 
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
box off

[a,b]= max(I)
new_x = x(b);
new_c = p(b)>rand;
[mu_c2,  mu_y2, sigma2_y2, Sigma2_y2, ~, ~,~,~, var_muc]= prediction_bin(theta, [xtrain,new_x], [ctrain, new_c], x_test, kernelfun, 'modeltype', modeltype);

mr = 1;
mc = 2;
i = 0;
fig=figure('units','centimeters','outerposition',1+[0 0 width height(mr)]);
fig.Color =  [1 1 1];
tiledlayout(mr,mc, 'TileSpacing' , 'tight', 'Padding', 'tight')
nexttile()
i=i+1;
H2 = h(mu_c2);
plot(x, H, 'color', C(1,:)); hold on;
plot(x, H2, 'color', C(2,:)); hold off;
xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$H(c|\mu_c)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
box off
title('Aleatoric uncertainty')

nexttile()

Hg2 = 0.5*log(2*pi*sigma2_y2)+0.5;
plot(x, Hg, 'color', C(1,:)); hold on;
plot(x, Hg2, 'color', C(2,:)); hold off;
xlabel('$x$', 'Fontsize', Fontsize)
ylabel('$H(c|\mu_c)$', 'Fontsize', Fontsize)
set(gca,'XTick',[0 0.5 1],'Fontsize', Fontsize)
ytick = get(gca,'YTick');
set(gca,'YTick', linspace(min(ytick), max(ytick), 3), 'Fontsize', Fontsize)
box off
title('Epistemic uncertainty')

% 
% % GP classification with the wrong hyperparameters
% theta = rand(size(theta_true));
% 
% [mu_c,  mu_y, sigma2_y]= prediction_bin(theta, xtrain, ctrain, x_test, kernelfun, 'modeltype', modeltype);
% 
% Ylim = [-5,5];
% h=figure(3);
% h.Color =  [1 1 1];
% h.Name = 'Activation function';
% plot(x,y,'LineWidth',1.5,'Color', C(1,:)); hold on;
% plot(x, mu_y,'LineWidth',1.5,'Color', colo(22,:)); hold off
% legend('True activation function', 'Inferred activation function')
% xlabel('x','Fontsize',Fontsize)
% ylabel('f(x)','Fontsize',Fontsize)
% set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',Ylim)
% grid off
% box off
% 
% h=figure(4);
% h.Color =  [1 1 1];
% h.Name = 'GP classification before hyperparameters optimization';
% scatter(xtrain, ctrain, 'MarkerFaceColor', C(1,:), 'MarkerEdgeColor', 'none') ; hold on;
% plot(x_test, mu_c, 'LineWidth',1.5,'Color', colo(22,:)) ; hold on;
% plot(x_test, p, 'LineWidth',1.5,'Color',  C(1,:)) ; hold off;
% legend('Data points', 'Inferred probability', 'True probability')
% xlabel('x','Fontsize',Fontsize)
% ylabel('f(x)','Fontsize',Fontsize)
% set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',[0,1])
% grid off
% box off
% 
% 
% %% Local optimization of hyperparameters
% options=[];
% theta = minFunc(@(hyp)negloglike_bin(hyp, xtrain, ctrain, kernelfun), theta, options);
% 
% %% Prediction with the new hyperparameters
% [mu_c,  mu_y, sigma2_y]= prediction_bin(theta, xtrain, ctrain, x_test, kernelfun);
% 
% h=figure(5);
% h.Color =  [1 1 1];
% h.Name = 'GP classification after hyperparameters optimization';
% scatter(xtrain, ctrain, 'MarkerFaceColor', C(1,:), 'MarkerEdgeColor', 'none') ; hold on;
% plot(x_test, mu_c, 'LineWidth',1.5,'Color', colo(22,:)) ; hold on;
% plot(x_test, p, 'LineWidth',1.5,'Color',  C(1,:)) ; hold off;
% legend('Data points', 'Inferred probability', 'True probability')
% xlabel('x','Fontsize',Fontsize)
% ylabel('f(x)','Fontsize',Fontsize)
% set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim,  'Ylim',[0,1])
% grid off
% box off
% 
% 
