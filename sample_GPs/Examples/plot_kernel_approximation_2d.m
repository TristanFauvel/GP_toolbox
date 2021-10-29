letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
graphics_style_paper;
figure_path = '/home/tfauvel/Documents/PhD/Figures/Thesis_figures/Chapter_1/';


approximation.method = 'RRGP'; %'RRGP'

kernelname = 'Gaussian';

if strcmp(kernelname,'ARD')
    kernelfun = @ARD_kernelfun;
    theta.cov = [-2*log(0.1),-2*log(0.1), 0];
    lengthscale = exp(-theta.cov(1)/2);
elseif strcmp(kernelname,'Gaussian')
    kernelfun = @Gaussian_kernelfun;
    theta.cov = [-2*log(0.1), 0];
    lengthscale = exp(-theta.cov(1)/2);
elseif strcmp(kernelname,'Matern52')
    kernelfun = @Matern52_kernelfun;
    lengthscale = 1/10;
    theta.cov = [log(lengthscale),0];
elseif strcmp(kernelname,'Gaussian')
    kernelfun = @Matern32_kernelfun;
    lengthscale = 1/10;
    theta.cov = [log(lengthscale),0];
end

theta.mean = 0;

D = 2;
n= 100;
lb = 0;
ub = 1;
xrange = linspace(lb, ub, n);
[xx,yy]= ndgrid(xrange,xrange);
xx= xx(:)';
yy = yy(:)';
x = [xx;yy];
x0 = [0.5;0.5];
meanfun = @constant_mean;
K = kernelfun(theta.cov,x0,x, 'false', 'nugget');

regularization = 'nugget';
type = 'regression';
hyps.ncov_hyp =numel(theta.cov); % number of hyperparameters for the covariance function
hyps.nmean_hyp =1; % number of hyperparameters for the mean function
hyps.hyp_lb = -10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
hyps.hyp_ub = 10*ones(hyps.ncov_hyp  + hyps.nmean_hyp,1);
model = gp_regression_model(D, meanfun, kernelfun, regularization, hyps, lb, ub, kernelname);



m = 64;
approximation.nfeatures = m;
phi = sample_features_GP(theta, model, approximation);
K1 = phi(x0)*phi(x)';

m = 128;
approximation.nfeatures = m;
phi = sample_features_GP(theta, model, approximation);
K2 = phi(x0)*phi(x)';

m = 256;
approximation.nfeatures = m;
phi = sample_features_GP(theta, model, approximation);
K3 = phi(x0)*phi(x)';

m = 512;
approximation.nfeatures = m;
phi = sample_features_GP(theta, model, approximation);
K4 = phi(x0)*phi(x)';

m = 1024;
approximation.nfeatures = m;
phi = sample_features_GP(theta, model, approximation);
K5 =phi(x0)*phi(x)';

mr = 1;
mc = 5;
legend_pos = [0.02,1];

clim= [0,1];

fig=figure('units','centimeters','outerposition',1+[0 0 fwidth fheight(mr)]);
fig.Color =  background_color;
layout = tiledlayout(mr,mc, 'TileSpacing', 'tight', 'padding','compact');
i = 0;
ticks = [-0.5,0.5]./lengthscale;
% ylim = [0,1];
nexttile();
i=i+1;
imagesc(xrange, xrange, reshape(K2,n,n), clim);
box off
xlabel('$x$')
ylim = get(gca, 'YLim');
set(gca,'YLim',ylim,'YTick',[0,0.5,1], 'Xtick',[0,0.5, 1])
set(gca,'xticklabel',{[num2str(ticks(1)),'$\lambda$'], '0', [num2str(ticks(2)),'$\lambda$']}, 'Fontsize', Fontsize)
% text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
title('$m = 128$', 'Fontsize', Fontsize)
pbaspect([1 1 1])

nexttile();
i=i+1;
imagesc(xrange, xrange,reshape(K3,n,n), clim);
box off
xlabel('$x$')
set(gca,'YLim',ylim,'YTick',[0,0.5,1], 'Xtick',[0,0.5, 1])
set(gca,'xticklabel',{[num2str(ticks(1)),'$\lambda$'], '0', [num2str(ticks(2)),'$\lambda$']}, 'Fontsize', Fontsize)
% text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
title('$m =256$', 'Fontsize', Fontsize)
pbaspect([1 1 1])

nexttile();
i=i+1;
imagesc(xrange, xrange,reshape(K4,n,n), clim);
box off
xlabel('$x$')
set(gca,'YLim',ylim,'YTick',[0,0.5,1], 'Xtick',[0,0.5, 1])
set(gca,'xticklabel',{[num2str(ticks(1)),'$\lambda$'], '0', [num2str(ticks(2)),'$\lambda$']}, 'Fontsize', Fontsize)
% text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
title('$m = 512$', 'Fontsize', Fontsize)
pbaspect([1 1 1])

nexttile();
i=i+1;
imagesc(xrange, xrange,reshape(K5,n,n), clim);
box off
xlabel('$x$')
set(gca,'YLim',ylim,'YTick',[0,0.5,1], 'Xtick',[0,0.5, 1])
set(gca,'xticklabel',{[num2str(ticks(1)),'$\lambda$'], '0', [num2str(ticks(2)),'$\lambda$']}, 'Fontsize', Fontsize)
% text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
title('$m = 1024$', 'Fontsize', Fontsize)
pbaspect([1 1 1])

nexttile();
i=i+1;
imagesc(xrange, xrange,reshape(K,n,n), clim);
box off
xlabel('$x$')
set(gca,'YLim',ylim,'YTick',[0,0.5,1], 'Xtick',[0,0.5, 1])
set(gca,'xticklabel',{[num2str(ticks(1)),'$\lambda$'], '0', [num2str(ticks(2)),'$\lambda$']}, 'Fontsize', Fontsize)
% text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
title('True', 'Fontsize', Fontsize)
pbaspect([1 1 1])


figname  = 'approximation';
folder = [figure_path,figname];
figname  = 'approximation_2D';

savefig(fig, [folder,'/', figname, '.fig'])
exportgraphics(fig, [folder,'/' , figname, '.pdf']);
exportgraphics(fig, [folder,'/' , figname, '.png'], 'Resolution', 300);

