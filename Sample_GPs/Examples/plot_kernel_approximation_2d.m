letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
graphics_style_paper;
figure_path = '/home/tfauvel/Documents/PhD/Figures/Thesis_figures/Chapter_1/';

% kernelfun = @ARD_kernelfun;
% kernelname = 'ARD';
% method = 'RRGP'; %'RRGP'
% theta = [-2*log(0.1),0];
% lengthscale = exp(-theta(1)/2);
kernelfun = @Matern52_kernelfun;
kernelname = 'Matern52';
method = 'RRGP'; %'RRGP'
lengthscale = 1/10;
theta = [log(lengthscale),0];

D = 2;
n= 100;
xrange = linspace(0,1,n);
[xx,yy]= ndgrid(xrange,xrange);
xx= xx(:)';
yy = yy(:)';
x = [xx;yy];
x0 = [0.5;0.5];

K = kernelfun(theta,x0,x);

m = 4;
phi = sample_features_GP(theta, D, model, approximation);
K1 = phi(x0)*phi(x)';

m = 8;
phi = sample_features_GP(theta, D, model, approximation);
K2 = phi(x0)*phi(x)';

m = 16;
phi = sample_features_GP(theta, D, model, approximation);
K3 = phi(x0)*phi(x)';

m = 32;
phi = sample_features_GP(theta, D, model, approximation);
K4 = phi(x0)*phi(x)';

m = 64;
phi = sample_features_GP(theta, D, model, approximation);
K5 = phi(x0)*phi(x)';

mr = 1;
mc = 5;
legend_pos = [0.02,1];

fig=figure('units','centimeters','outerposition',1+[0 0 width height(mr)]);
fig.Color =  [1 1 1];
layout = tiledlayout(mr,mc, 'TileSpacing', 'tight', 'padding','compact');
i = 0;
ticks = [-0.5,0.5]./lengthscale;
% ylim = [0,1];
nexttile();
i=i+1;
imagesc(reshape(K2,n,n));
box off
xlabel('$x$')
ylim = get(gca, 'YLim');
set(gca,'YLim',ylim,'YTick',[0,0.5,1], 'Xtick',[0,0.5, 1])
set(gca,'xticklabel',{[num2str(ticks(1)),'$\lambda$'], '0', [num2str(ticks(2)),'$\lambda$']}, 'Fontsize', Fontsize)
% text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
title('$m = 8$', 'Fontsize', Fontsize)
pbaspect([1 1 1])

nexttile();
i=i+1;
imagesc(reshape(K3,n,n));
box off
xlabel('$x$')
set(gca,'YLim',ylim,'YTick',[0,0.5,1], 'Xtick',[0,0.5, 1])
set(gca,'xticklabel',{[num2str(ticks(1)),'$\lambda$'], '0', [num2str(ticks(2)),'$\lambda$']}, 'Fontsize', Fontsize)
% text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
title('$m = 16$', 'Fontsize', Fontsize)
pbaspect([1 1 1])

nexttile();
i=i+1;
imagesc(reshape(K4,n,n));
box off
xlabel('$x$')
set(gca,'YLim',ylim,'YTick',[0,0.5,1], 'Xtick',[0,0.5, 1])
set(gca,'xticklabel',{[num2str(ticks(1)),'$\lambda$'], '0', [num2str(ticks(2)),'$\lambda$']}, 'Fontsize', Fontsize)
% text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
title('$m = 32$', 'Fontsize', Fontsize)
pbaspect([1 1 1])

nexttile();
i=i+1;
imagesc(reshape(K5,n,n));
box off
xlabel('$x$')
set(gca,'YLim',ylim,'YTick',[0,0.5,1], 'Xtick',[0,0.5, 1])
set(gca,'xticklabel',{[num2str(ticks(1)),'$\lambda$'], '0', [num2str(ticks(2)),'$\lambda$']}, 'Fontsize', Fontsize)
% text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
title('$m = 64$', 'Fontsize', Fontsize)
pbaspect([1 1 1])

nexttile();
i=i+1;
imagesc(reshape(K,n,n));
box off
xlabel('$x$')
set(gca,'YLim',ylim,'YTick',[0,0.5,1], 'Xtick',[0,0.5, 1])
set(gca,'xticklabel',{[num2str(ticks(1)),'$\lambda$'], '0', [num2str(ticks(2)),'$\lambda$']}, 'Fontsize', Fontsize)
% text(legend_pos(1), legend_pos(2),['$\bf{', letters(i), '}$'],'Units','normalized','Fontsize', letter_font)
title('$m = 64$', 'Fontsize', Fontsize)
pbaspect([1 1 1])


figname  = 'approximationimation';
folder = [figure_path,figname];
figname  = 'approximationimation_2D';

savefig(fig, [folder,'/', figname, '.fig'])
exportgraphics(fig, [folder,'/' , figname, '.pdf']);
exportgraphics(fig, [folder,'/' , figname, '.png'], 'Resolution', 300);

