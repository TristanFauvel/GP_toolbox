
clear all
close all

graphics_style_paper;

link = @normcdf; %inverse link function

%% Define the range of parameters
n = 1000;
x = linspace(0,1, n);

[p,q]= meshgrid(x);
x2d_test = [p(:), q(:)]';
colo= othercolor('GnBu7');


gfunc = @forretal08;
g = gfunc(x);

%% Find the true global optimum of g
[gmax, id_xmax] = max(g);
xmax = x(id_xmax);

fig=figure();
fig.Color =  [1 1 1];
fig.Name = 'Value function $g(x)$';



ax1 = subplot(1,3,1);
plot(x, g, 'Color', colo(end, :),'LineWidth',1.5); hold on;
v=vline(xmax, 'Color', 'k','LineWidth',1.5); hold off;
xlabel('$x$','Fontsize',Fontsize)
ylabel('$g(x)$','Fontsize',Fontsize)
title('Value function $g(x)$','Fontsize',Fontsize)
pbaspect(ax1, [1 1 1])
box off

%% Plot the corresponding preference map
% fig=figure(2);
% fig.Color =  [1 1 1];
ax2 =subplot(1,3,2);
imagesc(x, x, (g'-g))
xlabel('$x$','Fontsize',Fontsize)
ylabel('$x''$','Fontsize',Fontsize)
title('Preference function $f(x,x'')$','Fontsize',Fontsize)
set(gca,'YDir','normal')
pbaspect(ax2, [1 1 1])
colorbar
colormap(cmap)
ax = get(gca);

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% fig=figure(3);
% fig.Color =  [1 1 1];
ax3 =subplot(1,3,3);
imagesc(x, x, link(g'-g))
xlabel('$x$','Fontsize',Fontsize)
ylabel('$x''$','Fontsize',Fontsize)
set(gca,'YDir','normal')
title('$\Phi(f(x, x''))$','Fontsize',Fontsize)

pbaspect(ax3,[1 1 1])
colorbar
colormap(cmap)
ax = get(gca);

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% saveTightFigure(fig,'/home/tfauvel/Desktop/optim_retina/Figures/preference_plot.png')
saveas(fig, 'preference_plot','png')