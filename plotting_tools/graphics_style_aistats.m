fact = 1;
colo= othercolor('GnBu7');
markersize = fact*10;
Fontsize = fact*9;
linewidth = fact*2;
set(0,'defaulttextInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',Fontsize)
set(0,'defaulttextFontSize',Fontsize)
set(0,'defaulttextlinewidth',linewidth)
%close all
colo= othercolor('GnBu7');
RdBu=cbrewer('div', 'RdBu', 255, 'spline');
cmap = flipud(RdBu);
cmap(cmap<0) = 0;
letter_font = Fontsize;

fwidth = 17.145;
fheight = [0.9/2*16, 1.5/2*16, 1.5/2*16, 16];

% RdBu=cbrewer('div', 'RdBu', 255, 'spline');
% cmap = flipud(RdBu);
% cmap(cmap<0) = 0;
mycolors =  [cmap(11,:);cmap(end-11,:)];
C = cbrewer('seq', 'YlGnBu', 255, 'spline');
mycolors = [mycolors; C(85,:)];
C = cbrewer('div', 'Spectral', 255, 'spline');
mycolors = [mycolors; C(end,:)];
C = mycolors;
colororder(C)
% colors_chart =  [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 1, 0, 0];
% colors_chart =  [cmap(1,:);cmap(end,:)];
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
