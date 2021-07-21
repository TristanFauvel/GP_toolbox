Fontsize =44;
markersize = 50;
linewidth = 4;
set(0,'defaulttextInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',Fontsize)
set(0,'defaulttextFontSize',Fontsize)
set(0,'defaulttextlinewidth',linewidth)
colo= othercolor('GnBu7');
RdBu=cbrewer('div', 'RdBu', 255);
cmap = RdBu;
cmap(cmap<0) = 0;
