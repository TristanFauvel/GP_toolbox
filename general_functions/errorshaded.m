function errorshaded(x, y, yerr, varargin)
% prettyplot(x, y, varargin)
%
% plot, with shaded errbars
%
% INPUTS
% x [1, nx] (default = 1:nx):   x coordinate
% y [1, nx]/[ntrials, nx]     : y coordinate
%
% OPTIONAL INPUTS
% errbars (default = 'std_err)
% Xlim
% LineStyle
% Fontsize
% Opactiy

opts = namevaluepairtostruct(struct(...
    'Color', 'black', ...
    'LineWidth', 1.5, ...
    'Opacity', 0.3, ...
    'Fontsize', 14, ...   
    'Xlim', [], ...
    'LineStyle', '-' ...
    ), varargin);


UNPACK_STRUCT(opts, false);


if ischar(Color)
    Color = rgb(Color);
end

% standard error on the mean, or standard deviation
y = y(:); yerr = yerr(:); x = x(:);
errmat = [y-yerr 2*yerr];          %matrix for errorbars


plot(x, y, 'Color', Color, 'LineWidth', LineWidth, 'LineStyle', LineStyle)   %main plot
hold on
g=area(x, errmat);                                  %error bars
set(g(1),'FaceColor', 'none', 'EdgeColor', 'none')  %make lower area-plot invisible
set(g(2),'FaceColor', Color, 'EdgeColor', 'none')   %colour for main error fill-in
alpha(Opacity);                                       %make transparent (alpha 1 makes opaque,0 is completely see-through)

if isempty(Xlim)
    Xlim = [x(1), x(end)];
end
set(gca, 'Fontsize', Fontsize, 'Xlim', Xlim)