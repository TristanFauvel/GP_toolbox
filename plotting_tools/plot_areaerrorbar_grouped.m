
% ----------------------------------------------------------------------- %
% Function plot_areaerrorbar plots the mean and standard deviation of a   %
% set of data filling the space between the positive and negative mean    %
% error using a semi-transparent background, completely customizable.     %
%                                                                         %
%   Input parameters:                                                     %
%       - data:     Data matrix, with rows corresponding to observations  %
%                   and columns to samples.                               %
%       - options:  (Optional) Struct that contains the customized params.%
%           * options.handle:       Figure handle to plot the result.     %
%           * options.color_area:   RGB color of the filled area.         %
%           * options.color_line:   RGB color of the mean line.           %
%           * options.alpha:        Alpha value for transparency.         %
%           * options.line_width:   Mean line width.                      %
%           * options.x_axis:       X time vector.                        %
%           * options.error:        Type of error to plot (+/-).          %
%                   if 'std',       one standard deviation;               %
%                   if 'sem',       standard error mean;                  %
%                   if 'var',       one variance;                         %
%                   if 'c95',       95% confidence interval.              %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       data = repmat(sin(1:0.01:2*pi),100,1);                            %
%       data = data + randn(size(data));                                  %
%       plot_areaerrorbar(data);                                          %
% ----------------------------------------------------------------------- %
%   Author:  Victor Martinez-Cagigal                                      %
%   Date:    30/04/2018                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
function plots = plot_areaerrorbar(data, options)
% Default options
if(nargin<2)
    options.handle     = figure(1);
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'std';
    options.semilogy = false;
end
if(isfield(options,'x_axis')==0), options.x_axis = 1:size(data{1},2); end
options.x_axis = options.x_axis(:);

% Computing the mean and standard deviation of the data matrix
nplots = numel(data);
for i = 1:nplots
    data_mean{i} = mean(data{i},1);
    data_std{i}  = std(data{i},0,1);
    % Type of error plot
    switch(options.error)
        case 'std', error{i} = data_std{i};
        case 'sem', error{i} = (data_std{i}./sqrt(size(data{i},1)));
        case 'var', error{i} = (data_std{i}.^2);
        case 'c95', error{i} = (data_std{i}./sqrt(size(data{i},1))).*1.96;
    end

end
colors =  options.colors;
plots = [];
% Plotting the result
figure(options.handle);
x_vector = [options.x_axis', fliplr(options.x_axis')];
if options.semilogy
    error = cell2mat(error');
    data_mean = cell2mat(data_mean');

    h =  semilogy(options.x_axis, data_mean', ...
        'LineWidth', options.line_width);
    hold on;
    set(h, {'color'}, num2cell(colors(1:nplots,:),2));
    for i = 1:nplots
        patch = fill(x_vector, [(data_mean(i,:) + error(i,:)),fliplr(data_mean(i,:)-error(i,:))], colors(i,:));
        set(patch, 'edgecolor', 'none');
        set(patch, 'facecolor',colors(i,:));

        set(patch, 'FaceAlpha', options.alpha);
        hold on;
    end
    plots = h;
else
    for i = 1:nplots
        if nplots>size(colors,1)
            h =  plot(options.x_axis, data_mean{i}, ...
                'LineWidth', options.line_width, 'LineStyle', options.lines{i}); hold on;
            plots = [plots, h];
            patch = fill(x_vector, [data_mean{i}+error{i},fliplr(data_mean{i}-error{i})], h.Color);
            set(patch, 'edgecolor', 'none');
            set(patch, 'FaceAlpha', options.alpha);
            hold on;

        else
            patch = fill(x_vector, [data_mean{i}+error{i},fliplr(data_mean{i}-error{i})], colors(i,:));
            set(patch, 'edgecolor', 'none');
            set(patch, 'FaceAlpha', options.alpha);
            hold on;
            h =  plot(options.x_axis, data_mean{i}, 'color',colors(i,:), ...
                'LineWidth', options.line_width, 'LineStyle', options.lines{i});
            plots = [plots, h];

        end
    end
    set(gca, 'xlim', options.xlim)
end
