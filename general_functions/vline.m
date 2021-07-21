function line = vline(x, varargin)
% plots vertical line, preserving axis
%
% vline(x, varargin)
%
% OPTIONAL INPUTS:
% Color
% axis_handle
% LineStyle
% Linewidth

 
opts = namevaluepairtostruct(struct( ...    
    'Color', 'k', ...
    'axis_handle', gca,...    
    'LineStyle', '--', ...
    'Linewidth', 1, ...
    'ymax', [] ...
     ), varargin);

 
UNPACK_STRUCT(opts, false)



% 
% line = plot(axis_handle, [x x], max(abs(ymax), 1e3)*[-1 1]*1e3, ...
%     LineStyle,...
%     'Color', Color, ...
%     'Linewidth', Linewidth); hold on

 
Xlim = get(axis_handle, 'Xlim');
Ylim = get(axis_handle, 'Ylim');
if isempty(ymax)
    ymax= max(abs(Ylim));
end
 
%line = plot(axis_handle, [x x], [0,ymax], LineStyle, 'Color', Color,'Linewidth', Linewidth); hold on
line = plot(axis_handle, [x x], [Ylim(1),ymax], LineStyle, 'Color', Color,'Linewidth', Linewidth); hold on

%  line = plot(axis_handle, [x x], max(max(abs(Ylim)), 1e3)*[-1 1]*1e3, ...
%     LineStyle,...
%     'Color', Color, ...
%     'Linewidth', Linewidth); hold on

set(axis_handle, 'Xlim', Xlim, 'Ylim', Ylim)