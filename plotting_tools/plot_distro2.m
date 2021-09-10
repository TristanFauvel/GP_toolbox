function [p1,p2, h] = plot_distro2(x, mu_c, Y, col1, col2, groups)
num_quantiles= 10;
quantiles = linspace(0,0.5,num_quantiles+1);%0:0.05:0.5;
quantiles = quantiles(1:end);

median = quantile(Y, 0.5);
x= x(:);
xf = flipdim(x,1);

for i = 1:num_quantiles
    s= quantiles(i);
    du = quantile(Y, 0.5 + quantiles(end-i+1));
    du = du(:);
    dd = quantile(Y, 0.5 - quantiles(end-i+1));
    dd = flipdim(dd(:),1);
    
    for j = numel(groups)
    xx = [x(groups{j}); xf(groups{j})];
    edges = [du(groups{j}); dd(groups{j})];
    hc1{i} = fill(xx, edges, color_spectrum(2*s,col2), 'EdgeColor', 'none'); hold on;
    end
end


p2 = scatter(x, median, 10, col2, 'filled'); hold on;
p1 = scatter(x, mu_c, 10, col1, 'filled'); hold on;

% plot(x, mu_y+sqrt(sigma2_y), 'color', col, 'Linewidth', linewidth); hold on;
% plot(x, mu_y-sqrt(sigma2_y), 'color', col, 'Linewidth', linewidth); hold on;
alpha(1);
h = hc1{1};
end
function col = color_spectrum(p, col)
no_col = [1 1 1];
full_col = col;
col = (1 - p)*no_col + p*full_col;
end

