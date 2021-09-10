function [p1,p2, h] = plot_distro(x, mu_c, Y, col1, col2,linewidth)
num_quantiles= 10;
quantiles = linspace(0,0.5,num_quantiles+1);%0:0.05:0.5;
quantiles = quantiles(1:end);

% quantiles = linspace(0,0.5,num_quantiles+1);%0:0.05:0.5;
% quantiles = quantiles(2:end);
median = quantile(Y, 0.5);

for i = 1:num_quantiles 
    s= quantiles(i);
    %     du = norminv(s, 0, 1).*sqrt(sigma2_y);
%     dd = norminv(s, 0, 1).*sqrt(sigma2_y);
    du = quantile(Y, 0.5 + quantiles(end-i+1));
    dd = quantile(Y, 0.5 - quantiles(end-i+1));
%     edges = [median(:)+du(:); flipdim(median(:)-dd(:),1)];
    edges = [du(:); flipdim(dd(:),1)];

    hc1{i} = fill([x'; flipdim(x',1)], edges, color_spectrum(2*s,col2), 'EdgeColor', 'none'); hold on;
end
p2 = plot(x, median, '-', 'color', col2, 'Linewidth', linewidth); hold on;

p1 = plot(x, mu_c, '-', 'color', col1, 'Linewidth', linewidth); hold on;

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

