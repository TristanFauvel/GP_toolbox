function [p, h] = plot_gp(x, mu_y, sigma2_y, col,linewidth)
num_quantiles= 10;
quantiles = linspace(0,0.5,num_quantiles+1);%0:0.05:0.5;
quantiles = quantiles(2:end);
i= 0;
for s = quantiles
    i=i+1;
    edges = [mu_y+norminv(s, 0, 1).*sqrt(sigma2_y); ...
        flipdim(mu_y-norminv(s, 0, 1).*sqrt(sigma2_y),1)];
    hc1{i} = fill([x'; flipdim(x',1)], edges, color_spectrum(2*s,col), 'EdgeColor', 'none'); hold on;
end
p = plot(x, mu_y, '-', 'color', col, 'Linewidth', linewidth); hold on;
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

