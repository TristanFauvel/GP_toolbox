function [p, h] = plot_gp(x, mu_y, sigma2_y, col,linewidth, varargin)
opts = namevaluepairtostruct(struct( ...
    'background', [1 1 1] ...
    ), varargin);

UNPACK_STRUCT(opts, false)

num_quantiles= 10;
quantiles = linspace(0,0.5,num_quantiles+1);%0:0.05:0.5;
quantiles = quantiles(2:end);
i= 0;
for s = quantiles
    i=i+1;
    edges = [mu_y+norminv(s, 0, 1).*sqrt(sigma2_y); ...
        flipdim(mu_y-norminv(s, 0, 1).*sqrt(sigma2_y),1)];
    hc1{i} = fill([x'; flipdim(x',1)], edges, color_spectrum(2*s,col,background), 'EdgeColor', 'none'); hold on;
end
p = plot(x, mu_y, '-', 'color', col, 'Linewidth', linewidth); hold on;
alpha(1);
h = hc1{1};
end
function col = color_spectrum(p, col,background)
no_col = background(:);
full_col = col(:);
col = ((1 - p)*no_col + p*full_col)';
end

