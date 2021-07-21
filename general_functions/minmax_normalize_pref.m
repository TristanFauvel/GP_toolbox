function [x_norm, min_x, max_x] = minmax_normalize_pref(x) 
%% Min-max normalization of the training set for preference learning with duels
d= size(x,1)/2;
min_x = NaN(2*d,1);
max_x = NaN(2*d,1);

for i =1:d
    xi = x(i:d:end,:);
    min_x([i,i+d],1) = min(xi(:));
    max_x([i,i+d],1) = max(xi(:));
end
x_norm = (x-min_x)./(max_x - min_x);


return

d= size(x,1)/2;
x1 = x(1:d:end,:);
x2 = x(2:d:end,:);
xflat = [x1(:),x2(:)]';
min_xd= min(xflat, [], 2);
max_xd= max(xflat, [], 2);
min_x = [min_xd; min_xd];
max_x = [max_xd; max_xd];
x_norm = (x-min_x)./(max_x - min_x);
